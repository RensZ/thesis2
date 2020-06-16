
#include "MyApplications/customFunctions.h"



// convert calendar date and time to seconds after J2000
double secondsAfterJ2000(Eigen::Vector6i datetime){
    using namespace tudat;
    using namespace tudat::basic_astrodynamics;
    const unsigned int secondsPerDay = 60*60*24;
    const unsigned int julianDayJ2000 = 2451545;
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay(datetime[0],datetime[1],datetime[2],datetime[3],datetime[4],datetime[5]);
    return (julianDay - julianDayJ2000)*secondsPerDay;
}

double secondsAfterJ2000(std::vector<int> datetime){
    if (datetime.size() != 6){
        std::runtime_error("error, vector input of function secondsAfterJ2000 should have size 6 (Y,M,D,h,m,s)");
    }
    using namespace tudat;
    using namespace tudat::basic_astrodynamics;
    const unsigned int secondsPerDay = 60*60*24;
    const unsigned int julianDayJ2000 = 2451545;
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay(datetime.at(0),datetime.at(1),datetime.at(2),datetime.at(3),datetime.at(4),datetime.at(5));
    return (julianDay - julianDayJ2000)*secondsPerDay;
}


//! Function to calculate the Sun-Mercury-Earth angle
//! this is a simplified approach which simply uses the mean motion
//! w.r.t. an epoch where both planets alligned
//double mercurySunEarthAngle(const double time){

//    // from NASA fact sheet https://nssdc.gsfc.nasa.gov/planetary/factsheet/
//    const double yearMercury = 88.0  *24.0*60.0*60.0;
//    const double yearEarth   = 365.2 *24.0*60.0*60.0;

//    // Get mean motion of planets
//    const double pi = 3.14159265359;
//    const double meanMotionMercury = 2*pi/yearMercury;
//    const double meanMotionEarth   = 2*pi/yearEarth;

//    // get time since last eclipse, when Mercury and Earth perfectly aligned (MSE = 0)
//    Eigen::Vector6i timeAlignedVector;
//    timeAlignedVector << 2019, 11, 11, 15, 20, 0; // YYYY, MM, DD, hh, mm, ss
//    const double timeAligned = secondsAfterJ2000(timeAlignedVector);
//    const double timeDelta = time-timeAligned;

//    // get angular distance traveled since last eclipse (or in the past)
//    const double angleTravelledMercury = timeDelta*meanMotionMercury;
//    const double angleTravelledEarth   = timeDelta*meanMotionEarth;

//    // get value between 0-180 degrees
//    double relativeAngle = std::fmod((abs(angleTravelledMercury-angleTravelledEarth)),2*pi);
//    if (relativeAngle > pi){
//        relativeAngle = 2*pi-relativeAngle;
//    }
//    return relativeAngle;
//}

//! Function to calculate the angle between 2 bodies w.r.t. a central body
//! Based on spice data
double angleBetween2Bodies(const double time,
                           const std::string centralBody,
                           const std::string body1,
                           const std::string body2)
{
    using namespace tudat::spice_interface;

    Eigen::Vector3d posBody1 = getBodyCartesianPositionAtEpoch(
                body1,centralBody,"ECLIPJ2000","None",time);
    Eigen::Vector3d posBody2 = getBodyCartesianPositionAtEpoch(
                body2,centralBody,"ECLIPJ2000","None",time);

    double cosAngle = posBody1.dot(posBody2)/
            ( posBody1.norm() * posBody2.norm() );

    return acos(cosAngle);
}


//! Function to generate a list of observation times between an initial and final time, with a time step in between.
//! a maximum number of observations can be indicated, if the pattern generates more than that, a random selection will be deleted to meet the maximum.
//! a maximum Mercury-Sun-Earth angle can be indicated, observations will not be included above the requirement (see function above)
//! in addition times can be added to the list manually which will be added last
std::vector<double> makeObservationTimeList(const double initialTime,
                                            const double endTime,
                                            const double timeStep,
                                            const double arcDuration,
                                            const unsigned int maximumNumberOfTrackingDays, //if not applicable, set to a very large number.
                                            const double maxMSEangle, //radians. if not applicable, set to value greater than pi.
                                            const std::vector<double> flybyTimes)
{
    std::vector< std::vector<double> > dailyObservationListCollection;
    const double day = 86400.0; //assumed frequency of tracking arcs once per day

    // add times to observation list until the duration of a tracking arc has passed
    double currentTime = initialTime;
    std::vector<double> currentDayObservationList;
    while (currentTime <= endTime){

        // add data points that satisfy the MSE angle requirement
        double trackingTime = currentTime;
        while (trackingTime <= currentTime+arcDuration){
            double MSEangle = angleBetween2Bodies(currentTime, "Sun", "Mercury", "Earth");
            if (MSEangle < maxMSEangle){
                currentDayObservationList.push_back(trackingTime);
            }
            trackingTime += timeStep;
        }

        // save data if more than 90% of the arc satisfied the MSE requirement
        if (currentDayObservationList.size() > (arcDuration/timeStep)*0.90){
            dailyObservationListCollection.push_back(currentDayObservationList);
        }

        // move on to next day
        currentDayObservationList.clear();
        currentTime += day;
    }

    // if list is longer than the max number of days, delete random entries
    srand(0);
    while ( dailyObservationListCollection.size() > maximumNumberOfTrackingDays ){
        unsigned int numberOfObservations = dailyObservationListCollection.size();
        unsigned int randomIndex = rand()%(numberOfObservations-1);
        dailyObservationListCollection.erase( dailyObservationListCollection.begin() + randomIndex );
    }
    std::cout<<"total number of tracking days: "<<dailyObservationListCollection.size()<<std::endl;


    // append all entries of the collection to one big vector for output
    std::vector< double > finalObservationList;
    for( unsigned int i = 0; i < dailyObservationListCollection.size(); i++ ){
        currentDayObservationList = dailyObservationListCollection.at(i);
        for( unsigned int j = 0; j < currentDayObservationList.size(); j++ ){
            finalObservationList.push_back(currentDayObservationList.at(j));
        }
    }

    // add times which are manually indicated (flybys)
    for( unsigned int i = 0; i < flybyTimes.size(); i++ )
        finalObservationList.push_back(flybyTimes.at(i));

    std::cout<<"total number of observation epochs: "<<finalObservationList.size()<<std::endl;
    return finalObservationList;
}


double noiseLevelBasedOnMSEangle(const double time,
                                 const double noiseAtMinAngle,
                                 const double noiseAtMaxAngle,
                                 const double maxAngleDeg){ //in degrees

    double noise;

    if (noiseAtMaxAngle == noiseAtMinAngle){
        noise = noiseAtMaxAngle;
    } else{

        using namespace tudat::unit_conversions;

        // calculate MSE angle with previous function
        double relativeAngle = angleBetween2Bodies(time, "Sun", "Mercury", "Earth");

        // get noise level (linear function between a minimum and maximum point)
        const double maxAngle = convertDegreesToRadians(maxAngleDeg);
        const double minAngle = 0;

        const double slope = (noiseAtMaxAngle-noiseAtMinAngle)/(maxAngle-minAngle);
        const double intercept = noiseAtMaxAngle - slope*maxAngle;
        noise = slope*relativeAngle + intercept;

    }

    return noise;
}


//! Function to make the noise level dependent on the Mercury-Sun-Earth angle,
//! based on a linear relation constructed with the minimum and maximum noise.
double noiseSampleBasedOnMSEangle(const double time,
                            const double noiseAtMinAngle,
                            const double noiseAtMaxAngle,
                            const double maxAngleDeg){ //in degrees

    double noise = noiseLevelBasedOnMSEangle(time, noiseAtMinAngle, noiseAtMaxAngle, maxAngleDeg);

    // create a gaussian sample
    boost::random::mt19937 rng(time);
    boost::random::normal_distribution<> nd(0.0,noise);
    const double sample = nd(rng);

//    std::cout << "SPE (deg): " << relativeAngle*180/pi << " ---- "
//              << "noise level (m): " << noise << " ---- "
//              << "rng noise sample (m): " << sample << std::endl;

    return sample;
}

//! Function to make the noise level dependent on the Mercury-Sun-Earth angle,
//! based on a linear relation constructed with the minimum and maximum noise.
double noiseSampleBasedOnMSEangleForMultipleMissions(const double time,
                            std::vector< double > noiseAtMinAngleVector,
                            std::vector< double > noiseAtMaxAngleVector,
                            std::vector< double > maxAngleDegVector, //in degrees
                            std::vector< std::vector< double > > seperateBaseTimeLists){

    double noise;
    unsigned int check = 0;

    for (unsigned int m = 0; m<seperateBaseTimeLists.size(); m++){
        std::vector< double > currentBaseTimeList = seperateBaseTimeLists.at(m);

        if (std::find(currentBaseTimeList.begin(), currentBaseTimeList.end(), time)
                != currentBaseTimeList.end()) {

            noise = noiseLevelBasedOnMSEangle(
                        time, noiseAtMinAngleVector.at(m),
                        noiseAtMaxAngleVector.at(m), maxAngleDegVector.at(m));
            check += 1;

//            std::cout<<"m: "<<m<<" noise: "<<noise;
        }
    }

    if ( check == 0){
        std::cout<<"finding noise value failed for t = "<<time<<std::endl;
        std::runtime_error("time not found in any mission observation list");
    }
    if ( check > 1){
        std::cout<<"finding noise value failed for t = "<<time
                 <<" check = "<<check<<std::endl;
        std::runtime_error("time found in multiple mission observation lists");
    }

    // create a gaussian sample
    boost::random::mt19937 rng(time);
    boost::random::normal_distribution<> nd(0.0,noise);
    const double sample = nd(rng);

//    std::cout << "SPE (deg): " << relativeAngle*180/pi << " ---- "
//              << "noise level (m): " << noise << " ---- "
//              << "rng noise sample (m): " << sample << std::endl;

//    std::cout<<" sample: "<<sample<<std::endl;
    return sample;
}



double averageOfDoubleVector(std::vector<double> input){
    double average = 0.0;
    double vectorsize = input.size( );
    for(unsigned int i=0; i<vectorsize; i++){
        average += input.at( i )/vectorsize;
    }
    return average;
}


Eigen::MatrixXd interpolatePositionErrors(const Eigen::MatrixXd errorMatrix,
                                          const std::vector<double> timesAtWhichToInterpolate){

    std::vector<double> t(errorMatrix.col(0).data(), errorMatrix.col(0).data() + errorMatrix.rows());
    std::vector<double> x(errorMatrix.col(1).data(), errorMatrix.col(1).data() + errorMatrix.rows());
    std::vector<double> y(errorMatrix.col(2).data(), errorMatrix.col(2).data() + errorMatrix.rows());
    std::vector<double> z(errorMatrix.col(3).data(), errorMatrix.col(3).data() + errorMatrix.rows());

    std::cout<<x.at(0)<<" "<<x.at(t.size()-1)<<" "<<x.size();
    std::cout<<y.at(0)<<" "<<y.at(t.size()-1)<<" "<<y.size();
    std::cout<<z.at(0)<<" "<<z.at(t.size()-1)<<" "<<z.size();

    using namespace tudat::interpolators;
    LinearInterpolatorDouble xInterpolator(t, x, huntingAlgorithm, use_boundary_value);
    LinearInterpolatorDouble yInterpolator(t, y, huntingAlgorithm, use_boundary_value);
    LinearInterpolatorDouble zInterpolator(t, z, huntingAlgorithm, use_boundary_value);

    Eigen::MatrixXd interpolatedErrorMatrix(timesAtWhichToInterpolate.size(),3);
    for( unsigned int i=0; i<timesAtWhichToInterpolate.size(); i++){
        double currentTime = timesAtWhichToInterpolate.at(i);

        interpolatedErrorMatrix(i,0) = xInterpolator.interpolate(currentTime);
        interpolatedErrorMatrix(i,1) = yInterpolator.interpolate(currentTime);
        interpolatedErrorMatrix(i,2) = zInterpolator.interpolate(currentTime);

//        std::cout<<"t: " <<currentTime
//                 <<" x: "<<interpolatedErrorMatrix(i,0)
//                 <<" y: "<<interpolatedErrorMatrix(i,1)
//                 <<" y: "<<interpolatedErrorMatrix(i,2)<<std::endl;

    }

    return interpolatedErrorMatrix;
}





Eigen::MatrixXd interpolatePositionErrorsBasedOnTrueAnomaly(const Eigen::MatrixXd errorMatrix,
                                                            const std::vector<double> timesAtWhichToInterpolate,
                                                            const std::string vehicle,
                                                            const double mercuryGravitationalParameter){

    using namespace tudat::spice_interface;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::input_output;
    using namespace tudat::interpolators;
    using namespace tudat::unit_conversions;

    // load SPICE data of vehicles
    std::string vehicleName; std::string vehicleKernel;
    Eigen::Vector6i minKernelDatetime;
    if (vehicle == "BepiColombo"){
        vehicleName = "BEPICOLOMBO MPO";
        vehicleKernel = "bc_mpo_mlt_50037_20260314_20280529_v01.bsp";
        minKernelDatetime << 2026, 3, 14, 0, 0, 0;
        loadSpiceKernelInTudat(getSpiceKernelPath() + "bc_sci_v04.tf");
    }
    else if (vehicle == "MESSENGER"){
        vehicleName = "MESSENGER";
        vehicleKernel = "msgr_040803_150430_150430_od431sc_2.bsp";
        minKernelDatetime << 2011, 3, 17, 0, 0, 0;
    } else{
        std::runtime_error("ERROR: vehicle name not recognised");
    }

    loadSpiceKernelInTudat(getSpiceKernelPath() + vehicleKernel);
    double minKernelTime = secondsAfterJ2000(minKernelDatetime);

    // load True Anomaly data (output of python)
    std::string vehicleArcIndicesFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/arcindices_"+vehicle+".txt";
    std::string vehicleTrueAnomalyFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/trueanomaly_inputs_"+vehicle+".txt";
    Eigen::VectorXd arcIndices = readMatrixFromFile(vehicleArcIndicesFilename, ",");
    Eigen::MatrixXd trueAnomalyInput = readMatrixFromFile(vehicleTrueAnomalyFilename, ",");

    // iterate over times at which an error value is desired
    Eigen::MatrixXd interpolatedErrorMatrix(timesAtWhichToInterpolate.size(),3);
    for( unsigned int i=0; i<timesAtWhichToInterpolate.size(); i++){

        double currentTime = timesAtWhichToInterpolate.at(i);

        // find indices of previous and next arc
        bool inBetweenTwoArcsFlag = false;
        Eigen::VectorXd arcAverageTimes = errorMatrix.col(0);
        unsigned int numberOfArcs = arcAverageTimes.size();
        int previousArc = -1;

        for (unsigned int i=0; i<numberOfArcs-1; i++){
            if (arcAverageTimes(i) <= currentTime
                    && arcAverageTimes(i+1) >= currentTime ){
                previousArc = i;
                inBetweenTwoArcsFlag = true;
            }
        }

        // if not found, check if the asked time is before the first arc or after the last
        if (previousArc == -1){
            if (arcAverageTimes(0) >= currentTime){
                previousArc = 0;
            } else if (arcAverageTimes(numberOfArcs) <= currentTime){
                previousArc = numberOfArcs-1;
            } else{
                std::runtime_error("ERROR: the previous and next arc could not be found");
            }
        }

        // find true anomaly of vehicle at current time
        double vehicleTrueAnomaly;
        if (currentTime < minKernelTime){ // if flyby, there is no orbit/true anomaly to speak of. Error at pericenter is taken of closest arc.
            vehicleTrueAnomaly = 0.0;
        } else{
            Eigen::Vector6d vehicleState = getBodyCartesianStateAtEpoch(vehicleName,"Mercury","ECLIPJ2000","None",currentTime);
            vehicleTrueAnomaly = convertRadiansToDegrees(convertCartesianToKeplerianElements(vehicleState, mercuryGravitationalParameter)(5));
        }

        // if before the first arc or after the last arc,
        // find error corresponding to the given TA of the nearest arc
        if (inBetweenTwoArcsFlag == false){

            unsigned int arcStart = arcIndices(previousArc);
            unsigned int arcFinish = arcIndices(previousArc+1);

            Eigen::VectorXd arcDataTA = trueAnomalyInput.block(arcStart,2,arcFinish-arcStart,1);
            Eigen::MatrixXd arcDataPosition = trueAnomalyInput.block(arcStart,3,arcFinish-arcStart,3);

            Eigen::VectorXd arcDataTADif = (arcDataTA - vehicleTrueAnomaly*Eigen::VectorXd::Ones(arcDataTA.size())).cwiseAbs();
            std::ptrdiff_t ind; arcDataTADif.minCoeff(&ind);

            interpolatedErrorMatrix.block(i,0,1,3) = arcDataPosition.block(ind,0,1,3);

        }

        // if in between two arcs,
        // find error corresponding to the given TA for both arcs and interpolate based on time
        else{

            std::vector<double> arct, arcx, arcy, arcz;

            for (int j = previousArc; j<previousArc+2; j++){
                unsigned int arcStart = arcIndices(j);
                unsigned int arcFinish = arcIndices(j+1);

                Eigen::VectorXd arcDataTA = trueAnomalyInput.block(arcStart,2,arcFinish-arcStart,1);
                Eigen::MatrixXd arcDataPosition = trueAnomalyInput.block(arcStart,3,arcFinish-arcStart,3);

                Eigen::VectorXd arcDataTADif = (arcDataTA - vehicleTrueAnomaly*Eigen::VectorXd::Ones(arcDataTA.size())).cwiseAbs();
                std::ptrdiff_t ind; arcDataTADif.minCoeff(&ind);

                arct.push_back(arcAverageTimes(j));
                arcx.push_back(arcDataPosition(ind,0));
                arcy.push_back(arcDataPosition(ind,1));
                arcz.push_back(arcDataPosition(ind,2));

            }

            LinearInterpolatorDouble xInterpolator2(arct, arcx, huntingAlgorithm, use_boundary_value);
            LinearInterpolatorDouble yInterpolator2(arct, arcy, huntingAlgorithm, use_boundary_value);
            LinearInterpolatorDouble zInterpolator2(arct, arcz, huntingAlgorithm, use_boundary_value);

            interpolatedErrorMatrix(i,0) = xInterpolator2.interpolate(currentTime);
            interpolatedErrorMatrix(i,1) = yInterpolator2.interpolate(currentTime);
            interpolatedErrorMatrix(i,2) = zInterpolator2.interpolate(currentTime);

            arct.clear(); arcx.clear(); arcy.clear(); arcz.clear();
        }

    if (interpolatedErrorMatrix(i,0) == 0.0 && interpolatedErrorMatrix(i,1) == 0.0 && interpolatedErrorMatrix(i,2) == 0.0){
        std::runtime_error("ERROR: satellite error cannot be equal to 0");
        std::cout<<"prev arc: "<<previousArc<<" flag: "<<inBetweenTwoArcsFlag<<std::endl;
        std::cout<<interpolatedErrorMatrix.block(i,0,1,3)<<std::endl;
    }

    }

    return interpolatedErrorMatrix;
}

Eigen::Matrix3d transformAngularMomentumFromLocalToGlobalFrame(
        const Eigen::Matrix3d angularMomentumInLocalFrame,
        const std::string localFrame,
        const std::string globalFrame,
        const double currentTime)
{
    using namespace tudat::spice_interface;

    Eigen::Matrix3d transformation =
            computeRotationMatrixBetweenFrames(globalFrame, localFrame, currentTime);

    Eigen::Matrix3d angularMomentumInGlobalFrame =
            transformation * angularMomentumInLocalFrame;

    return angularMomentumInGlobalFrame;
}

double phaseAccordingToSolarMinimum(
        const double solarMinimumEpoch,
        const double period)
{
    using namespace tudat::mathematical_constants;
    return (2.0*PI/period)*(solarMinimumEpoch+period/4.0);
}

double simpleSine(
        const double amplitude,
        const double period,
        const double phase,
        const double time){
    using namespace tudat::mathematical_constants;
    return amplitude*sin(2.0*PI*time/period+phase);
}

std::map< double, Eigen::MatrixXd > tabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime,
        std::function< double() > amplitudeFunction,
        std::function< double() > periodFunction,
        std::function< double() > phaseFunction)
{

    const double amplitude = amplitudeFunction();
    const double period = periodFunction();
    const double phase = phaseFunction();

    std::cout<<"Making table of Spherical Harmonics corrections with "
             <<"amplitude: "<<amplitude<<" / "
             <<"period: "<<period<<" / "
             <<"phase: "<<phase<<std::endl;

    std::map< double, Eigen::MatrixXd > coefficientCorrections;

    const double interval = 3600.0;
    double currentTime = initialTime - interval;

    Eigen::Matrix<double, 1, 1> currentCorrection = Eigen::Matrix<double, 1, 1>::Zero();

    while (currentTime <= finalTime + interval){
        currentCorrection(0,0) = simpleSine(amplitude, period, phase, currentTime);
        coefficientCorrections.insert(std::make_pair(currentTime, currentCorrection));
        currentTime += interval;

    }

//    std::string outputPath = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/SHtable/";
//    std::string amplitudeString = std::to_string(amplitude);

//    using namespace tudat::input_output;
//    writeDataMapToTextFile( coefficientCorrections, "cosineCoefficientCorrectionsAmp"+amplitudeString+".dat",
//                                          outputPath, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

    return coefficientCorrections;

}


std::map< double, Eigen::MatrixXd > zeroTabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime)
{

    std::map< double, Eigen::MatrixXd > coefficientCorrections;

    const double interval = 3600.0;
    double currentTime = initialTime - interval;

    Eigen::Matrix<double, 1, 1> currentCorrection = Eigen::Matrix<double, 1, 1>::Zero();

    while (currentTime <= finalTime + interval){
        currentCorrection(0,0) = 0.0;
        coefficientCorrections.insert(std::make_pair(currentTime, currentCorrection));
        currentTime += interval;
    }

    return coefficientCorrections;

}

// Montenbruck & Gill eq 8.42, executed slightly different to avoid square matrices with size = number of observations (result is tested to be identical)
Eigen::MatrixXd calculateConsiderCovarianceMatrix(
        const Eigen::MatrixXd P,
        const Eigen::VectorXd W_diagonal,
        const Eigen::MatrixXd C,
        const Eigen::MatrixXd Hx,
        const Eigen::MatrixXd Hc)
{

    const unsigned int numberOfObservations = Hx.rows();

    Eigen::MatrixXd term1 = P * Hx.transpose();
    std::cout<<"P*Hx calculated... ";

    for (unsigned int i = 0; i < numberOfObservations; i++){
        term1.col(i) *= W_diagonal(i);
    }
    std::cout<<"P*Hx multiplied with W... ";

    Eigen::MatrixXd term2a = term1 * Hc;
    std::cout<<"first brackets multiplied with Hc...";

    Eigen::MatrixXd term2c = Hc.transpose() * term1.transpose();
    std::cout<<"Hc^T multiplied with last brackets...";

    return P + term2a * C * term2c;
}


