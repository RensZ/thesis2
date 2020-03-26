
#include "MyApplications/customFunctions.h"

// Get path for output directory.

static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{
    // Declare file path string assigned to filePath.
    // __FILE__ only gives the absolute path of the header file!
    std::string filePath_( __FILE__ );

    // Strip filename from temporary string and return root-path string.
    std::string outputPath = filePath_.substr( 0, filePath_.length( )
                                - std::string( "MyApplications" ).length( ) )
                                + std::string( "/" );
    if( extraDirectory != "" ){outputPath += extraDirectory;}
    if( outputPath.at( outputPath.size( ) - 1 ) != '/' ){outputPath += "/";}

    return outputPath;
}


// convert calendar date and time to seconds after J2000
double secondsAfterJ2000(Eigen::Vector6i datetime){
    using namespace tudat;
    using namespace tudat::basic_astrodynamics;
    const unsigned int secondsPerDay = 60*60*24;
    const unsigned int julianDayJ2000 = 2451545;
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay(datetime[0],datetime[1],datetime[2],datetime[3],datetime[4],datetime[5]);
    return (julianDay - julianDayJ2000)*secondsPerDay;
}


//! Function to calculate the Sun-Mercury-Earth angle
//! this is a simplified approach which simply uses the mean motion
//! w.r.t. an epoch where both planets alligned
double mercurySunEarthAngle(const double time){

    // from NASA fact sheet https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    const double yearMercury = 88.0  *24.0*60.0*60.0;
    const double yearEarth   = 365.2 *24.0*60.0*60.0;

    // Get mean motion of planets
    const double pi = 3.14159265359;
    const double meanMotionMercury = 2*pi/yearMercury;
    const double meanMotionEarth   = 2*pi/yearEarth;

    // get time since last eclipse, when Mercury and Earth perfectly aligned (MSE = 0)
    Eigen::Vector6i timeAlignedVector;
    timeAlignedVector << 2019, 11, 11, 15, 20, 0; // YYYY, MM, DD, hh, mm, ss
    const double timeAligned = secondsAfterJ2000(timeAlignedVector);
    const double timeDelta = time-timeAligned;

    // get angular distance traveled since last eclipse (or in the past)
    const double angleTravelledMercury = timeDelta*meanMotionMercury;
    const double angleTravelledEarth   = timeDelta*meanMotionEarth;

    // get value between 0-180 degrees
    double relativeAngle = std::fmod((abs(angleTravelledMercury-angleTravelledEarth)),2*pi);
    if (relativeAngle > pi){
        relativeAngle = 2*pi-relativeAngle;
    }
    return relativeAngle;
}


//! Function to generate a list of observation times between an initial and final time, with a time step in between.
//! a maximum number of observations can be indicated, if the pattern generates more than that, a random selection will be deleted to meet the maximum.
//! a maximum Mercury-Sun-Earth angle can be indicated, observations will not be included above the requirement (see function above)
//! in addition times can be added to the list manually which will be added last
std::vector<double> makeObservationTimeList(const double initialTime,
                                            const double endTime,
                                            const double timeStep,
                                            const unsigned int maximumNumberOfObservations, //if not applicable, set to a very large number.
                                            const double maxMSEangle, //radians. if not applicable, set to value greater than pi.
                                            const std::vector<double> flybyTimes)
{
    std::vector< double > observationTimeList;
    double currentTime = initialTime;

    // add times to observation list which satisfy MSE requirements
    while (currentTime <= endTime){
        double MSEangle = mercurySunEarthAngle(currentTime);
        if (MSEangle < maxMSEangle){
            observationTimeList.push_back(currentTime);
        }
            currentTime += timeStep;
    }

    // if list is longer than the max number of observations, delete random entries
    srand(0);
    while ( observationTimeList.size() > maximumNumberOfObservations ){
        unsigned int numberOfObservations = observationTimeList.size();
        unsigned int randomIndex = rand()%(numberOfObservations-1);
        observationTimeList.erase( observationTimeList.begin() + randomIndex );
    }

    // add times which are manually indicated
    for( unsigned int i = 0; i < flybyTimes.size(); i++ )
        observationTimeList.push_back(flybyTimes.at(i));
    return observationTimeList;
}


//! Function to make the noise level dependent on the Mercury-Sun-Earth angle,
//! based on a linear relation constructed with the minimum and maximum noise.
double noiseBasedOnMSEangle(const double time){

    // calculate MSE angle with previous function
    double relativeAngle = mercurySunEarthAngle(time);

    // get noise level (linear function between a minimum and maximum point)
    const double pi = 3.14159265359;
    const double maxAngle = pi-35.0*pi/180.0;
    const double minAngle = 0;
    const double noiseAtMaxAngle = 3.0;
    const double noiseAtMinAngle = 0.5;

    const double slope = (noiseAtMaxAngle-noiseAtMinAngle)/(maxAngle-minAngle);
    const double intercept = noiseAtMaxAngle - slope*maxAngle;
    const double noise = slope*relativeAngle + intercept;

    // create a gaussian sample
    boost::random::mt19937 rng(time);
    boost::random::normal_distribution<> nd(0.0,noise);
    const double sample = nd(rng);

//    std::cout << "SPE (deg): " << relativeAngle*180/pi << " ---- "
//              << "noise level (m): " << noise << " ---- "
//              << "rng noise sample (m): " << sample << std::endl;

    return sample;
}

