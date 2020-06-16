 /*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved.
 *
 *    Author: Rens van der Zwaard
 *
 *    Changelog
 *      19-02-19    Created based on TUDAT template application
 *
 *    Notes
 *
 */

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>
#include <random>

#include "tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/timeVaryingGravitationalParameterPartial.h"

#include "json.hpp"

// custom functions written for the main application are placed here to save space:
#include "tudatApplications/thesis/MyApplications/customFunctions.h"


// Get path for output directory.
static inline std::string getOutputPath(
        const std::string& extraDirectory = "" )
{

    std::string outputPath = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/" + extraDirectory;
    return outputPath;
}

int main( )
{

    // Import namespaces
    using namespace tudat;
    using namespace tudat::observation_models;
    using namespace tudat::orbit_determination;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::propagators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::coordinate_conversions;
    using namespace tudat::ground_stations;
    using namespace tudat::observation_models;
    using namespace tudat::statistics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;
    using namespace tudat::physical_constants;
    using json = nlohmann::json;


    ////////////////////////////
    //// APPLICATION INPUTS ////
    ////////////////////////////

    // Acceleration settings
    const bool calculateDeSitterCorrection = false; //need to implement partials
    const bool includeTimeVaryingGravitationalMomentsSun = false;
    const unsigned int maximumDegreeSunGravitationalMomentVariation = 2;

    // Use multiple arcs for Mercury per mission
    const bool useMultipleMercuryArcs = true;

    // Parameter apriori values and uncertainties
    const bool useAprioriValues = true;

    // Planet propagation settings
    const bool propogatePlanets = false; // Propogate the other planets besides Mercury (NOTE: need observations for other planets, or LS can't find solutions for other planets)

    // Parameter estimation settings
    const unsigned int maximumNumberOfIterations = 5;
    const double sigmaPosition = 1000.0; //educated guess
    const double sigmaVelocity = 1.0; //educated guess

    // ABM integrator settings (if RK4 is used instead, initialstepsize is taken)
    const double initialTimeStep = 3600.0/2.0;
    const double minimumStepSize = 3600.0/2.0;
    const double maximumStepSize = 3600.0/2.0;
    const double relativeErrorTolerence = 1.0;
    const double absoluteErrorTolerence = 1.0;
    const unsigned int minimumOrder = 8;
    const unsigned int maximumOrder = 8;

    // Observation settings
    const bool includeSpacecraftPositionError = true;

    // Other planetary parameters, currently not included in json
    const double mercuryGravitationalParameter = (2.2031870798779644e+04)*(1E9); //m3/s2, from https://pgda.gsfc.nasa.gov/products/71
    const double sunRadius = 695.7E6; //m, from nasa fact sheet

    // Solar cycle
    Eigen::Vector6i solarMinimumDatetime;
    solarMinimumDatetime << 2008, 12, 15, 0, 0, 0;
    const double solarMinimumEpoch = secondsAfterJ2000(solarMinimumDatetime);
    const double solarCycleDuration = 11.0*physical_constants::JULIAN_YEAR;


    ////////////////////////
    //// MISSION INPUTS ////
    ////////////////////////

    // Load json input

    std::string input_filename = "inputs_MESSENGER_and_BepiColombo.json";
    std::cout<<"---- RUNNING SIMULATION FOR INPUTS WITH FILENAME: "<<input_filename<<" ----"<<std::endl;

    std::string json_directory = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/";
    std::ifstream json_file(json_directory + input_filename);
    json json_input;
    json_file >> json_input;

    std::cout<<"input values imported: "<<json_input<<std::endl;

    // Acceleration settings
    const bool calculateSchwarzschildCorrection = json_input["calculateSchwarzschildCorrection"];
    const bool calculateLenseThirringCorrection = json_input["calculateLenseThirringCorrection"];
    const bool includeSEPViolationAcceleration = json_input["includeSEPViolationAcceleration"];
    const bool includeTVGPAcceleration = json_input["includeTVGPAcceleration"];
    const bool estimateJ2Amplitude = json_input["estimateJ2Amplitude"];
    const bool estimateJ2Period = json_input["estimateJ2Period"];
    const bool estimateJ2Phase = json_input["estimateJ2Phase"];

    // Retrieve input parameters including uncertainties and apriori values
    const double sunAngularMomentum = json_input["sunAngularMomentum"];
    const double sunGravitationalParameter = json_input["sunGravitationalParameter"];
    const double timeVaryingGravitationalParameter = json_input["timeVaryingGravitationalParameter"];

    const double sigmaGamma = json_input["sigma_gamma"];
    const double sigmaBeta = json_input["sigma_beta"];
    const double sigmaAlpha1 = json_input["sigma_alpha1"];
    const double sigmaAlpha2 = json_input["sigma_alpha2"];
    const double sigmaNordtvedt = json_input["sigma_Nordtvedt"];
    const double sigmaSunGM = json_input["sigma_mu_Sun"];
    const double sigmaTVGP = json_input["sigma_TVGP"];

    const double unnormalisedSunJ2 = json_input["sunJ2"];
    const double sunJ2 = unnormalisedSunJ2 / calculateLegendreGeodesyNormalizationFactor(2,0);
    const double sigmaUnnormalisedSunJ2 = json_input["sigma_J2_Sun"];
    const double sigmaSunJ2 = sigmaUnnormalisedSunJ2 / calculateLegendreGeodesyNormalizationFactor(2,0);
    const double unnormalisedSunJ4 = json_input["sunJ4"];
    const double sunJ4 = unnormalisedSunJ4 / calculateLegendreGeodesyNormalizationFactor(4,0);
    const double unnormalisedSigmaSunJ4 = json_input["sigma_J4_Sun"];
    const double sigmaSunJ4 = unnormalisedSigmaSunJ4 / calculateLegendreGeodesyNormalizationFactor(4,0);


    // Use constraint nordtvedt=4*beta-gamma-3
    const bool useNordtvedtConstraint = json_input["useNordtvedtConstraint"];
    const bool estimatePPNalphas = json_input["estimatePPNalphas"];

    // Location of simulation output
    std::string outputSubFolderName = json_input["outputSubFolderName"];
    std::string outputSubFolder = getOutputPath( ) + outputSubFolderName;
    if (useMultipleMercuryArcs) {outputSubFolder = outputSubFolder + "_multiarc";}

    // retreive parameters per mission
    std::vector< double > initialTimeVector;
    std::vector< double > finalTimeVector;
    std::vector< std::string > vehicleVector;

    std::vector< double > observationInitialTimeVector;
    std::vector< double > observationFinalTimeVector;
    std::vector< double > observationTimeStepVector;
    std::vector< double > trackingArcDurationVector;
    std::vector< unsigned int > maximumNumberOfTrackingDaysVector;
    std::vector< std::vector<double> > flybyListVector;

    std::vector< double > noiseAtMinAngleVector;
    std::vector< double > noiseAtMaxAngleVector;
    std::vector< double > maxMSEAngleDegVector;

    std::vector< std::string > missionJSONFiles = json_input["missionJSONFiles"];
    const unsigned int numberOfMissions = missionJSONFiles.size();
    for (unsigned int m=0; m<numberOfMissions; m++){

        // retreive mission JSON data
        std::ifstream json_file_mission(json_directory + missionJSONFiles.at(m));
        json json_input_mission;
        json_file_mission >> json_input_mission;

        vehicleVector.push_back(json_input_mission["vehicle"]);

        // Retrieve simulation start and end time
        std::vector<int> json1 = json_input_mission["initialTime"];
        Eigen::Vector6i initialTime(json1.data()); json1.clear();
        initialTimeVector.push_back(secondsAfterJ2000(initialTime));
        std::vector<int> json2 = json_input_mission["finalTime"];
        Eigen::Vector6i finalTime(json2.data()); json2.clear();
        finalTimeVector.push_back(secondsAfterJ2000(finalTime));

        // retreive regular observation schedule
        std::vector<int> json3 = json_input_mission["observationInitialTime"];
        Eigen::Vector6i observationInitialTime(json3.data()); json3.clear();
        observationInitialTimeVector.push_back(secondsAfterJ2000(observationInitialTime));
        std::vector<int> json4 = json_input_mission["observationFinalTime"];
        Eigen::Vector6i observationFinalTime(json4.data()); json4.clear();
        observationFinalTimeVector.push_back(secondsAfterJ2000(observationFinalTime));
        observationTimeStepVector.push_back(json_input_mission["observationTimeStep"]);
        trackingArcDurationVector.push_back(json_input_mission["trackingArcDuration"]);
        maximumNumberOfTrackingDaysVector.push_back(json_input_mission["maximumNumberOfTrackingDays"]);

        // retreive flybys
        std::vector<int> flybyObject = json_input_mission["flybyObject"];
        std::vector<double> flybyList;
        std::vector<int> currentFlyby;
        for (unsigned int i=0; i<flybyObject.size()/6; i++ ){

            for (unsigned int c=0; c<6; c++ ){
                currentFlyby.push_back(flybyObject.at(6*i+c));
            }
            flybyList.push_back(secondsAfterJ2000(currentFlyby));
            currentFlyby.clear();
        }
        flybyListVector.push_back(flybyList);

        // observation Noise
        noiseAtMinAngleVector.push_back(json_input_mission["noiseAtMinAngle"]);
        noiseAtMaxAngleVector.push_back(json_input_mission["noiseAtMaxAngle"]);
        maxMSEAngleDegVector.push_back(json_input_mission["maxMSEAngleDeg"]);

    }


    // settings time variable gravitational moments
    double J2amplitude = sunJ2/2.0;
//        double J2amplitude = 0.005E-7 / calculateLegendreGeodesyNormalizationFactor(2,0); // currently from Antia et al 2008
    double J2period = solarCycleDuration;
    double J2phase = phaseAccordingToSolarMinimum(solarMinimumEpoch, J2period);

    double J4amplitude = 0.06E-9 / calculateLegendreGeodesyNormalizationFactor(4,0); // currently from Antia et al 2008
    double J4period = solarCycleDuration;
    double J4phase = phaseAccordingToSolarMinimum(solarMinimumEpoch, J4period); //antiphase?

    const unsigned int maximumDegreeWithInput = 4;

    Eigen::VectorXd valuesSunGravitationalMoments = Eigen::VectorXd::Zero(maximumDegreeWithInput+1);
    valuesSunGravitationalMoments(0) = 1.0; //central gravity
    valuesSunGravitationalMoments(2) = sunJ2; valuesSunGravitationalMoments(4) = sunJ4;

    Eigen::VectorXd sigmaValuesSunGravitationalMoments = Eigen::VectorXd::Zero(maximumDegreeWithInput+1);
    sigmaValuesSunGravitationalMoments(2) = sigmaSunJ2; sigmaValuesSunGravitationalMoments(4) = sigmaSunJ4;

    Eigen::VectorXd amplitudesSunGravitationalMomentsVariation = Eigen::VectorXd::Zero(maximumDegreeWithInput+1);
    amplitudesSunGravitationalMomentsVariation(2) = J2amplitude; amplitudesSunGravitationalMomentsVariation(4) = J4amplitude;

    Eigen::VectorXd periodsSunGravitationalMomentsVariation = Eigen::VectorXd::Zero(maximumDegreeWithInput+1);
    periodsSunGravitationalMomentsVariation(2) = J2period; periodsSunGravitationalMomentsVariation(4) = J4period;

    Eigen::VectorXd phasesSunGravitationalMomentsVariation = Eigen::VectorXd::Zero(maximumDegreeWithInput+1);
    phasesSunGravitationalMomentsVariation(2) = J2phase; phasesSunGravitationalMomentsVariation(4) = J4phase;

    const std::shared_ptr< InterpolatorSettings > interpolatorSettings =
            std::make_shared< InterpolatorSettings >( cubic_spline_interpolator ); //interpolator of Sun SH coefficient variation



    /////////////////////
    //// ENVIRONMENT ////
    /////////////////////

    // Load spice kernels.
    std::cout << "loading SPICE kernels..." << std::endl;

    std::string kernelsPath = input_output::getSpiceKernelPath( );

    spice_interface::loadStandardSpiceKernels( );

//    std::vector< std::string > customKernels;
//    customKernels.push_back( kernelsPath + "tudat_merged_spk_kernel_thesis2.bsp" );
//    spice_interface::loadStandardSpiceKernels( customKernels );

    std::cout << "building environment..." << std::endl;


    // Define bodies
    unsigned int totalNumberOfBodies = 8;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Sun";
    bodyNames[ 1 ] = "Mercury";
    bodyNames[ 2 ] = "Venus";
    bodyNames[ 3 ] = "Earth";
    bodyNames[ 4 ] = "Mars";
    bodyNames[ 5 ] = "Jupiter";
    bodyNames[ 6 ] = "Saturn";
//    bodyNames[ 7 ] = "Uranus";
//    bodyNames[ 8 ] = "Neptune";
    bodyNames[ 7 ] = "Moon";

    // load SPICE settings
    double initialSimulationTime = *std::min_element(initialTimeVector.begin(), initialTimeVector.end());
    double finalSimulationTime = *std::max_element(finalTimeVector.begin(), finalTimeVector.end());
    double buffer = maximumOrder*maximumStepSize; //see Tudat libraries 1.1.3.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;

    // Default body settings
    bodySettings = getDefaultBodySettings( bodyNames, initialSimulationTime - buffer, finalSimulationTime + buffer);

    std::cout << "setting custom settings Sun..." << std::endl;

    // Nominal values spherical harmonics Sun
    Eigen::MatrixXd normalizedSineCoefficients
            = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMomentVariation + 1,
                                    maximumDegreeSunGravitationalMomentVariation + 1);
    Eigen::MatrixXd normalizedCosineCoefficients
            = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMomentVariation + 1,
                                    maximumDegreeSunGravitationalMomentVariation + 1);

    if (normalizedCosineCoefficients.col(0).size() != valuesSunGravitationalMoments.size()){
        std::runtime_error("error: vector sizes incompatible");
    }

    normalizedCosineCoefficients.col(0) = valuesSunGravitationalMoments.transpose();

    bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                sunGravitationalParameter, sunRadius,
                normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );


    // Time varying spherical harmonics coefficients Sun
    if (includeTimeVaryingGravitationalMomentsSun){

        std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

        relativity::variableJ2Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(2));
        relativity::variableJ2Interface->setPeriod(periodsSunGravitationalMomentsVariation(2));
        relativity::variableJ2Interface->setPhase(phasesSunGravitationalMomentsVariation(2));

        std::function< double() > amplitudeFunction =
                std::bind( &relativity::VariableJ2Interface::getAmplitude, relativity::variableJ2Interface );
        std::function< double() > periodFunction =
                std::bind( &relativity::VariableJ2Interface::getPeriod, relativity::variableJ2Interface );
        std::function< double() > phaseFunction =
                std::bind( &relativity::VariableJ2Interface::getPhase, relativity::variableJ2Interface );


        std::function< std::map< double, Eigen::MatrixXd >( ) > tabulatedSphericalHarmonicsCoefficientCorrectionsFunction =
                std::bind(tabulatedSphericalHarmonicsCoefficientCorrections,
                          initialSimulationTime, finalSimulationTime,
                          amplitudeFunction, periodFunction, phaseFunction);

        std::map< double, Eigen::MatrixXd > sineCoefficientCorrections =
                zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

        const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettings =
                std::make_shared< TabulatedGravityFieldVariationSettingsWithCosineFunction >(
                    tabulatedSphericalHarmonicsCoefficientCorrectionsFunction,
                    sineCoefficientCorrections,
                    2, 0, interpolatorSettings,
                    initialSimulationTime, finalSimulationTime, 3600.0);

        gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettings );


//            for (unsigned int d=2; d<=maximumDegreeSunGravitationalMomentVariation; d+=2){

//                std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections =
//                        tabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime,
//                                                                          amplitudesSunGravitationalMomentsVariation(d),
//                                                                          periodsSunGravitationalMomentsVariation(d),
//                                                                          phasesSunGravitationalMomentsVariation(d));
//                std::map< double, Eigen::MatrixXd > sineCoefficientCorrections =
//                        zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

//                input_output::writeDataMapToTextFile( cosineCoefficientCorrections, "cosineCoefficientCorrectionsJ"+std::to_string(d)+".dat",
//                                                      outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );
//                input_output::writeDataMapToTextFile( sineCoefficientCorrections, "sineCoefficientCorrectionsJ"+std::to_string(d)+".dat",
//                                                      outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

//                const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettings =
//                        std::make_shared< TabulatedGravityFieldVariationSettings >( cosineCoefficientCorrections,
//                                                                                    sineCoefficientCorrections,
//                                                                                    2, 0, interpolatorSettings,
//                                                                                    initialSimulationTime, finalSimulationTime, 3600.0);

//                gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettings );
//            }

        bodySettings[ "Sun" ]->gravityFieldVariationSettings = gravityFieldVariationSettings;
    }

    // Prepare angular momentum vector Sun
    const Eigen::Vector3d sunAngularMomentumVectorInSunFrame(0.0, 0.0, sunAngularMomentum);
    const Eigen::Vector3d sunAngularMomentumVectorPerUnitMassInSunFrame =
            sunAngularMomentumVectorInSunFrame /
            (sunGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT);

//    std::cout<<computeRotationMatrixBetweenFrames("ECLIPJ2000","IAU_Sun",initialSimulationTime+100000.0)<<std::endl;
//    std::cout<<computeRotationMatrixBetweenFrames("ECLIPJ2000","IAU_Sun",(initialSimulationTime+finalSimulationTime)/2.0)<<std::endl;
//    std::cout<<computeRotationMatrixBetweenFrames("ECLIPJ2000","IAU_Sun",finalSimulationTime)<<std::endl;


    // Create body map
    NamedBodyMap bodyMap = createBodies( bodySettings );

    if (useMultipleMercuryArcs){
        bodyMap[ "Mercury" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                                std::map< double, std::shared_ptr< Ephemeris > >( ), "SSB", "ECLIPJ2000" ) );
    }

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////
    //// ACCELERATIONS ////
    ///////////////////////

    std::cout << "defining accelerations..." << std::endl;

    // Define list of bodies to propagate
    std::vector< std::string > bodiesToPropagate;
    if (propogatePlanets == true){
        bodiesToPropagate = bodyNames;
    } else{
        bodiesToPropagate.push_back("Mercury");
    }
    unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );


    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.resize( numberOfNumericalBodies );

    // Set central body as Solar System Barycenter for Sun and each planet, for Moon use the Earth
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        if (bodiesToPropagate[i] == "Moon"){
            centralBodies[i] = "Earth";
        } else{
            centralBodies[i] = "SSB";
        }
    }

    relativity::ppnParameterSet->setParameterGamma(1.0);
    relativity::ppnParameterSet->setParameterBeta(1.0);
    relativity::ppnParameterSet->setParameterAlpha1(0.0);
    relativity::ppnParameterSet->setParameterAlpha2(0.0);

    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerations;
        for( unsigned int j = 0; j < totalNumberOfBodies; j++ )
        {
            if( bodiesToPropagate.at( i ) != bodyNames.at( j ) )
            {
                if (bodyNames.at( j ) == "Sun"){
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< SphericalHarmonicAccelerationSettings > (maximumDegreeSunGravitationalMomentVariation,0));

                    if (calculateSchwarzschildCorrection || calculateLenseThirringCorrection || calculateDeSitterCorrection){
                        currentAccelerations[ bodyNames.at( j ) ].push_back(
                                    std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                        calculateSchwarzschildCorrection,
                                        calculateLenseThirringCorrection,
                                        calculateDeSitterCorrection,
                                        "",
                                        sunAngularMomentumVectorPerUnitMassInSunFrame));
                    }

                    if (includeSEPViolationAcceleration == true){
                        currentAccelerations[ bodyNames.at( j ) ].push_back(
                                    std::make_shared< SEPViolationAccelerationSettings >(
                                        bodyNames, useNordtvedtConstraint));
                    }

                    if (includeTVGPAcceleration == true){
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< TimeVaryingGravitationalParameterAccelerationSettings >(
                                    timeVaryingGravitationalParameter));
                    }

                }
                else{
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( central_gravity ) );
                }
            }
        }
        accelerationMap[ bodiesToPropagate.at( i ) ] = currentAccelerations;
    }


    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );




    /////////////////////////////
    //// DEPENDENT VARIABLES ////
    /////////////////////////////

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    for( unsigned int i = 0; i < totalNumberOfBodies; i++ ){
        if (!(bodyNames.at( i ) == "Sun") && !(bodyNames.at( i ) == "Mercury")){
            dependentVariablesList.push_back(
                      std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            central_gravity, "Mercury" , bodyNames.at( i ) ) );
        }
    }

    dependentVariablesList.push_back(
                std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                "Mercury", "Sun", maximumDegreeSunGravitationalMomentVariation, 0 ) );

    if (calculateSchwarzschildCorrection){
        dependentVariablesList.push_back(
                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 1 ) );
        dependentVariablesList.push_back(
                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 2 ) );
    }

    if (calculateLenseThirringCorrection){
        dependentVariablesList.push_back(
                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 3 ) );
    }

    if (calculateDeSitterCorrection){
        dependentVariablesList.push_back(
                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 4 ) );
    }

    if (includeSEPViolationAcceleration){
    dependentVariablesList.push_back(
              std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    sep_violation_acceleration, "Mercury" , "Sun" ) );
    }

    if (includeTVGPAcceleration){
    dependentVariablesList.push_back(
              std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    time_varying_gravitational_parameter_acceleration, "Mercury" , "Sun" ) );
    }

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
    std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


    //////////////////////////////
    //// PROPAGATION SETTINGS ////
    //////////////////////////////

    std::cout << "defining propagation settings..." << std::endl;

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    Eigen::VectorXd systemInitialState;
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings;

    if (useMultipleMercuryArcs){
        std::vector< Eigen::VectorXd > arcInitialStates;

        for (unsigned int m=0; m<numberOfMissions; m++){

            Eigen::VectorXd currentArcInitialState = getBodyCartesianStateAtEpoch(
                        bodiesToPropagate.at(0), centralBodies.at(0),
                        "ECLIPJ2000", "None", initialTimeVector.at(m) );

            std::cout<<"initial state: "<<currentArcInitialState.transpose()<<std::endl;
            arcInitialStates.push_back( currentArcInitialState );

            std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                  std::make_shared< propagators::PropagationTimeTerminationSettings >(
                        finalTimeVector.at(m), true );

            propagatorSettingsList.push_back(
                        std::make_shared< TranslationalStatePropagatorSettings< double > >(
                            centralBodies, accelerationModelMap, bodiesToPropagate,
                            currentArcInitialState, terminationSettings, cowell, dependentVariablesToSave ) );
        }

        // Create propagator settings
        propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );

    } else{

        // Get initial state of bodies to be propagated
        systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodyMap, initialSimulationTime );

        // Define propagator settings.
        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
              std::make_shared< propagators::PropagationTimeTerminationSettings >( finalSimulationTime, true );

        propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate,
                  systemInitialState, terminationSettings, cowell, dependentVariablesToSave);

    }


    std::shared_ptr< AdamsBashforthMoultonSettings< double > > integratorSettings =
            std::make_shared< AdamsBashforthMoultonSettings< double > > (
                initialSimulationTime, initialTimeStep,
                minimumStepSize, maximumStepSize,
                relativeErrorTolerence, absoluteErrorTolerence,
                minimumOrder, maximumOrder);

//    std::shared_ptr< IntegratorSettings< > > integratorSettings =
//            std::make_shared< IntegratorSettings< > >
//            ( rungeKutta4, initialSimulationTime, initialTimeStep );




    ////////////////////////////
    //// DYNAMICS SIMULATOR ////
    ////////////////////////////

    std::cout << "running dynamics simulator..." << std::endl;

    std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
    std::vector< std::map< double, Eigen::VectorXd > > spiceStatesAtPropagationTimes;
    allBodiesPropagationHistory.resize( bodiesToPropagate.size() );
    spiceStatesAtPropagationTimes.resize( bodiesToPropagate.size() );

    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    if (useMultipleMercuryArcs){

        MultiArcDynamicsSimulator <> dynamicsSimulator (bodyMap,
                                                        integratorSettings,
                                                        propagatorSettings,
                                                        initialTimeVector,
                                                        true,false,true);

        std::vector< std::map< double, Eigen::VectorXd > > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::vector< std::map< double, Eigen::VectorXd > > dependentVariablesResult = dynamicsSimulator.getDependentVariableHistory();

        for( unsigned int m = 0; m < numberOfMissions; m++ ){

            std::map< double, Eigen::VectorXd > intermediateIntegrationResult = integrationResult.at( m );

            // Retrieve numerically integrated state for each body.
            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
                 stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
            {
                for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
                {
                    allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
                    spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
                            getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);

                }
            }

            std::map< double, Eigen::VectorXd > intermediateDependentVariablesResult = dependentVariablesResult.at( m );

            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateDependentVariablesResult.begin( );
                 stateIterator != intermediateDependentVariablesResult.end( ); stateIterator++ )
            {
                dependentVariablesHistory[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
            }
        }

    } else{

        SingleArcDynamicsSimulator <> dynamicsSimulator (bodyMap,integratorSettings,propagatorSettings,true,false,true);
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory();

        // Retrieve numerically integrated state for each body.
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( ); stateIterator++ )
        {
            for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
            {
                allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
                spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
                        getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);

            }
        }

    }

    for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
    {
        // Write propagation history to file.
        input_output::writeDataMapToTextFile(
                    allBodiesPropagationHistory[ i ],
                    "StatePropagationHistory" + bodiesToPropagate.at( i ) + ".dat",
                    outputSubFolder,
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );

        input_output::writeDataMapToTextFile(
                    spiceStatesAtPropagationTimes[ i ],
                    "spiceStatesAtPropagationTimes" + bodiesToPropagate.at( i ) + ".dat",
                    outputSubFolder,
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );
    }

    // Write dependent variables history to file.
    input_output::writeDataMapToTextFile(
                dependentVariablesHistory,
                "DependentVariablesHistory.dat",
                outputSubFolder,
                "",
                std::numeric_limits< double >::digits10,
                std::numeric_limits< double >::digits10,
                "," );







    ////////////////////////////////
    //// CREATE GROUND STATIONS ////
    ////////////////////////////////

    std::cout << "creating ground stations..." << std::endl;


    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;

    groundStationNames.push_back( "Station1" );
    createGroundStation( bodyMap.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );


    ///////////////////////////////////////////
    //// DEFINE LINK ENDS FOR OBSERVATIONS ////
    ///////////////////////////////////////////

    std::cout << "defining link ends for observations..." << std::endl;


    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    LinkEnds linkEnds;
    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {

        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Mercury", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );
        linkEnds.clear( );

        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Mercury", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
        linkEnds.clear( );
    }

    LinkEnds twoWayRangeLinkEnds;
    twoWayRangeLinkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( 0 ) );
    twoWayRangeLinkEnds[ reflector1 ] = std::make_pair( "Mercury", "" );
    twoWayRangeLinkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( 0 ) );

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;

    linkEndsPerObservable[ n_way_range ].push_back( twoWayRangeLinkEnds );




    /////////////////////////////////////
    //// CREATE OBSERVATION SETTINGS ////
    /////////////////////////////////////

    std::cout << "creating observation settings..." << std::endl;

    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                               lightTimePerturbingBodies ) );

    observation_models::ObservationSettingsMap observationSettingsMap;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define settings for observable, with light-time corrections, no biases for selected 1-way range links
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           std::make_shared< ObservationSettings >( currentObservable, lightTimeCorrectionSettings ) ) );
        }
    }






    ////////////////////////////////
    //// ESTIMATABLE PARAMETERS ////
    ////////////////////////////////

    std::cout << "defining parameters to estimate..." << std::endl;

    std::vector< std::shared_ptr < EstimatableParameterSettings > > parameterNames;
    std::vector<double> varianceVector;

    // Add bodies that will be propagated to the parameters to be estimated
    if (useMultipleMercuryArcs){
        Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * initialTimeVector.size( ) );
        for( unsigned int m = 0; m < numberOfMissions; m++ )
        {
            systemInitialState.segment( m * 6, 6 ) = propagatorSettingsList.at( m )->getInitialStates( );

            for( unsigned int i = 0; i < 3; i++ ){
                varianceVector.push_back(sigmaPosition*sigmaPosition);
            }
            for( unsigned int i = 3; i < 6; i++ ){
                varianceVector.push_back(sigmaVelocity*sigmaVelocity);
            }
        }

        parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                      "Mercury", systemInitialState, initialTimeVector, "SSB" ) );
    }
    else{
        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
            int j = 6*i;
            parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                          bodiesToPropagate[i], systemInitialState.segment(j,6), centralBodies[i] ) );

            for( unsigned int i = 0; i < 3; i++ ){
                varianceVector.push_back(sigmaPosition*sigmaPosition);
            }
            for( unsigned int i = 3; i < 6; i++ ){
                varianceVector.push_back(sigmaVelocity*sigmaVelocity);
            }
        }
    }


    // relativistic parameters
    if (calculateSchwarzschildCorrection == true
        || includeSEPViolationAcceleration == true){
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_gamma ) );
    varianceVector.push_back(sigmaGamma*sigmaGamma);

    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_beta ) );
    varianceVector.push_back(sigmaBeta*sigmaBeta);
    }

    // Nordtvedt parameter
    if (includeSEPViolationAcceleration == true && useNordtvedtConstraint == false){
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_nordtvedt_parameter ) );
    varianceVector.push_back(sigmaNordtvedt*sigmaNordtvedt);
    }

    if (estimatePPNalphas == true){
        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("global_metric", ppn_parameter_alpha1 ) );
        varianceVector.push_back(sigmaAlpha1*sigmaAlpha1);

        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("global_metric", ppn_parameter_alpha2 ) );
        varianceVector.push_back(sigmaAlpha2*sigmaAlpha2);
    }

    // time varying gravitational parameter
    if (includeTVGPAcceleration == true){
        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("global_metric", time_varying_gravitational_parameter));
        varianceVector.push_back(sigmaTVGP*sigmaTVGP);
    }


    // gravity field Sun (mu and spherical harmonics)
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("Sun", gravitational_parameter));
    varianceVector.push_back(sigmaSunGM*sigmaSunGM);

    std::vector< std::pair< int, int > > blockIndices;
    for (unsigned int d=2; d<=maximumDegreeSunGravitationalMomentVariation; d+=2){
        blockIndices.push_back(std::make_pair(d,0));
        varianceVector.push_back(sigmaValuesSunGravitationalMoments(d)*sigmaValuesSunGravitationalMoments(d));
    }
    parameterNames.push_back(std::make_shared<SphericalHarmonicEstimatableParameterSettings>
                             (blockIndices,"Sun",spherical_harmonics_cosine_coefficient_block));


    // time varying J2 Sun
    double variableJ2parametersVariance = 0.5; //sigma taken as a percentage of the mean
    if (estimateJ2Amplitude){
        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("Sun", variable_J2_amplitude));
        double sigmaJ2Amplitude = variableJ2parametersVariance * relativity::variableJ2Interface->getAmplitude();
        varianceVector.push_back(sigmaJ2Amplitude*sigmaJ2Amplitude);
    }
    if (estimateJ2Period){
        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("Sun", variable_J2_period));
        double sigmaJ2Period = variableJ2parametersVariance * relativity::variableJ2Interface->getAmplitude();
        varianceVector.push_back(sigmaJ2Period*sigmaJ2Period);
    }
    if (estimateJ2Phase){
        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                 ("Sun", variable_J2_phase));
        double sigmaJ2Phase = variableJ2parametersVariance * relativity::variableJ2Interface->getAmplitude();
        varianceVector.push_back(sigmaJ2Phase*sigmaJ2Phase);
    }



    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal
    printEstimatableParameterEntries( parametersToEstimate );


    ////////////////////////////////////////
    //// RUN ORBITDETERMINATION MANAGER ////
    ////////////////////////////////////////

    std::cout << "Running OD manager..." << std::endl;

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );



    ////////////////////////////////////////
    //// DEFINING OD SIMULATOR SETTINGS ////
    ////////////////////////////////////////

    std::cout << "Defining observation settings..." << std::endl;



    // generate list of observations


    std::vector< double > baseTimeList;
    std::vector< std::vector< double > > seperateBaseTimeLists;
    for( unsigned int m = 0; m < numberOfMissions; m++ ){

        std::vector< double > currentMissionBaseTimeList =
                makeObservationTimeList(observationInitialTimeVector.at(m),
                                        observationFinalTimeVector.at(m),
                                        observationTimeStepVector.at(m),
                                        trackingArcDurationVector.at(m),
                                        maximumNumberOfTrackingDaysVector.at(m),
                                        unit_conversions::convertDegreesToRadians(maxMSEAngleDegVector.at(m)),
                                        flybyListVector.at(m));

        seperateBaseTimeLists.push_back(currentMissionBaseTimeList);

        baseTimeList.insert( baseTimeList.end(),
                             currentMissionBaseTimeList.begin(),
                             currentMissionBaseTimeList.end() );

    }


    // Create measurement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        // Define observable type and link ends
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        // Define observation times and reference link ends

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            if (currentObservable == position_observable){
                measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                        std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( observed_body, baseTimeList );
            }
            else{
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( receiver, baseTimeList );
            }
        }
    }


    // Create noise functions per observable

    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;

    std::function< double( const double ) > mercuryOrbiterNoiseFunction;
    mercuryOrbiterNoiseFunction = [noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists](const double time){
        return noiseSampleBasedOnMSEangleForMultipleMissions
                (time, noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists);
    };

    noiseFunctions[ one_way_range ] = mercuryOrbiterNoiseFunction;
    noiseFunctions[ n_way_range ] = mercuryOrbiterNoiseFunction;



    ///////////////////////////////
    //// SIMULATE OBSERVATIONS ////
    ///////////////////////////////

    std::cout << "Simulating observations..." << std::endl;

    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
                measurementSimulationInput,
                orbitDeterminationManager.getObservationSimulators( ),
                noiseFunctions);

    // similar container, but the "observation" will be the noise value instead of the actual observation
    PodInputDataType observationWeightsAndTimes = observationsAndTimes;

    // add spacecraft initial position error to the observations
    std::map<double, Eigen::Vector3d> saveSatelliteError;
    if (includeSpacecraftPositionError == true){

        std::cout << "Adding satellite estimation initial position error..." << std::endl;
//        Eigen::Vector3d constantSatelliteError; constantSatelliteError << 10.0, 10.0, 10.0;

        std::vector< Eigen::MatrixXd > interpolatedErrorMatrices;
        for (unsigned int m=0; m<numberOfMissions; m++){

            std::string vehicleErrorFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/error_inputs_"+vehicleVector.at(m)+".txt";
            Eigen::MatrixXd error_input = input_output::readMatrixFromFile(vehicleErrorFilename, ",");
            interpolatedErrorMatrices.push_back( interpolatePositionErrorsBasedOnTrueAnomaly(
                        error_input, seperateBaseTimeLists.at(m),
                        vehicleVector.at(m), mercuryGravitationalParameter) );
        }

        std::random_device rd;
        std::mt19937 gen(rd());

        // get to the location in the map where we can find the range observables
        PodInputDataType::iterator podInputIterator = observationsAndTimes.begin();
        PodInputDataType::iterator weightsIterator = observationWeightsAndTimes.begin();
        while (podInputIterator != observationsAndTimes.end()){

            if (podInputIterator->first == one_way_range || podInputIterator->first == n_way_range){

                SingleObservablePodInputType::iterator singleObservableIterator = podInputIterator->second.begin();
                SingleObservablePodInputType::iterator weightsIterator2 = weightsIterator->second.begin();
                while (singleObservableIterator != podInputIterator->second.end()){

                    // retrieve observations and their respective times for current observable type
                    ObservationVectorType allObservations = singleObservableIterator->second.first;
                    std::vector< double > allObservationTimes = singleObservableIterator->second.second.first;

                    ObservationVectorType newObservations = Eigen::VectorXd(allObservations.size());
                    ObservationVectorType observationWeights = Eigen::VectorXd(allObservations.size());

                    // for every observation, retrieve and add the range bias that should be added
                    for (unsigned int i=0; i<allObservationTimes.size(); i++){

                        double observationTime = allObservationTimes.at( i );

                        unsigned int check = 0;
                        Eigen::MatrixXd interpolatedErrorMatrix;

                        for (unsigned int m = 0; m<numberOfMissions; m++){
                            std::vector< double > currentBaseTimeList = seperateBaseTimeLists.at(m);

                            if (std::find(currentBaseTimeList.begin(), currentBaseTimeList.end(), observationTime)
                                    != currentBaseTimeList.end()) {

                                interpolatedErrorMatrix = interpolatedErrorMatrices.at(m);
                                check += 1;
                            }
                        }

                        if ( check == 0){
                            std::cout<<"finding noise value failed for t = "<<observationTime<<std::endl;
                            std::runtime_error("time not found in any mission observation list");
                        }
                        if ( check > 1){
                            std::cout<<"finding noise value failed for t = "<<observationTime
                                     <<" check = "<<check<<std::endl;
                            std::runtime_error("time found in multiple mission observation lists");
                        }

                        Eigen::Vector3d currentSatelliteError = interpolatedErrorMatrix.row(i);
                        Eigen::Vector3d randomErrorSample;
                        for (unsigned int j=0; j<3; j++){
                            std::normal_distribution<double> d(0.0, currentSatelliteError( j ));
                            randomErrorSample( j ) = d(gen);
                        }

                        Eigen::Vector3d mercuryPositionWrtEarth = -getBodyCartesianStateAtEpoch("Earth","Mercury","IAU_Mercury","None",observationTime).segment(0,3);
                        Eigen::Vector3d rangeUnitVector = mercuryPositionWrtEarth / mercuryPositionWrtEarth.norm( );

                        double rangeCorrection = randomErrorSample.dot(rangeUnitVector);
                        if (podInputIterator->first == n_way_range){rangeCorrection *= 2.0;}

                        newObservations(i) = allObservations(i) + rangeCorrection;

                        double noiseLevel = abs(rangeCorrection) + noiseSampleBasedOnMSEangleForMultipleMissions (observationTime,
                             noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists);

                        observationWeights(i) = 1.0/(noiseLevel*noiseLevel);

//                        std::cout<<observationTime<<" // "<<
//                                   currentSatelliteError.transpose()<<" // "<<
//                                   randomErrorSample.transpose()<<" // "<<
//                                   rangeCorrection<<" // "<<
//                                   noiseLevel<<std::endl;

                        saveSatelliteError.insert(std::make_pair(observationTime,currentSatelliteError));
                    }

                    singleObservableIterator->second.first = newObservations;
                    weightsIterator2->second.first = observationWeights;
                    singleObservableIterator++; weightsIterator2++;
                }
            }
            podInputIterator++; weightsIterator++;
        }
    }



    /////////////////////////////
    //// ESTIMATE PARAMETERS ////
    /////////////////////////////

    std::cout << "Estimating parameters..." << std::endl;

    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

////    srand(time(NULL));  //warning: not reproducable
    srand(0);

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ){
        // perturb body positions by  1 meters
        parameterPerturbation.segment(i*6,3) = Eigen::Vector3d::Constant( 1.0 );
        // perturb body velocities by 0.001 m/s
        parameterPerturbation.segment(i*6+3,3) = Eigen::Vector3d::Constant( 0.001 );
    }

//    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ){
//        // perturb body positions by random value between -10 and 10 meters
//        parameterPerturbation.segment(i*6,3) = Eigen::Vector3d( (rand()%200-100.0)/10.0,
//                                                                (rand()%200-100.0)/10.0,
//                                                                (rand()%200-100.0)/10.0 );
//        // perturb body velocities by random value between -0.01 and 0.01 m/s
//        parameterPerturbation.segment(i*6+3,3) = Eigen::Vector3d( (rand()%20-10.0)/100.0,
//                                                                (rand()%20-10.0)/100.0,
//                                                                (rand()%20-10.0)/100.0 );
//    }

//    // perturb parameters with a value between -0.001 and +0.001 times the given apriori sigma
////    for( unsigned int i = 0; i < truthParameters.size()-6*numberOfNumericalBodies; i++ ){
//    for( unsigned int i = 0; i < truthParameters.size(); i++ ){
////        unsigned int j = i + 6*numberOfNumericalBodies;
//        double randomPercentage = (rand()%20-10.0)/(1000.0);
//        std::cout << "randomPercentage: " << randomPercentage;
//        parameterPerturbation( i ) = sqrt(varianceVector.at( i ))*randomPercentage;
//        std::cout<< " // std: "<<sqrt(varianceVector.at( i )) <<" // perturbation: "<<parameterPerturbation( i ) <<std::endl;
//    }

    initialParameterEstimate += parameterPerturbation;


    std::cout << "True parameter values:" << std::endl;
    std::cout << truthParameters.transpose() << std::endl;

    std::cout << "Initial guesses:" << std::endl;
    std::cout << initialParameterEstimate.transpose() << std::endl;


    // Define a priori covariance matrix
    Eigen::MatrixXd aprioriCovariance =
        Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));

    Eigen::MatrixXd inverseOfAprioriCovariance;

    if (useAprioriValues == true){
        for( unsigned int i = 0; i < truthParameters.size(); i++ ){
            aprioriCovariance( i,i ) = varianceVector.at( i );
        }
        std::cout << "a priori covariance matrix:" << std::endl << aprioriCovariance << std::endl;
        inverseOfAprioriCovariance = aprioriCovariance.inverse();
    } else{
        inverseOfAprioriCovariance = Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));
    }


    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                inverseOfAprioriCovariance,
                initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, false, false, true );
//    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->manuallySetObservationWeightMatrix(observationWeightsAndTimes);

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( maximumNumberOfIterations ) );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "provide output..." << std::endl;

    // propagate according to integration history. earlier result is separated here as the times are needed on their own.
    std::vector<double> fullStateHistoryTimes;
    std::map<double, Eigen::VectorXd> propagationHistory = allBodiesPropagationHistory.at( 0 );
    std::map<double, Eigen::VectorXd>::iterator historyit = propagationHistory.begin();

    while (historyit != propagationHistory.end()){
        fullStateHistoryTimes.push_back(historyit->first);
        historyit++;
    }


    // Propagate covariance matrix
    std::cout<<"propagating covariance matrix..."<<std::endl;

    Eigen::MatrixXd initialCovarianceMatrix = podOutput->getUnnormalizedCovarianceMatrix( );

    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(propagatedCovariance, initialCovarianceMatrix,
                        orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                        fullStateHistoryTimes);
                        //baseTimeList
                        //60.0, initialEphemerisTime + 3600.0, finalEphemerisTime - 3600.0 );

    std::map<double, double> trueAnomalyMap;
    std::map<double, Eigen::Vector6d> propagatedErrorUsingCovMatrix;
    std::map<double, Eigen::Vector6d> propagatedRSWErrorUsingCovMatrix;

    std::map<double, Eigen::MatrixXd>::iterator it = propagatedCovariance.begin();

    while (it != propagatedCovariance.end()){

        Eigen::Vector6d currentCartesianState = propagationHistory.at(it->first);

        // save true anomaly of MESSENGER around Mercury at this time step
        Eigen::Vector6d currentKeplerianState = convertCartesianToKeplerianElements(
                    currentCartesianState, sunGravitationalParameter);
        double currentTrueAnomaly = currentKeplerianState(5);

        trueAnomalyMap.insert(std::make_pair( it->first, currentTrueAnomaly ) );

        // save propagated error (using covariance) in Cartesian and RSW frame
        Eigen::MatrixXd currentCovariance = it->second;
        Eigen::Matrix3d currentTransformationToRSW =
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(currentCartesianState);

        Eigen::Vector6d currentSigmaVector = currentCovariance.diagonal( ).cwiseSqrt( );

        Eigen::Vector6d currentSigmaVectorRSW;
        currentSigmaVectorRSW.segment(0,3) =
                ( currentTransformationToRSW * currentCovariance.block(0,0,3,3) * currentTransformationToRSW.transpose( ) ).diagonal( ).cwiseSqrt( );
        currentSigmaVectorRSW.segment(3,3) =
                ( currentTransformationToRSW * currentCovariance.block(3,3,3,3) * currentTransformationToRSW.transpose( ) ).diagonal( ).cwiseSqrt( );

        propagatedErrorUsingCovMatrix.insert(std::make_pair( it->first, currentSigmaVector ) );
        propagatedRSWErrorUsingCovMatrix.insert(std::make_pair( it->first, currentSigmaVectorRSW ) );

        it++;
    }


    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    Eigen::VectorXd formalError = podOutput->getFormalErrorVector( );

    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;
    std::cout << "True to form estimation error ratio is: " << std::endl <<
                 ( podOutput->getFormalErrorVector( ).cwiseQuotient( estimationError ) ).transpose( ) << std::endl;

    // Save data in files
    std::cout<<"writing output to files..."<<std::endl;


    input_output::writeDataMapToTextFile( propagatedCovariance,
                                          "PropagatedCovariance.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedErrorUsingCovMatrix,
                                          "propagatedErrorUsingCovMatrix.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedRSWErrorUsingCovMatrix,
                                          "propagatedRSWErrorUsingCovMatrix.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( saveSatelliteError,
                                          "vehicleErrorBasedOnTrueAnomaly.dat",
                                          outputSubFolder );

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "EstimationInformationMatrix.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( truthParameters,
                                     "TruthParameters.dat", 16,
                                     outputSubFolder );


    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "EstimationInformationMatrixNormalization.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "EstimationWeightsDiagonal.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "EstimationResiduals.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "EstimationCorrelations.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "ResidualHistory.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "ParameterHistory.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "ObservationMeasurements.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationTimes.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "ObservationLinkEnds.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationObservableTypes.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "ObservationMeasurements.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( estimationError,
                                     "ObservationTrueEstimationError.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getFormalErrorVector(),
                                     "ObservationFormalEstimationError.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podInput->getInverseOfAprioriCovariance( ),
                                     "InverseAPrioriCovariance.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_,
                                     "InverseNormalizedCovariance.dat", 16,
                                     outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ),
                                     "UnnormalizedCovariance.dat", 16,
                                     outputSubFolder );

    std::cout << "done!" << std::endl;

// Final statement.
// The exit code EXIT_SUCCESS indicates that the program was successfully executed.
return EXIT_SUCCESS;

}
