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

    std::cout<<"importing namespaces..."<<std::endl;

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

    // Load spice kernels.
    std::cout<<"loading spice kernels..."<<std::endl;

    std::string kernelsPath = input_output::getSpiceKernelPath( );

    std::vector< std::string > customKernels;
    customKernels.push_back( kernelsPath + "tudat_merged_spk_kernel_thesis4.bsp" );
    spice_interface::loadStandardSpiceKernels( customKernels );
//            spice_interface::loadStandardSpiceKernels( );


    ////////////////////////////
    //// APPLICATION INPUTS ////
    ////////////////////////////

    std::cout<<"processing inputs"<<std::endl;

    // Acceleration settings
    const bool calculateDeSitterCorrection = false; //no partials implemented
    const bool estimateJ2Amplitude = true;
    const bool estimateJ2Period = false;
    const bool estimateJ2Phase = false;
    const bool estimateJ4Amplitude = true;
    const bool estimateJ4Period = false;
    const bool estimateJ4Phase = false;

    // Parameter estimation settings
    const unsigned int maximumNumberOfIterations = 10;
    const double sigmaPosition = 1000.0; //educated guess
    const double sigmaVelocity = 1.0; //educated guess
    const bool useMultipleMercuryArcs = false;
    const bool reintegrateVariationalEquations = true; //warning: when setting this false and estimating eta with a true value of 0, eta won't improve as the partials are not calculated after perturbing the initial state
    const bool ignoreNordtvedtConstraintInEstimation = false;
    const bool includeSpacecraftPositionError = true;
    const bool includeLightTimeCorrections = false;
    const double observationReductionFactor = 1.0; //decreases amount of observations by a factor n to speed up the algorithm considerably
    const bool testCeres = false; // to perform tests on the asteroid consider covariance
    const bool testWithoutFlybys = true;


    // integrator settings
    const double initialTimeStep = 3600.0/2.0;
    const double minimumStepSize = 3600.0/2.0;
    const double maximumStepSize = 3600.0/2.0;
    const double relativeErrorTolerence = 1.0;
    const double absoluteErrorTolerence = 1.0;
    const unsigned int minimumOrder = 12;
    const unsigned int maximumOrder = 12;

    // Other planetary parameters, currently not included in json
    const double mercuryGravitationalParameter = (2.2031870798779644e+04)*(1E9); //m3/s2, from https://pgda.gsfc.nasa.gov/products/71
    const double sunRadius = 695.7E6; //m, from nasa fact sheet

    // Solar cycle
    Eigen::Vector6i solarMinimumDatetime;
    solarMinimumDatetime << 2008, 12, 15, 0, 0, 0;
    const double solarMinimumEpoch = secondsAfterJ2000(solarMinimumDatetime);
    const double solarCycleDuration = 11.0*physical_constants::JULIAN_YEAR;
    const double solarDay = 25.38*physical_constants::JULIAN_DAY; //carrington sidereal rotation period

    // output settings
    int onlyEveryXthValue = 20;

    // conversion for the asteroids
    double convertAsteroidGMtoSI = 1E-18
            * physical_constants::ASTRONOMICAL_UNIT*physical_constants::ASTRONOMICAL_UNIT*physical_constants::ASTRONOMICAL_UNIT
            / (physical_constants::JULIAN_DAY*physical_constants::JULIAN_DAY);



    ////////////////////////
    //// MISSION INPUTS ////
    ////////////////////////


    // parameter inputs
    std::vector< std::string > filenames;
    filenames.push_back("inputs_multiplemissions_Fienga2019.json");
//    filenames.push_back("inputs_multiplemissions_Park2017.json");
//    filenames.push_back("inputs_multiplemissions_Antia2008.json");

    // scenario pairs
    std::vector< std::pair<int,int> > scenarioPairs;
    scenarioPairs.push_back(std::make_pair(1,1));
//    scenarioPairs.push_back(std::make_pair(2,2));
//    scenarioPairs.push_back(std::make_pair(3,3));
//    scenarioPairs.push_back(std::make_pair(4,4));
//    scenarioPairs.push_back(std::make_pair(3,1));
//    scenarioPairs.push_back(std::make_pair(1,3));


    for (unsigned int f = 0; f<filenames.size(); f++){
        for (unsigned int s = 0; s<scenarioPairs.size(); s++){

            std::string input_filename = filenames.at(f);
            std::cout<<"---- RUNNING SIMULATION FOR INPUTS WITH FILENAME: "<<input_filename<<" ----"<<std::endl;

            std::string json_directory = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/";
            std::ifstream json_file(json_directory + input_filename);
            json json_input;
            json_file >> json_input;
    //        std::cout<<"input values imported: "<<json_input<<std::endl;

            // Scenario
            const int realityScenario = scenarioPairs.at(s).first;
            const int estimationScenario = scenarioPairs.at(s).second;
            std::cout<<"---- AND SCENARIOS:"<<std::endl;
            std::cout<<"---- reality scenario: "<<printScenario(realityScenario)<<std::endl;
            std::cout<<"---- estimation scenario: "<<printScenario(estimationScenario)<<std::endl;

            unsigned int maximumDegreeSunGravitationalMomentInReality = 2;
            if (realityScenario == 2 || realityScenario == 4 || realityScenario == 6){ maximumDegreeSunGravitationalMomentInReality = 4; }
            unsigned int maximumDegreeSunGravitationalMomentInEstimation = 2;
            if (estimationScenario == 2 || estimationScenario == 4 || estimationScenario == 6){ maximumDegreeSunGravitationalMomentInEstimation = 4; }

            bool includeTimeVaryingGravitationalMomentsSunInReality = false;
            if (realityScenario > 2){includeTimeVaryingGravitationalMomentsSunInReality = true;}
            bool includeTimeVaryingGravitationalMomentsSunInEstimation = false;
            if (estimationScenario > 2){includeTimeVaryingGravitationalMomentsSunInEstimation = true;}

            unsigned int maximumDegreeSunGravitationalMomentVariationInReality = 2;
            if (realityScenario == 4 || realityScenario == 6){maximumDegreeSunGravitationalMomentVariationInReality = 4;}
            unsigned int maximumDegreeSunGravitationalMomentVariationInEstimation = 2;
            if (estimationScenario == 4 || estimationScenario == 6){maximumDegreeSunGravitationalMomentVariationInEstimation = 4;}

            // Location of simulation output
            std::string subFolderName = json_input["outputSubFolderName"];
            std::string outputSubFolderName = subFolderName
                    + "_reality" + std::to_string(realityScenario)
                    + "_estimation" + std::to_string(estimationScenario);
            std::string outputSubFolder = getOutputPath( ) + outputSubFolderName;
            if (useMultipleMercuryArcs) {outputSubFolder = outputSubFolder + "_multiarc";}
            if (testCeres){
                outputSubFolder += "_testCeres";
            }
            if (testWithoutFlybys){
                outputSubFolder += "_testWithoutFlybys";
            }

            // Acceleration settings
            const bool calculateSchwarzschildCorrection = json_input["calculateSchwarzschildCorrection"];
            const bool calculateLenseThirringCorrection = json_input["calculateLenseThirringCorrection"];
            const bool includeSEPViolationAcceleration = json_input["includeSEPViolationAcceleration"];
            const bool includeTVGPAcceleration = json_input["includeTVGPAcceleration"];

            // Retrieve input parameters including uncertainties and apriori values
            const double sunAngularMomentum = json_input["sunAngularMomentum"];
            const double sunGravitationalParameter = json_input["sunGravitationalParameter"];
            const double timeVaryingGravitationalParameter = json_input["timeVaryingGravitationalParameter"];
            const double sigmaSunAngularMomentum = json_input["sigma_S_Sun"];
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

            // Parameter settings
            const bool gammaIsAConsiderParameter = json_input["gammaIsAConsiderParameter"];
            const bool nordtvedtConstraintTrueOrFalse = json_input["useNordtvedtConstraint"];

            const bool estimateSunAngularMomentum = json_input["estimateSunAngularMomentum"];
            bool considerSunAngularMomentum = false;
            if (estimateSunAngularMomentum == false){ considerSunAngularMomentum = true; };
            const bool estimatePPNalphas = json_input["estimatePPNalphas"];

            bool considerPPNalphas = false;
            if (estimatePPNalphas == false){ considerPPNalphas = true; };


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
                    if (testWithoutFlybys == false){
                        flybyList.push_back(secondsAfterJ2000(currentFlyby));
                    }
                    currentFlyby.clear();
                }
                flybyListVector.push_back(flybyList);


                // observation Noise
                noiseAtMinAngleVector.push_back(json_input_mission["noiseAtMinAngle"]);
                noiseAtMaxAngleVector.push_back(json_input_mission["noiseAtMaxAngle"]);
                maxMSEAngleDegVector.push_back(json_input_mission["maxMSEAngleDeg"]);

            }


            // settings time variable gravitational moments
            double J2amplitudeUnnormalized = json_input["sunJ2_A"];
            double J2amplitude = J2amplitudeUnnormalized / calculateLegendreGeodesyNormalizationFactor(2,0);
        //        double J2amplitude = 0.005E-7 / calculateLegendreGeodesyNormalizationFactor(2,0); // currently from Antia et al 2008
            double J2period = solarCycleDuration;
            double J2phase = phaseAccordingToSolarMinimum(solarMinimumEpoch, J2period);

            double J4amplitudeUnnormalized = json_input["sunJ4_A"];
            double J4amplitude = J4amplitudeUnnormalized / calculateLegendreGeodesyNormalizationFactor(4,0);
            double J4period = solarCycleDuration;
            double J4phase = phaseAccordingToSolarMinimum(solarMinimumEpoch, J4period); //antiphase?

            double J2phaseInEstimation = J2phase;
            double J4phaseInEstimation = J4phase;
            if (realityScenario > 4){
                J2phase += mathematical_constants::PI;
                J4phase += mathematical_constants::PI;
            }
            if (estimationScenario > 4){
                J2phaseInEstimation += mathematical_constants::PI;
                J4phaseInEstimation += mathematical_constants::PI;
            }


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


            // some checks to prevent incompatible inputs
            if ((calculateLenseThirringCorrection == false) && (estimateSunAngularMomentum)){
                std::runtime_error("cannot estimate sun angular momentum when lense thirring acceleration is not calculated");
            }
            if ((estimateSunAngularMomentum) && (considerSunAngularMomentum)){
                std::runtime_error("ppn alphas cannot both be estimatable parameters and consider parameters");
            }
            if ((estimatePPNalphas) && (considerPPNalphas)){
                std::runtime_error("ppn alphas cannot both be estimatable parameters and consider parameters");
            }

            if ((includeTimeVaryingGravitationalMomentsSunInReality) && (maximumDegreeSunGravitationalMomentVariationInReality < 2)){
                std::runtime_error("varying graviational moments cannot be included in reality, check settings");
            }
            if ((includeTimeVaryingGravitationalMomentsSunInEstimation) && (maximumDegreeSunGravitationalMomentVariationInEstimation < 2)){
                std::runtime_error("varying graviational moments cannot be included in estimation, check settings");
            }



            /////////////////////
            //// ENVIRONMENT ////
            /////////////////////

            std::cout << "building environment..." << std::endl;



            std::vector< std::string > bodyNames;
            bodyNames.push_back("Sun");
            bodyNames.push_back("Mercury");
            bodyNames.push_back("Venus");
            bodyNames.push_back("Earth");
            bodyNames.push_back("Mars");
            bodyNames.push_back("Jupiter");
            bodyNames.push_back("Saturn");
            bodyNames.push_back("Uranus");
            bodyNames.push_back("Neptune");
            bodyNames.push_back("Moon");

            if (testCeres){
                spice_interface::loadSpiceKernelInTudat( kernelsPath + "codes_300ast_20100725.tf" );
                spice_interface::loadSpiceKernelInTudat( kernelsPath + "codes_300ast_20100725.bsp" );
                bodyNames.push_back("2000001");
            }

            const unsigned int totalNumberOfBodies = bodyNames.size();

            // load SPICE settings
            double initialSimulationTime;
            if (testWithoutFlybys){
                initialSimulationTime = *std::min_element(observationInitialTimeVector.begin(), observationInitialTimeVector.end());
            } else{
                initialSimulationTime = *std::min_element(initialTimeVector.begin(), initialTimeVector.end());
            }
            const double finalSimulationTime = *std::max_element(finalTimeVector.begin(), finalTimeVector.end());
            const double buffer = maximumOrder*maximumStepSize; //see Tudat libraries 1.1.3.
            std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;

            // Default body settings
            bodySettings = getDefaultBodySettings( bodyNames,
                                                   initialSimulationTime - buffer,
                                                   finalSimulationTime + buffer,
                                                   initialTimeStep);

//            // Simplification, for outer solar sytem bodies, take tabulated values
            bodySettings["Jupiter"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Saturn"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Uranus"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Neptune"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );

            if (testCeres){
                // set ephemeris settings to tabulated with 1 hour time step
                bodySettings[ "2000001" ]->ephemerisSettings
                        = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                            initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );

                // set gravitational parameter as given in INPOP19a
                bodySettings[ "2000001" ]->gravityFieldSettings
                        = std::make_shared< CentralGravityFieldSettings >(139643.532 * convertAsteroidGMtoSI);
            }


            std::cout << "setting custom settings Sun..." << std::endl;

            // Nominal values spherical harmonics Sun
            Eigen::MatrixXd normalizedSineCoefficients
                    = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMomentInReality + 1,
                                            maximumDegreeSunGravitationalMomentInReality + 1);
            Eigen::MatrixXd normalizedCosineCoefficients
                    = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMomentInReality + 1,
                                            maximumDegreeSunGravitationalMomentInReality + 1);

            if (normalizedCosineCoefficients.col(0).size() != valuesSunGravitationalMoments.size()){
                std::runtime_error("error: vector sizes incompatible");
            }

            normalizedCosineCoefficients.col(0) = valuesSunGravitationalMoments.transpose();

            bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        sunGravitationalParameter, sunRadius,
                        normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );

            // angular momentum
    //        const Eigen::Vector3d sunAngularMomentumVectorInSunFrame(0.0, 0.0, sunAngularMomentum);
    //        const Eigen::Vector3d sunAngularMomentumVectorPerUnitMassInSunFrame =
    //                sunAngularMomentumVectorInSunFrame /
    //                (sunGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT);
    //        const double sunAngularMomentumPerUnitMass =
    //                sunAngularMomentum / (sunGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT);


            bodySettings[ "Sun" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Sun",
                        spice_interface::computeRotationQuaternionBetweenFrames("ECLIPJ2000", "IAU_Sun", initialSimulationTime ),
                        initialSimulationTime, 2.0 * mathematical_constants::PI / solarDay,
                        sunAngularMomentum );



            // Time varying spherical harmonics coefficients Sun
            if (includeTimeVaryingGravitationalMomentsSunInReality || estimationScenario > 2){

                // J2
                std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

                if ((includeTimeVaryingGravitationalMomentsSunInReality == false)
                    && (estimationScenario > 2)){ // if in reality there is no time variable J2 but we want to estimate it later, include it anyway in reality with amplitude practically 0, such that environment/acceleration models are properly loaded
                    relativity::variableJ2Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(2)/1.0E5);
                } else{
                    relativity::variableJ2Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(2));
                }
                relativity::variableJ2Interface->setPeriod(periodsSunGravitationalMomentsVariation(2));
                relativity::variableJ2Interface->setPhase(phasesSunGravitationalMomentsVariation(2));

                std::function< double() > amplitudeFunctionJ2 =
                        std::bind( &relativity::VariableJ2Interface::getAmplitude, relativity::variableJ2Interface );
                std::function< double() > periodFunctionJ2 =
                        std::bind( &relativity::VariableJ2Interface::getPeriod, relativity::variableJ2Interface );
                std::function< double() > phaseFunctionJ2 =
                        std::bind( &relativity::VariableJ2Interface::getPhase, relativity::variableJ2Interface );


                std::function< std::map< double, Eigen::MatrixXd >( ) > tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ2 =
                        std::bind(tabulatedSphericalHarmonicsCoefficientCorrections,
                                  initialSimulationTime, finalSimulationTime,
                                  amplitudeFunctionJ2, periodFunctionJ2, phaseFunctionJ2);

                std::map< double, Eigen::MatrixXd > sineCoefficientCorrectionsJ2 =
                        zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

                const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettingsJ2 =
                        std::make_shared< TabulatedGravityFieldVariationSettingsWithCosineFunction >(
                            tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ2,
                            sineCoefficientCorrectionsJ2,
                            2, 0, interpolatorSettings,
                            initialSimulationTime, finalSimulationTime, 3600.0);

                gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettingsJ2 );

                //J4
                if (maximumDegreeSunGravitationalMomentVariationInReality >= 4 || estimationScenario == 4 || estimationScenario == 6){

                    if (includeTimeVaryingGravitationalMomentsSunInReality == false){ // if in reality there is no time variable J4 but we want to estimate it later, include it anyway in reality with amplitude practically 0, such that environment/acceleration models are properly loaded
                        relativity::variableJ4Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(4)/1.0E5);
                    } else{
                         relativity::variableJ4Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(4));
                    }
                    relativity::variableJ4Interface->setPeriod(periodsSunGravitationalMomentsVariation(4));
                    relativity::variableJ4Interface->setPhase(phasesSunGravitationalMomentsVariation(4));

                    std::function< double() > amplitudeFunctionJ4 =
                            std::bind( &relativity::VariableJ4Interface::getAmplitude, relativity::variableJ4Interface );
                    std::function< double() > periodFunctionJ4 =
                            std::bind( &relativity::VariableJ4Interface::getPeriod, relativity::variableJ4Interface );
                    std::function< double() > phaseFunctionJ4 =
                            std::bind( &relativity::VariableJ4Interface::getPhase, relativity::variableJ4Interface );


                    std::function< std::map< double, Eigen::MatrixXd >( ) > tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ4 =
                            std::bind(tabulatedSphericalHarmonicsCoefficientCorrections,
                                      initialSimulationTime, finalSimulationTime,
                                      amplitudeFunctionJ4, periodFunctionJ4, phaseFunctionJ4);

                    std::map< double, Eigen::MatrixXd > sineCoefficientCorrectionsJ4 =
                            zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

                    const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettingsJ4 =
                            std::make_shared< TabulatedGravityFieldVariationSettingsWithCosineFunction >(
                                tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ4,
                                sineCoefficientCorrectionsJ4,
                                2, 0, interpolatorSettings,
                                initialSimulationTime, finalSimulationTime, 3600.0);

                    gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettingsJ4 );
                }

                bodySettings[ "Sun" ]->gravityFieldVariationSettings = gravityFieldVariationSettings;
            }

            std::cout << "creating environment:" << std::endl;

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
            bodiesToPropagate.push_back("Mercury");
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

            // Setting all PPN parameters to GR values (if multiple simulations are done in sequence, estmiated parameters are not reset automatically)
            relativity::ppnParameterSet->setParameterGamma(1.0);
            relativity::ppnParameterSet->setParameterBeta(1.0);
            relativity::ppnParameterSet->setParameterAlpha1(0.0);
            relativity::ppnParameterSet->setParameterAlpha2(0.0);
            relativity::ppnParameterSet->setNordtvedtParameter(0.0);

    //        std::shared_ptr< ephemerides::RotationalEphemeris > sunRotationalEphemeris
    //                = bodyMap.at("Sun")->getRotationalEphemeris();
    //        Eigen::Matrix3d localSunFrameToGlobalFrame =
    //                sunRotationalEphemeris->getRotationToTargetFrame(initialSimulationTime).normalized().toRotationMatrix();


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
                                        std::make_shared< SphericalHarmonicAccelerationSettings > (maximumDegreeSunGravitationalMomentInReality,0));

                            if (calculateSchwarzschildCorrection || calculateLenseThirringCorrection || calculateDeSitterCorrection){
                                currentAccelerations[ bodyNames.at( j ) ].push_back(
                                            std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                calculateSchwarzschildCorrection,
                                                calculateLenseThirringCorrection,
                                                calculateDeSitterCorrection));
                            }

                            if (includeSEPViolationAcceleration == true){
                                currentAccelerations[ bodyNames.at( j ) ].push_back(
                                            std::make_shared< SEPViolationAccelerationSettings >(
                                                bodyNames, nordtvedtConstraintTrueOrFalse, ignoreNordtvedtConstraintInEstimation));
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
                        "Mercury", "Sun", maximumDegreeSunGravitationalMomentInReality, 0 ) );

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

            Eigen::Matrix<long double, Eigen::Dynamic, 1> systemInitialState;
            std::shared_ptr< PropagatorSettings< long double > > propagatorSettings;
//            std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;

//            if (useMultipleMercuryArcs){

//                std::vector< Eigen::VectorXd > arcInitialStates;

//                for (unsigned int m=0; m<numberOfMissions; m++){

//                    Eigen::VectorXd currentArcInitialState = getBodyCartesianStateAtEpoch(
//                                bodiesToPropagate.at(0), centralBodies.at(0),
//                                "ECLIPJ2000", "None", initialTimeVector.at(m) );

//                    std::cout<<"initial state: "<<currentArcInitialState.transpose()<<std::endl;
//                    arcInitialStates.push_back( currentArcInitialState );

//                    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
//                          std::make_shared< propagators::PropagationTimeTerminationSettings >(
//                                finalTimeVector.at(m), true );

//                    propagatorSettingsList.push_back(
//                                std::make_shared< TranslationalStatePropagatorSettings< long double > >(
//                                    centralBodies, accelerationModelMap, bodiesToPropagate,
//                                    currentArcInitialState, terminationSettings, cowell, dependentVariablesToSave ) );
//                }

//                // Create propagator settings
//                propagatorSettings = std::make_shared< MultiArcPropagatorSettings< long double > >( propagatorSettingsList );

//            } else{

                // Get initial state of bodies to be propagated
                systemInitialState = getInitialStatesOfBodies(
                            bodiesToPropagate, centralBodies, bodyMap, initialSimulationTime )
                                .cast<long double>();

                // Define propagator settings.
                std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                      std::make_shared< propagators::PropagationTimeTerminationSettings >( finalSimulationTime, true );

                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< long double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate,
                          systemInitialState, terminationSettings, cowell, dependentVariablesToSave);

//            }


            std::shared_ptr< AdamsBashforthMoultonSettings< double > > integratorSettings =
                    std::make_shared< AdamsBashforthMoultonSettings< double > > (
                        initialSimulationTime, initialTimeStep,
                        minimumStepSize, maximumStepSize,
                        relativeErrorTolerence, absoluteErrorTolerence,
                        minimumOrder, maximumOrder);
            std::shared_ptr< AdamsBashforthMoultonSettings< double > > backwardIntegratorSettings =
                    std::make_shared< AdamsBashforthMoultonSettings< double > > (
                        finalSimulationTime, -1.0*initialTimeStep,
                        -1.0*minimumStepSize, -1.0*maximumStepSize,
                        relativeErrorTolerence, absoluteErrorTolerence,
                        minimumOrder, maximumOrder);

        //    std::shared_ptr< IntegratorSettings< > > integratorSettings =
        //            std::make_shared< IntegratorSettings< > >
        //            ( rungeKutta4, initialSimulationTime, initialTimeStep );

//            std::shared_ptr< IntegratorSettings< > > integratorSettings =
//                    std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
//                        initialSimulationTime, initialTimeStep,
//                        RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
//                        minimumStepSize, maximumStepSize,
//                        relativeErrorTolerence, absoluteErrorTolerence);
//            std::shared_ptr< IntegratorSettings< > > backwardIntegratorSettings =
//                    std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
//                        finalSimulationTime, -1.0*initialTimeStep,
//                        RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
//                        -1.0*minimumStepSize, -1.0*maximumStepSize,
//                        relativeErrorTolerence, absoluteErrorTolerence);



            ////////////////////////////
            //// DYNAMICS SIMULATOR ////
            ////////////////////////////

            std::cout << "running dynamics simulator..." << std::endl;

            std::vector< std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > > allBodiesPropagationHistory;
            std::vector< std::map< double, Eigen::VectorXd > > spiceStatesAtPropagationTimes;
            std::vector< Eigen::VectorXd > arcFinalStates;
            allBodiesPropagationHistory.resize( bodiesToPropagate.size() );
            spiceStatesAtPropagationTimes.resize( bodiesToPropagate.size() );

            std::map< double, Eigen::VectorXd > dependentVariablesHistory;

//            if (useMultipleMercuryArcs){


//                MultiArcDynamicsSimulator <long double> dynamicsSimulator (bodyMap,
//                                                                integratorSettings,
//                                                                propagatorSettings,
//                                                                initialTimeVector,
//                                                                true,false,true);

//                std::vector< std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//                std::vector< std::map< double, Eigen::VectorXd > > dependentVariablesResult = dynamicsSimulator.getDependentVariableHistory();

//                for( unsigned int m = 0; m < numberOfMissions; m++ ){

//                    std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > intermediateIntegrationResult = integrationResult.at( m );

//                    // Retrieve numerically integrated state for each body.
//                    for( std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
//                         stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
//                    {
//                        for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
//                        {
//                            allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
//                            spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
//                                    getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);

//                        }
//                    }

//                    arcFinalStates.push_back(allBodiesPropagationHistory[0].at(finalTimeVector.at(m)));

//                    std::map< double, Eigen::VectorXd > intermediateDependentVariablesResult = dependentVariablesResult.at( m );

//                    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateDependentVariablesResult.begin( );
//                         stateIterator != intermediateDependentVariablesResult.end( ); stateIterator++ )
//                    {
//                        dependentVariablesHistory[ stateIterator->first ] = stateIterator->second;
//                    }
//                }

//            } else{

                SingleArcDynamicsSimulator <long double> dynamicsSimulator (bodyMap,integratorSettings,propagatorSettings,true,false,true);
                std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory();

                // Retrieve numerically integrated state for each body.
                for( std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> >::const_iterator stateIterator = integrationResult.begin( );
                     stateIterator != integrationResult.end( ); stateIterator++ )
                {
                    for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
                    {
                        allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
                        spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
                                getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);

                    }
                }

//            }

            for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
            {
                // Write propagation history to file.
                input_output::writeDataMapToTextFile(
                            onlyEveryXthValueFromDataMap( allBodiesPropagationHistory[ i ], onlyEveryXthValue),
                            "StatePropagationHistory" + bodiesToPropagate.at( i ) + ".dat",
                            outputSubFolder,
                            "",
                            std::numeric_limits< long double >::digits10,
                            std::numeric_limits< long double >::digits10,
                            "," );

                input_output::writeDataMapToTextFile(
                            onlyEveryXthValueFromDataMap( spiceStatesAtPropagationTimes[ i ], onlyEveryXthValue),
                            "spiceStatesAtPropagationTimes" + bodiesToPropagate.at( i ) + ".dat",
                            outputSubFolder,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );
            }

            // Write dependent variables history to file.
            input_output::writeDataMapToTextFile(
                        onlyEveryXthValueFromDataMap(dependentVariablesHistory, onlyEveryXthValue),
                        "DependentVariablesHistoryReality.dat",
                        outputSubFolder,
                        "",
                        std::numeric_limits< double >::digits10,
                        std::numeric_limits< double >::digits10,
                        "," );

            // clear memory
    //        allBodiesPropagationHistory.clear();
            dependentVariablesHistory.clear();


            ////////////////////////////
            //// BACKWARD PROPAGATION ////
            ////////////////////////////

            std::cout << "performing backward propagation..." << std::endl;

//            std::shared_ptr< AdamsBashforthMoultonSettings< double > > backwardIntegratorSettings =
//                    std::make_shared< AdamsBashforthMoultonSettings< double > > (
//                        finalSimulationTime, -1.0*initialTimeStep,
//                        -1.0*minimumStepSize, -1.0*maximumStepSize,
//                        relativeErrorTolerence, absoluteErrorTolerence,
//                        minimumOrder, maximumOrder);

            std::shared_ptr< PropagationTimeTerminationSettings > backwardTerminationSettings;
            std::shared_ptr< PropagatorSettings< long double > > backwardPropagatorSettings;

            std::vector< std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > > allBodiesBackwardPropagationHistory;
            allBodiesBackwardPropagationHistory.resize( bodiesToPropagate.size() );

//            if (useMultipleMercuryArcs){

//                std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > backwardPropagatorSettingsList;

//                for (unsigned int m=0; m<numberOfMissions; m++){

//                    backwardTerminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >(
//                                initialTimeVector.at(m), true );

//                    backwardPropagatorSettingsList.push_back(
//                                std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                                    centralBodies, accelerationModelMap, bodiesToPropagate,
//                                    arcFinalStates.at(m), backwardTerminationSettings, cowell ) );
//                }

//                // Create propagator settings
//                backwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< double > >( backwardPropagatorSettingsList );


//                MultiArcDynamicsSimulator <> backwardDynamicsSimulator (bodyMap,
//                                                                backwardIntegratorSettings,
//                                                                backwardPropagatorSettings,
//                                                                finalTimeVector,
//                                                                true,false,true);

//                std::vector< std::map< double, Eigen::VectorXd > > backwardIntegrationResult = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

//                for( unsigned int m = 0; m < numberOfMissions; m++ ){

//                    std::map< double, Eigen::VectorXd > intermediateIntegrationResult = backwardIntegrationResult.at( m );

//                    // Retrieve numerically integrated state for each body.
//                    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
//                         stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
//                    {
//                        for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
//                        {
//                            allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
//                            spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
//                                    getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);

//                        }
//                    }

//                }

//            } else{

                backwardTerminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( initialSimulationTime, true );

                Eigen::Matrix<long double, Eigen::Dynamic, 1> systemFinalState = integrationResult.at(finalSimulationTime);

                backwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< long double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate,
                          systemFinalState, backwardTerminationSettings);


                SingleArcDynamicsSimulator <long double> backwardDynamicsSimulator
                        (bodyMap,backwardIntegratorSettings,backwardPropagatorSettings,true,false,true);
                std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > backwardIntegrationResult
                        = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

                // Retrieve numerically integrated state for each body.
                for( std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> >::const_iterator stateIterator = backwardIntegrationResult.begin( );
                     stateIterator != backwardIntegrationResult.end( ); stateIterator++ )
                {
                    for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
                    {
                        allBodiesBackwardPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );

                    }
                }

//            }

            for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
            {
                // Write propagation history to file.
                input_output::writeDataMapToTextFile(
                            onlyEveryXthValueFromDataMap( allBodiesBackwardPropagationHistory[ i ], onlyEveryXthValue),
                            "StatePropagationHistory" + bodiesToPropagate.at( i ) + "Backwards.dat",
                            outputSubFolder,
                            "",
                            std::numeric_limits< long double >::digits10,
                            std::numeric_limits< long double >::digits10,
                            "," );

            }

            // clear memory
            allBodiesBackwardPropagationHistory.clear();



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

            std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
            if (includeLightTimeCorrections){
                std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
                lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                           lightTimePerturbingBodies ) );
            }

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

            // In order to get the partial derivatives of consider parameters wrt observations, run an estimation of the consider parameters
            std::vector< std::shared_ptr < EstimatableParameterSettings > > parameterNames;
            std::vector< double > varianceVector;

            // Add bodies that will be propagated to the parameters to be estimated
//            if (useMultipleMercuryArcs){
//                Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * initialTimeVector.size( ) );
//                for( unsigned int m = 0; m < numberOfMissions; m++ ){
//                    systemInitialState.segment( m * 6, 6 ) = propagatorSettingsList.at( m )->getInitialStates( );
//                    for( unsigned int i = 0; i < 3; i++ ){ varianceVector.push_back(sigmaPosition*sigmaPosition); }
//                    for( unsigned int i = 3; i < 6; i++ ){ varianceVector.push_back(sigmaVelocity*sigmaVelocity); }
//                }
//                parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
//                                              "Mercury", systemInitialState, initialTimeVector, "SSB" ) );
//            }
//            else{
                for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ){
                    int j = 6*i;
                    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< long double > >(
                                                  bodiesToPropagate[i], systemInitialState.segment(j,6), centralBodies[i] ) );
                    for( unsigned int i = 0; i < 3; i++ ){ varianceVector.push_back(sigmaPosition*sigmaPosition); }
                    for( unsigned int i = 3; i < 6; i++ ){ varianceVector.push_back(sigmaVelocity*sigmaVelocity); }
                }
//            }

            if (testCeres){
                // set GM of the asteroid as estimatable parameter, get apriori sigma from INPOP19a
                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("2000001", gravitational_parameter));
                double asteroidSigma = 340.331 * convertAsteroidGMtoSI;
                varianceVector.push_back(asteroidSigma*asteroidSigma);

                std::cout<<"GM Ceres uncertainty [m3/s2]: "<<340.331 * convertAsteroidGMtoSI<<std::endl;
            }


            bool gammaIsEstimated = false;
            bool betaIsEstimated = false;
            bool nordtvedtParameterIsEstimated = false;

            // relativistic parameters
            if (calculateSchwarzschildCorrection == true
                || includeSEPViolationAcceleration == true){

                if ( gammaIsAConsiderParameter == false){
                    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                             ("global_metric", ppn_parameter_gamma ) );
                    varianceVector.push_back(sigmaGamma*sigmaGamma);
                    gammaIsEstimated = true;
                }


                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_parameter_beta ) );
                varianceVector.push_back(sigmaBeta*sigmaBeta);
                betaIsEstimated = true;
            }

            // Nordtvedt parameter
            int nordtvedtParameterIndex = -1;
            if (includeSEPViolationAcceleration &&
                    (nordtvedtConstraintTrueOrFalse == false || ignoreNordtvedtConstraintInEstimation == false)){
                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_nordtvedt_parameter ) );
                varianceVector.push_back(sigmaNordtvedt*sigmaNordtvedt);
                nordtvedtParameterIsEstimated = true;
                nordtvedtParameterIndex = varianceVector.size()-1;
            }


            if (estimatePPNalphas == true){
                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_parameter_alpha1 ) );
                varianceVector.push_back(sigmaAlpha1*sigmaAlpha1);

                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_parameter_alpha2 ) );
                varianceVector.push_back(sigmaAlpha2*sigmaAlpha2);
            }

            // angular momentum
            if (calculateLenseThirringCorrection && estimateSunAngularMomentum){
                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("Sun", angular_momentum));
                varianceVector.push_back(sigmaSunAngularMomentum*sigmaSunAngularMomentum);
            }


            // time varying gravitational parameter
            if (includeTVGPAcceleration == true){
                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", time_varying_gravitational_parameter));
                varianceVector.push_back(sigmaTVGP*sigmaTVGP);
            }

            // gravitational parameter Sun
            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                     ("Sun", gravitational_parameter));
            varianceVector.push_back(sigmaSunGM*sigmaSunGM);

            if (includeTimeVaryingGravitationalMomentsSunInEstimation){

                // time varying J2 Sun
                if (estimationScenario >= 2){
                    double variableJ2parametersVariance = 1.0; //sigma taken as a percentage of the mean
                    if (estimateJ2Amplitude){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J2_amplitude));
                        double sigmaJ2Amplitude = variableJ2parametersVariance * relativity::variableJ2Interface->getAmplitude();
                        varianceVector.push_back(sigmaJ2Amplitude*sigmaJ2Amplitude);
                    }
                    if (estimateJ2Period){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J2_period));
                        double sigmaJ2Period = variableJ2parametersVariance * relativity::variableJ2Interface->getPeriod();
                        varianceVector.push_back(sigmaJ2Period*sigmaJ2Period);
                    }
                    if (estimateJ2Phase){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J2_phase));
                        double sigmaJ2Phase = variableJ2parametersVariance * relativity::variableJ2Interface->getPhase();
                        varianceVector.push_back(sigmaJ2Phase*sigmaJ2Phase);
                    }
                }

                // time varying J4 Sun
                if (estimationScenario == 4 || estimationScenario == 6){
                    double variableJ4parametersVariance = 1.0; //sigma taken as a percentage of the mean
                    if (estimateJ4Amplitude){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J4_amplitude));
                        double sigmaJ4Amplitude = variableJ4parametersVariance * relativity::variableJ4Interface->getAmplitude();
                        varianceVector.push_back(sigmaJ4Amplitude*sigmaJ4Amplitude);
                    }
                    if (estimateJ4Period){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J4_period));
                        double sigmaJ4Period = variableJ4parametersVariance * relativity::variableJ4Interface->getPeriod();
                        varianceVector.push_back(sigmaJ4Period*sigmaJ4Period);
                    }
                    if (estimateJ4Phase){
                        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                                 ("Sun", variable_J4_phase));
                        double sigmaJ4Phase = variableJ4parametersVariance * relativity::variableJ4Interface->getPhase();
                        varianceVector.push_back(sigmaJ4Phase*sigmaJ4Phase);
                    }
                }
            }


            // conventional spherical harmonics (always needs to be last estimatable parameter)
            std::vector< std::pair< int, int > > blockIndices;
            for (unsigned int d=2; d<=maximumDegreeSunGravitationalMomentInEstimation; d+=2){
                blockIndices.push_back(std::make_pair(d,0));
                varianceVector.push_back(sigmaValuesSunGravitationalMoments(d)*sigmaValuesSunGravitationalMoments(d));
            }
            parameterNames.push_back(std::make_shared<SphericalHarmonicEstimatableParameterSettings>
                                     (blockIndices,"Sun",spherical_harmonics_cosine_coefficient_block));


            // Create parameters
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
                    createParametersToEstimate( parameterNames, bodyMap, propagatorSettings );

            // Print identifiers and indices of parameters to terminal
            printEstimatableParameterEntries( parametersToEstimate );


            // If gamma, beta, eta, are estimatable parameters and nordtvedt constraint is set to true, enforce it in the estimation
            bool enforceNordtvedtConstraintInEstimation = false;
            if ( (gammaIsEstimated || gammaIsAConsiderParameter )
                && betaIsEstimated && nordtvedtParameterIsEstimated
                && nordtvedtConstraintTrueOrFalse
                && ignoreNordtvedtConstraintInEstimation == false){
                enforceNordtvedtConstraintInEstimation = true;
            }


            ////////////////////////////////////////
            //// RUN ORBITDETERMINATION MANAGER ////
            ////////////////////////////////////////

            std::cout << "Running OD manager..." << std::endl;

            // Create orbit determination object (propagate orbit, create observation models)
            OrbitDeterminationManager< long double, double > orbitDeterminationManager =
                    OrbitDeterminationManager< long double, double >(
                        bodyMap, parametersToEstimate, observationSettingsMap,
                        integratorSettings, propagatorSettings );



            ////////////////////////////////////////
            //// DEFINING OD SIMULATOR SETTINGS ////
            ////////////////////////////////////////

            std::cout << "Defining observation settings..." << std::endl;


            // generate list of observations
            if (observationReductionFactor > 1){
                std::cout<<"WARNING: amount of observations is being reduced by a factor "<<observationReductionFactor<<std::endl;
            }

            std::vector< double > baseTimeList;
            std::vector< std::vector< double > > seperateBaseTimeLists;
            for( unsigned int m = 0; m < numberOfMissions; m++ ){

                std::vector< double > currentMissionBaseTimeList =
                        makeObservationTimeList(observationInitialTimeVector.at(m),
                                                observationFinalTimeVector.at(m),
                                                observationTimeStepVector.at(m)*observationReductionFactor,
                                                trackingArcDurationVector.at(m),
                                                maximumNumberOfTrackingDaysVector.at(m),
                                                unit_conversions::convertDegreesToRadians(maxMSEAngleDegVector.at(m)),
                                                flybyListVector.at(m));

                seperateBaseTimeLists.push_back(currentMissionBaseTimeList);

                baseTimeList.insert( baseTimeList.end(),
                                     currentMissionBaseTimeList.begin(),
                                     currentMissionBaseTimeList.end() );

            }
            int baseTimeListSize = baseTimeList.size();
            std::cout << "Size total observation time list: "<<baseTimeListSize<< std::endl;


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

                for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ ){
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


//            // Create noise functions per observable

//            std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;

//            std::function< double( const double ) > mercuryOrbiterNoiseFunction;
//            mercuryOrbiterNoiseFunction = [noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists](const double time){
//                return noiseSampleBasedOnMSEangleForMultipleMissions
//                        (time, noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists);
//            };

//            noiseFunctions[ one_way_range ] = mercuryOrbiterNoiseFunction;
//            noiseFunctions[ n_way_range ] = mercuryOrbiterNoiseFunction;



            ///////////////////////////////
            //// SIMULATE OBSERVATIONS ////
            ///////////////////////////////

            std::cout << "Simulating observations..." << std::endl;

            // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
            // reference link ends.
            typedef Eigen::Matrix< long double, Eigen::Dynamic, 1 > ObservationVectorType;
            typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
                    SingleObservablePodInputType;
            typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

            // Simulate perfect observations
            PodInputDataType observationsAndTimes = simulateObservations< long double, double >(
                        measurementSimulationInput,
                        orbitDeterminationManager.getObservationSimulators( ) );
                        //,noiseFunctions);

            // similar container, but the "observation" will be the noise value instead of the actual observation
            PodInputDataType observationWeightsAndTimes = observationsAndTimes;

            // add spacecraft initial position error to the observations
            std::map< double, Eigen::Vector3d > interpolatedErrorMatrix;
            if (includeSpacecraftPositionError == true){

                std::cout << "Adding satellite estimation initial position error..." << std::endl;
        //        Eigen::Vector3d constantSatelliteError; constantSatelliteError << 10.0, 10.0, 10.0;

                // get interpolated error maps for lists of vehicle observation times, and merge them

                for (unsigned int m=0; m<numberOfMissions; m++){

                    std::string vehicleErrorFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/error_inputs_"+vehicleVector.at(m)+".txt";
                    Eigen::MatrixXd error_input = input_output::readMatrixFromFile(vehicleErrorFilename, ",");

                    std::map< double, Eigen::Vector3d > interpolatedErrorMatrixPerVehicle =
                            interpolatePositionErrorsBasedOnTrueAnomaly(
                                error_input, seperateBaseTimeLists.at(m),
                                vehicleVector.at(m), mercuryGravitationalParameter);

                    if (interpolatedErrorMatrix.size() == 0){
                        interpolatedErrorMatrix = interpolatedErrorMatrixPerVehicle;
                    } else{
                        std::map< double, Eigen::Vector3d >::iterator it = interpolatedErrorMatrixPerVehicle.begin();
                        while (it != interpolatedErrorMatrixPerVehicle.end()){
                            if ( interpolatedErrorMatrix.find(it->first) == interpolatedErrorMatrix.end() ){
                                interpolatedErrorMatrix.insert(std::make_pair(it->first, it->second));
                            }
                            else{
                                std::cout<<"ERROR at time "<<it->first<<std::endl;
                                std::runtime_error("duplicate observation times!");
                                std::cout<<std::endl;
                            }
                            it++;
                        }
                    }
                    interpolatedErrorMatrixPerVehicle.clear();
                }
                int interpolatedErrorMatrixSize = interpolatedErrorMatrix.size();
                std::cout << "Size interpolated vehicle error map: "<< interpolatedErrorMatrixSize << std::endl;
                if (baseTimeListSize != interpolatedErrorMatrixSize){
                    std::runtime_error("observation lists unequal!");
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

                            ObservationVectorType newObservations = Eigen::Matrix< long double, Eigen::Dynamic, 1 >(allObservations.size());
                            ObservationVectorType observationWeights = Eigen::Matrix< long double, Eigen::Dynamic, 1 >(allObservations.size());

                            // for every observation, retrieve and add the range bias that should be added
                            for (unsigned int i=0; i<allObservationTimes.size(); i++){

                                double observationTime = allObservationTimes.at( i );

                                double rangeNoiseLevel = noiseLevelBasedOnMSEangleForMultipleMissions(
                                            observationTime, noiseAtMinAngleVector, noiseAtMaxAngleVector, maxMSEAngleDegVector, seperateBaseTimeLists);

                                double totalErrorLevel = combinedRangeAndSatelliteErrorLevel(
                                            observationTime, interpolatedErrorMatrix.at(observationTime), rangeNoiseLevel);

                                std::normal_distribution<double> d(0.0, totalErrorLevel);
                                long double rangeErrorSample = static_cast<long double>(d(gen));

                                if (podInputIterator->first == n_way_range){rangeErrorSample *= 2.0;}

                                if (rangeErrorSample > 10000.0){
                                    std::cout<<"ERROR: unusually large range correction at time "<<observationTime<<std::endl;
                                    std::cout<<" total error level: "<<totalErrorLevel<<std::endl;
                                    std::cout<<" range correction (from random sample): "<<rangeErrorSample<<std::endl;
                                }

                                newObservations(i) = allObservations(i) + rangeErrorSample;
                                observationWeights(i) = 1.0/(totalErrorLevel*totalErrorLevel);

        //                        std::cout<<observationTime<<" // "<<
        //                                   currentSatelliteError.transpose()<<" // "<<
        //                                   randomErrorSample.transpose()<<" // "<<
        //                                   rangeCorrection<<" // "<<
        //                                   noiseLevel<<std::endl;
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

            Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialParameterEstimate =
                    parametersToEstimate->template getFullParameterValues< long double >( );
            Eigen::Matrix< long double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

            // Perturb parameter estimate
            Eigen::Matrix< long double, Eigen::Dynamic, 1 > parameterPerturbation =
                    Eigen::Matrix< long double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ){
                // perturb body positions by  1 meters
                parameterPerturbation.segment(i*6,3) = Eigen::Vector3d::Constant( 1.0 ).cast<long double>();
                // perturb body velocities by 0.001 m/s
                parameterPerturbation.segment(i*6+3,3) = Eigen::Vector3d::Constant( 0.001 ).cast<long double>();
            }

            // perturb eta slightly to prevent partials from being 0 due to which the estimation is unable to estimate eta
            if (nordtvedtParameterIndex >= 0){
                parameterPerturbation(nordtvedtParameterIndex) = 1.0E-10;
            }

            initialParameterEstimate += parameterPerturbation;



            // the phase of J2/4 variations is usually not an estimatable parameter
            // however a different phase than the real one can be set here
            // this represents that a sine curve with a phase is estimated that is different from the reality
            relativity::variableJ2Interface->setPhase(J2phaseInEstimation);
            relativity::variableJ4Interface->setPhase(J4phaseInEstimation);


            std::cout << "True parameter values:" << std::endl;
            std::cout << truthParameters.transpose() << std::endl;
            std::cout << "Parameter perturbations:" << std::endl;
            std::cout << parameterPerturbation.transpose() << std::endl;
            std::cout << "Initial guesses:" << std::endl;
            std::cout << initialParameterEstimate.transpose() << std::endl;


            // Define a priori covariance matrix
            Eigen::MatrixXd aprioriCovariance =
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));

            for( unsigned int i = 0; i < truthParameters.size(); i++ ){
                aprioriCovariance( i,i ) = varianceVector.at( i );
            }
            std::cout << "a priori covariance matrix:" << std::endl
                      << aprioriCovariance.diagonal().transpose() << std::endl;
            Eigen::MatrixXd inverseOfAprioriCovariance = aprioriCovariance.inverse();

            // Define estimation input
            std::shared_ptr< PodInput< long double, double > > podInput =
                    std::make_shared< PodInput< long double, double > >(
                        observationsAndTimes, initialParameterEstimate.rows( ),
                        inverseOfAprioriCovariance,
                        initialParameterEstimate - truthParameters );
            podInput->defineEstimationSettings( true,
                                                reintegrateVariationalEquations,
                                                true, true, true, true,
                                                enforceNordtvedtConstraintInEstimation );
            podInput->manuallySetObservationWeightMatrix(observationWeightsAndTimes);

            // Perform estimation
            std::shared_ptr< PodOutput< long double > > podOutput = orbitDeterminationManager.estimateParameters(
                        podInput, std::make_shared< EstimationConvergenceChecker >( maximumNumberOfIterations ) );

            // Print true estimation error, limited mostly by numerical error
            Eigen::VectorXd estimationError = (podOutput->parameterEstimate_ - truthParameters).cast<double>();
            Eigen::VectorXd formalError = podOutput->getFormalErrorVector( ).cast<double>();

            std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
            std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;
            std::cout << "True to form estimation error ratio is: " << std::endl <<
                         ( estimationError.cwiseQuotient(  podOutput->getFormalErrorVector( ) ) ).transpose( ) << std::endl;

            std::cout << "improvement ratio formal w.r.t apriori error: " << std::endl <<
                         ( aprioriCovariance.diagonal().cwiseSqrt().cwiseQuotient( podOutput->getFormalErrorVector( ) ) ).transpose( ) << std::endl;


            std::cout<< "writing output to files that is not needed for the next steps..." << std::endl;

            // errors
            input_output::writeMatrixToFile( truthParameters, "TruthParameters.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( estimationError, "ObservationTrueEstimationError.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( formalError, "ObservationFormalEstimationError.dat", 16, outputSubFolder );

            // residuals, parameters
            input_output::writeMatrixToFile( podOutput->residuals_, "EstimationResiduals.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ), "ResidualHistory.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ), "ParameterHistory.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ), "EstimationCorrelations.dat", 16, outputSubFolder );

            // dependent variables
            std::map< double, Eigen::VectorXd > dependentVariablesHistoryFinalIteration;
            std::vector< std::map< double, Eigen::VectorXd > > dependentVariablesHistoryFinalIterationAllArcs =
                    podOutput->dependentVariableHistoryFinalIteration_;

            if (useMultipleMercuryArcs){
                for (unsigned int a=0; a<dependentVariablesHistoryFinalIterationAllArcs.size(); a++){
                    std::map< double, Eigen::VectorXd > dependentVariablesHistoryCurrentArc = dependentVariablesHistoryFinalIterationAllArcs.at(a);
                    if (dependentVariablesHistoryFinalIteration.size() == 0){
                        dependentVariablesHistoryFinalIteration = dependentVariablesHistoryCurrentArc;
                    } else{
                        std::map< double, Eigen::VectorXd >::iterator it = dependentVariablesHistoryCurrentArc.begin();
                        while (it != dependentVariablesHistoryCurrentArc.end()){
                            dependentVariablesHistoryFinalIteration.insert(std::make_pair(it->first, it->second));
                            it++;
                        }
                    }
                }
            } else{
                dependentVariablesHistoryFinalIteration = dependentVariablesHistoryFinalIterationAllArcs.at(0);
            }

            input_output::writeDataMapToTextFile(
                        onlyEveryXthValueFromDataMap(dependentVariablesHistoryFinalIteration, onlyEveryXthValue), "DependentVariablesHistoryFinalIteration.dat",
                        outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );


            /////////////////////////////
            //// CONSIDER PARAMETERS ////
            /////////////////////////////

            // get other required matrices for the calculation
            Eigen::MatrixXd initialCovarianceMatrix = podOutput->getUnnormalizedCovarianceMatrix( ); // P
            Eigen::VectorXd observationWeightDiagonal = podOutput->weightsMatrixDiagonal_; // diagonal of W
            Eigen::MatrixXd unnormalizedPartialDerivatives = podOutput->getUnnormalizedPartialDerivatives( ); // Hx
            Eigen::MatrixXd considerCovarianceMatrix;
            Eigen::MatrixXd considerCovarianceMatrixIncludingAsteroids;
            unsigned int maxcovtype = 1;

            // clear memory
            podOutput = nullptr;

            // start consider covariance analysis
            if (gammaIsAConsiderParameter || considerPPNalphas || considerSunAngularMomentum){

                std::cout<< "calculating covariance due to consider parameters..."<< std::endl;

                // In order to get the partial derivatives of consider parameters wrt observations, run an estimation of the consider parameters
                std::vector< std::shared_ptr < EstimatableParameterSettings > > considerParameterNames;
                std::vector<double> considerVarianceVector;

//                if (useMultipleMercuryArcs){
//                    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * initialTimeVector.size( ) );
//                    for( unsigned int m = 0; m < numberOfMissions; m++ ) {
//                        systemInitialState.segment( m * 6, 6 ) = propagatorSettingsList.at( m )->getInitialStates( );
//                        for( unsigned int i = 0; i < 3; i++ ){ considerVarianceVector.push_back(sigmaPosition*sigmaPosition); }
//                        for( unsigned int i = 3; i < 6; i++ ){ considerVarianceVector.push_back(sigmaVelocity*sigmaVelocity); }
//                    }
//                    considerParameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
//                                                  "Mercury", systemInitialState, initialTimeVector, "SSB" ) );
//                }
//                else{
                    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ) {
                        int j = 6*i;
                        considerParameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< long double > >(
                                                      bodiesToPropagate[i], systemInitialState.segment(j,6), centralBodies[i] ) );
                        for( unsigned int i = 0; i < 3; i++ ){ considerVarianceVector.push_back(sigmaPosition*sigmaPosition); }
                        for( unsigned int i = 3; i < 6; i++ ){ considerVarianceVector.push_back(sigmaVelocity*sigmaVelocity); }
                    }
//                }

                if ( gammaIsAConsiderParameter == true ){
                    considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                             ("global_metric", ppn_parameter_gamma ) );
                    considerVarianceVector.push_back(sigmaGamma*sigmaGamma);
                }

                std::cout<<calculateLenseThirringCorrection<<considerSunAngularMomentum<<std::endl;
                if ( calculateLenseThirringCorrection && considerSunAngularMomentum ){
                    considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                             ("Sun", angular_momentum));
                    considerVarianceVector.push_back(sigmaSunAngularMomentum*sigmaSunAngularMomentum);
                }

                considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_parameter_alpha1 ) );
                considerVarianceVector.push_back(sigmaAlpha1*sigmaAlpha1);

                considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                         ("global_metric", ppn_parameter_alpha2 ) );
                considerVarianceVector.push_back(sigmaAlpha2*sigmaAlpha2);

                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > considerParametersToEstimate =
                        createParametersToEstimate( considerParameterNames, bodyMap, propagatorSettings );
                printEstimatableParameterEntries( considerParametersToEstimate );

                OrbitDeterminationManager< long double, double > orbitDeterminationManagerConsiderParameters =
                        OrbitDeterminationManager< long double, double >(
                            bodyMap, considerParametersToEstimate, observationSettingsMap,
                            integratorSettings, propagatorSettings );

                Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialConsiderParameterEstimate =
                        considerParametersToEstimate->template getFullParameterValues< long double >( );
                Eigen::Matrix< long double, Eigen::Dynamic, 1 > truthConsiderParameters = initialConsiderParameterEstimate;

                Eigen::MatrixXd considerParameterAprioriCovariance = // C
                    Eigen::MatrixXd::Zero( truthConsiderParameters.rows( ), truthConsiderParameters.rows( ));

                for( unsigned int i = 0; i < truthConsiderParameters.size(); i++ ){
                    considerParameterAprioriCovariance( i,i ) = considerVarianceVector.at( i );
                }
                std::cout << "consider parameter a priori covariance matrix:" << std::endl
                          << considerParameterAprioriCovariance.diagonal().transpose() << std::endl;

                std::cout<<"input consider parameter estimation:"<<std::endl;
                std::cout<<initialConsiderParameterEstimate.transpose()<<std::endl;
                std::cout<<truthConsiderParameters.transpose()<<std::endl;
                std::cout<<(initialConsiderParameterEstimate-truthConsiderParameters).transpose()<<std::endl;

                std::shared_ptr< PodInput< long double, double > > podInputConsiderParameters =
                        std::make_shared< PodInput< long double, double > >(
                            observationsAndTimes, initialConsiderParameterEstimate.rows( ),
                            considerParameterAprioriCovariance.inverse(),
                            initialConsiderParameterEstimate - truthConsiderParameters );
                podInputConsiderParameters->defineEstimationSettings( true, false, true, true );
                podInputConsiderParameters->manuallySetObservationWeightMatrix(observationWeightsAndTimes);

                std::shared_ptr< PodOutput< long double > > podOutputConsiderParameters = orbitDeterminationManagerConsiderParameters.estimateParameters(
                            podInputConsiderParameters, std::make_shared< EstimationConvergenceChecker >( 1 ) );

                Eigen::MatrixXd partialDerivativesOfConsiderParameters = // H_c
                        (podOutputConsiderParameters->getUnnormalizedPartialDerivatives()).rightCols(
                            truthConsiderParameters.size()-6*numberOfNumericalBodies);

                // calculate covariance matrix including contribution from consider parameters
                considerCovarianceMatrix = calculateConsiderCovarianceMatrix(
                            initialCovarianceMatrix, observationWeightDiagonal,
                            considerParameterAprioriCovariance.block(6,6,considerParameterAprioriCovariance.rows()-6, considerParameterAprioriCovariance.cols()-6),
                            unnormalizedPartialDerivatives, partialDerivativesOfConsiderParameters,
                            outputSubFolder);

                // get consider correlation matrix
                Eigen::VectorXd formalErrorWithConsiderParameters = considerCovarianceMatrix.diagonal( ).cwiseSqrt( );
                Eigen::MatrixXd considerCorrelationMatrix = considerCovarianceMatrix.cwiseQuotient(
                            formalErrorWithConsiderParameters * formalErrorWithConsiderParameters.transpose() );

                input_output::writeMatrixToFile( considerParameterAprioriCovariance, "ConsiderParameterAprioriCovariance.dat", 16, outputSubFolder );
                input_output::writeMatrixToFile( considerCorrelationMatrix, "EstimationConsiderCorrelations.dat", 16, outputSubFolder );
                input_output::writeMatrixToFile( formalErrorWithConsiderParameters, "ObservationFormalEstimationErrorWithConsiderParameters.dat", 16, outputSubFolder );

                maxcovtype = 2;

                // include asteroids
                if (testCeres == false){
                    considerCovarianceMatrixIncludingAsteroids =
                            considerCovarianceMatrix + calculateConsiderCovarianceOfAsteroids(
                                initialCovarianceMatrix, observationWeightDiagonal,
                                unnormalizedPartialDerivatives,
                                json_directory,
                                json_directory + "/asteroids_multiplemissions",
                                outputSubFolder);

                    Eigen::VectorXd formalErrorWithConsiderParametersIncludingAsteroids = considerCovarianceMatrixIncludingAsteroids.diagonal( ).cwiseSqrt( );
                    Eigen::MatrixXd considerCorrelationMatrixIncludingAsteroids = considerCovarianceMatrixIncludingAsteroids.cwiseQuotient(
                                formalErrorWithConsiderParametersIncludingAsteroids * formalErrorWithConsiderParametersIncludingAsteroids.transpose() );

                    input_output::writeMatrixToFile( considerCorrelationMatrixIncludingAsteroids, "EstimationConsiderIncludingAsteroidsCorrelations.dat", 16, outputSubFolder );
                    input_output::writeMatrixToFile( formalErrorWithConsiderParametersIncludingAsteroids, "ObservationFormalEstimationErrorWithConsiderIncludingAsteroidsParameters.dat", 16, outputSubFolder );

                    maxcovtype = 3;
                }
            }

            // clear memory
            bodyMap.clear();



            /////////////////////////////////////////////
            //// PROVIDE OUTPUT TO CONSOLE AND FILES ////
            /////////////////////////////////////////////


            // propagate according to integration history. earlier result is separated here as the times are needed on their own.
            std::vector<double> fullStateHistoryTimes;
            std::map<double, Eigen::Matrix<long double, Eigen::Dynamic, 1>> propagationHistory = allBodiesPropagationHistory.at( 0 );
            std::map<double, Eigen::Matrix<long double, Eigen::Dynamic, 1>>::iterator historyit = propagationHistory.begin();

            while (historyit != propagationHistory.end()){
                fullStateHistoryTimes.push_back(historyit->first);
                historyit++;
            }

            // Propagate covariance matrix (twice, both with and without the consider parameters included in the covariance)
            std::map<double, Eigen::MatrixXd > propagatedCovariance;
            std::map<double, Eigen::Vector6d> propagatedErrorUsingCovMatrix;
            std::map<double, Eigen::Vector6d> propagatedRSWErrorUsingCovMatrix;

            for (unsigned int covtype = 0; covtype<maxcovtype; covtype++){

                Eigen::MatrixXd initialCovariance;
                std::string saveString;
                if (covtype == 0){
                    std::cout<<"propagating covariance matrix..."<<std::endl;
                    initialCovariance = initialCovarianceMatrix;
                    saveString = "";
                } else if (covtype == 1){
                    std::cout<<"propagating consider covariance matrix..."<<std::endl;
                    initialCovariance = considerCovarianceMatrix;
                    saveString = "Consider";
                } else{
                    std::cout<<"propagating consider covariance matrix including asteroids..."<<std::endl;
                    initialCovariance = considerCovarianceMatrixIncludingAsteroids;
                    saveString = "ConsiderIncludingAsteroids";
                }

                propagateCovariance(propagatedCovariance, initialCovariance,
                                    orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                                    fullStateHistoryTimes);
                                    //baseTimeList
                                    //60.0, initialEphemerisTime + 3600.0, finalEphemerisTime - 3600.0 );

                std::map<double, double> trueAnomalyMap;


                std::map<double, Eigen::MatrixXd>::iterator it = propagatedCovariance.begin();

                while (it != propagatedCovariance.end()){

                    Eigen::Vector6d currentCartesianState = propagationHistory.at(it->first).cast<double>();

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


                input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedCovariance, onlyEveryXthValue),
                                                      "Propagated"+saveString+"Covariance.dat", outputSubFolder );
                input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedErrorUsingCovMatrix, onlyEveryXthValue),
                                                                                   "propagatedErrorUsing"+saveString+"CovMatrix.dat", outputSubFolder );
                input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedRSWErrorUsingCovMatrix, onlyEveryXthValue),
                                                      "propagatedRSWErrorUsing"+saveString+"CovMatrix.dat", outputSubFolder );

                // clear memory
                propagatedCovariance.clear();
                propagatedErrorUsingCovMatrix.clear();
                propagatedRSWErrorUsingCovMatrix.clear();

            }

            // Save data in files
            std::cout<<"writing output to files..."<<std::endl;

            // covariance and correlation
            input_output::writeMatrixToFile( considerCovarianceMatrix, "ConsiderCovarianceMatrix.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( considerCovarianceMatrixIncludingAsteroids, "ConsiderIncludingAsteroidsCovarianceMatrix.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( initialCovarianceMatrix, "InitialCovarianceMatrix.dat", 16, outputSubFolder );
            input_output::writeMatrixToFile( observationWeightDiagonal, "ObservationWeightDiagonal.dat", 16, outputSubFolder );

        //        input_output::writeMatrixToFile( unnormalizedPartialDerivatives, "test_unnormalizedPartialDerivatives.dat", 16, outputSubFolder );
        //        input_output::writeMatrixToFile( partialDerivativesOfConsiderParameters, "test_partialDerivativesOfConsiderParameters.dat", 16, outputSubFolder );

        //        input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_, "EstimationInformationMatrix.dat", 16, outputSubFolder );
        //        input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_, "EstimationInformationMatrixNormalization.dat", 16, outputSubFolder );
        //        input_output::writeMatrixToFile( podInput->getInverseOfAprioriCovariance( ), "InverseAPrioriCovariance.dat", 16, outputSubFolder );
        //        input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_, "InverseNormalizedCovariance.dat", 16, outputSubFolder );
        //        input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ), "UnnormalizedCovariance.dat", 16, outputSubFolder );

            // observations
            input_output::writeDataMapToTextFile( interpolatedErrorMatrix, "interpolatedErrorMatrix.dat", outputSubFolder );
//            input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ), "ObservationMeasurements.dat", 16, outputSubFolder );
//            input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
//                                             "ObservationTimes.dat", 16, outputSubFolder );
//            input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
//                                             "ObservationLinkEnds.dat", 16, outputSubFolder );
//            input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
//                                             "ObservationObservableTypes.dat", 16, outputSubFolder );
//            input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
//                                             "ObservationMeasurements.dat", 16, outputSubFolder );



            std::cout << "done!" << std::endl;       

        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;

}
