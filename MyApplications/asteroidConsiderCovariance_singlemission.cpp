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
#include "Tudat/Astrodynamics/Relativity/variableJ2Interface.h"
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
    const bool calculateDeSitterCorrection = false;
    const bool includeTimeVaryingGravitationalMomentsSun = false;
    unsigned int maximumDegreeSunGravitationalMoment = 2;
    unsigned int maximumDegreeSunGravitationalMomentVariation = 2;

    // Parameter estimation settings
    const unsigned int maximumNumberOfIterations = 1;
    const double sigmaPosition = 1000.0; //educated guess
    const double sigmaVelocity = 1.0; //educated guess
    const bool ignoreNordtvedtConstraintInEstimation = false;
    const bool includeSpacecraftPositionError = true;
    const bool includeLightTimeCorrections = false;


    // ABM integrator settings (if RK4 is used instead, initialstepsize is taken)
    const double initialTimeStep = 3600.0;
    const double minimumStepSize = 3600.0;
    const double maximumStepSize = 3600.0;
    const double relativeErrorTolerence = 1.0;
    const double absoluteErrorTolerence = 1.0;
    const unsigned int minimumOrder = 8;
    const unsigned int maximumOrder = 8;

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

    // Load json input
    std::vector< std::string > filenames;
    filenames.push_back("inputs_Genova2018.json"); // 0
    filenames.push_back("inputs_Schettino2015.json"); // 1
    filenames.push_back("inputs_Schettino2015_alphas.json"); // 2
    filenames.push_back("inputs_Imperi2018_nvtrue_flybys.json"); // 3
    filenames.push_back("inputs_Imperi2018_nvtrue_flybys_alphas.json"); // 4
    filenames.push_back("inputs_Imperi2018_nvfalse_flybys.json"); // 5
    filenames.push_back("inputs_Imperi2018_nvfalse_flybys_alphas.json"); // 6

    for (unsigned int f = 1; f<2; f++){
//    for (unsigned int f = 0; f<filenames.size(); f++){

        std::string input_filename = filenames.at(f);
        std::cout<<"---- RUNNING SIMULATION FOR INPUTS WITH FILENAME: "<<input_filename<<" ----"<<std::endl;

        std::string json_directory = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/";
        std::ifstream json_file(json_directory + input_filename);
        json json_input;
        json_file >> json_input;

        std::cout<<"input values imported: "<<json_input<<std::endl;


        // Retrieve simulation start and end time
        std::vector<int> json1 = json_input["initialTime"];
        Eigen::Vector6i initialTime(json1.data()); json1.clear();
        std::vector<int> json2 = json_input["finalTime"];
        Eigen::Vector6i finalTime(json2.data()); json2.clear();

        std::string vehicle = json_input["vehicle"];

        // Acceleration settings
        const bool calculateSchwarzschildCorrection = json_input["calculateSchwarzschildCorrection"];
        const bool calculateLenseThirringCorrection = json_input["calculateLenseThirringCorrection"];
        const bool includeSEPViolationAcceleration = json_input["includeSEPViolationAcceleration"];
        const bool includeTVGPAcceleration = json_input["includeTVGPAcceleration"];
        const bool estimateJ2Amplitude = json_input["estimateJ2Amplitude"];
        const bool estimateJ2Period = json_input["estimateJ2Period"];
        const bool estimateJ2Phase = json_input["estimateJ2Phase"];
        const bool estimateJ4Amplitude = json_input["estimateJ4Amplitude"];
        const bool estimateJ4Period = json_input["estimateJ4Period"];
        const bool estimateJ4Phase = json_input["estimateJ4Phase"];


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
        const bool nordtvedtConstraintTrueOrFalse = json_input["useNordtvedtConstraint"];
        const bool estimateSunAngularMomentum = json_input["estimateSunAngularMomentum"];
        const bool estimatePPNalphas = json_input["estimatePPNalphas"];
        bool ppnAlphasAreConsiderParameters = json_input["ppnAlphasAreConsiderParameters"];
        const bool gammaIsAConsiderParameter = json_input["gammaIsAConsiderParameter"];

        // retreive regular observation schedule
        std::vector<int> json3 = json_input["observationInitialTime"];
        Eigen::Vector6i observationInitialTime(json3.data()); json3.clear();
        std::vector<int> json4 = json_input["observationFinalTime"];
        Eigen::Vector6i observationFinalTime(json4.data()); json4.clear();
        const double observationTimeStep = json_input["observationTimeStep"];
        const double trackingArcDuration = json_input["trackingArcDuration"];
        const unsigned int maximumNumberOfTrackingDays = json_input["maximumNumberOfTrackingDays"];

        // retreive flybys
        std::vector<int> flybyObject = json_input["flybyObject"];
        std::vector<double> flybyList;
        std::vector<int> currentFlyby;
        for (unsigned int i=0; i<flybyObject.size()/6; i++ ){

            for (unsigned int c=0; c<6; c++ ){
                currentFlyby.push_back(flybyObject.at(6*i+c));
            }
            flybyList.push_back(secondsAfterJ2000(currentFlyby));
            currentFlyby.clear();
        }

        // observation Noise
        const double noiseAtMinAngle = json_input["noiseAtMinAngle"];
        const double noiseAtMaxAngle = json_input["noiseAtMaxAngle"];
        const double maxMSEAngleDeg = json_input["maxMSEAngleDeg"];

        // Location of simulation input and output
        std::string outputSubFolderName = json_input["outputSubFolderName"];
        std::string inputSubFolder = getOutputPath( ) + outputSubFolderName;
        std::string outputSubFolder = inputSubFolder + "_asteroids";

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


        // some checks to prevent incompatible inputs
        if ((includeTimeVaryingGravitationalMomentsSun == false) &&
                (estimateJ2Amplitude || estimateJ2Period || estimateJ2Phase || estimateJ4Amplitude || estimateJ4Period || estimateJ4Phase)){
            std::runtime_error("cannot estimate time varying gravitational parameters when includeTimeVaryingGravitationalMoments is set to false");
        }
        if ((calculateLenseThirringCorrection == false) && (estimateSunAngularMomentum)){
            std::runtime_error("cannot estimate sun angular momentum when lense thirring acceleration is not calculated");
        }
        if ((estimatePPNalphas) && (ppnAlphasAreConsiderParameters)){
            std::runtime_error("ppn alphas cannot both be estimatable parameters and consider parameters");
        }
        if ((estimateJ4Amplitude || estimateJ4Period || estimateJ4Phase)
                && maximumDegreeSunGravitationalMomentVariation < 4){
            std::runtime_error("incompatible settings, cannot estimate sine parameters of J4 when J4 is not included.");
        }




        /////////////////////
        //// ENVIRONMENT ////
        /////////////////////

        const unsigned int increment = 10;
        for (unsigned int inc=0; inc<10; inc += increment){

            unsigned int minIndex = inc;
            unsigned int maxIndex = inc+increment-1;
            std::cout << ">>>> RUNNING ASTEROID ANALYSIS FOR mu.txt INDICES "<<minIndex<<"-"<<maxIndex<<std::endl;

            std::cout << "building environment..." << std::endl;

            // Convert initial and final time to Julian
            double initialSimulationTime = secondsAfterJ2000(initialTime);
            double finalSimulationTime = secondsAfterJ2000(finalTime);

            std::cout<<"initial time: "<<initialSimulationTime<<" final time: "<<finalSimulationTime<<std::endl;

            // Load spice kernels.
            std::string kernelsPath = input_output::getSpiceKernelPath( );

            std::vector< std::string > customKernels;
            customKernels.push_back( kernelsPath + "tudat_merged_spk_kernel_thesis4.bsp" );
            spice_interface::loadStandardSpiceKernels( customKernels );

            spice_interface::loadSpiceKernelInTudat( kernelsPath + "codes_300ast_20100725.tf" );
            spice_interface::loadSpiceKernelInTudat( kernelsPath + "codes_300ast_20100725.bsp" );

    //        Eigen::Vector6i testEpochVector; testEpochVector << 2020, 1, 1, 0, 0, 0;
    //        double testEpoch = secondsAfterJ2000(testEpochVector);
    //        std::cout<<"testing asteroid kernels at 2020/1/1 00:00:00"<<std::endl;
    //        std::cout<<"    asteroid 1 state: "<<spice_interface::getBodyCartesianStateAtEpoch("2000001","SSB","ECLIPJ2000","None",testEpoch).transpose()<<std::endl;
    //        std::cout<<"    asteroid 2 state: "<<spice_interface::getBodyCartesianStateAtEpoch("2000002","SSB","ECLIPJ2000","None",testEpoch).transpose()<<std::endl;
    //        std::cout<<"    asteroid 3 state: "<<spice_interface::getBodyCartesianStateAtEpoch("2000003","SSB","ECLIPJ2000","None",testEpoch).transpose()<<std::endl;
    //        std::cout<<"    asteroid 4 state: "<<spice_interface::getBodyCartesianStateAtEpoch("2000004","SSB","ECLIPJ2000","None",testEpoch).transpose()<<std::endl;
    //        std::cout<<"    asteroid 5 state: "<<spice_interface::getBodyCartesianStateAtEpoch("2000005","SSB","ECLIPJ2000","None",testEpoch).transpose()<<std::endl;

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

            // load stuff of asteroids
            Eigen::VectorXd asteroidNumbersDouble = input_output::readMatrixFromFile(json_directory + "mu.txt", " ").col(0).segment(minIndex, increment);
            Eigen::VectorXi asteroidNumbers = asteroidNumbersDouble.cast<int>();
            std::cout<<"asteroid list from JPL: "<<asteroidNumbers.transpose()<<std::endl;

            const unsigned int maxAsteroid = asteroidNumbers.size();
            std::map< unsigned int, std::pair< double, double > > asteroidsINPOP19a
                    = readAsteroidsFile(json_directory + "AsteroidsINPOP19a.txt", " ");

            // loop over asteroids
            for (unsigned int a=0; a<maxAsteroid; a++){
                // check if key exists in map (i.e. asteroid is included in INPOP file)
                if (asteroidsINPOP19a.find(asteroidNumbers(a)) != asteroidsINPOP19a.end()){
                    // add SPICE number of asteroid to bodies
                    bodyNames.push_back(std::to_string(2000000 + asteroidNumbers(a)));
                } else{
                    std::cout<<"asteroid "<<asteroidNumbers(a)<<" not found in INPOP19a list"<<std::endl;
                }

            }
            const unsigned int totalNumberOfBodies = bodyNames.size();


            // Default body settings
            double buffer = maximumOrder*maximumStepSize; //see Tudat libraries 1.1.3.
            std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
            bodySettings = getDefaultBodySettings( bodyNames, initialSimulationTime - buffer, finalSimulationTime + buffer );


            // Simplification, for outer solar sytem bodies, take tabulated values
            bodySettings["Jupiter"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Saturn"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Uranus"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );
            bodySettings["Neptune"]->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );


            // asteroid settings
            std::cout << "setting custom settings Asteroids..." << std::endl;

            for (unsigned int a=0; a<maxAsteroid; a++){
                // check if key exists in map (i.e. asteroid is included in file)
                if (asteroidsINPOP19a.find(asteroidNumbers(a)) != asteroidsINPOP19a.end()){

                    // set ephemeris settings to tabulated with 1 hour time step
                    bodySettings[ std::to_string(2000000 + asteroidNumbers(a)) ]->ephemerisSettings
                            = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                                initialSimulationTime - buffer, finalSimulationTime + buffer, 3600.0, "SSB", "ECLIPJ2000" );

                    // set gravitational parameter as given in INPOP19a
                    bodySettings[ std::to_string(2000000 + asteroidNumbers(a)) ]->gravityFieldSettings
                            = std::make_shared< CentralGravityFieldSettings >(asteroidsINPOP19a.at(asteroidNumbers(a)).first * convertAsteroidGMtoSI);

                    std::cout<<"asteroid "<<asteroidNumbers(a)
                             <<" GM: "<<asteroidsINPOP19a.at(asteroidNumbers(a)).first<<" *10^-18 AU3/d2, converted: "
                             <<asteroidsINPOP19a.at(asteroidNumbers(a)).first * convertAsteroidGMtoSI<<" m3/s2"
                             <<"     "
                             <<" sigma: "<<asteroidsINPOP19a.at(asteroidNumbers(a)).second<<" *10^-18 AU3/d2, converted: "
                             <<asteroidsINPOP19a.at(asteroidNumbers(a)).second * convertAsteroidGMtoSI<<" m3/s2"<<std::endl;
                }
            }



            std::cout << "setting custom settings Sun..." << std::endl;

            // Nominal values spherical harmonics Sun
            Eigen::MatrixXd normalizedSineCoefficients
                    = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMoment + 1,
                                            maximumDegreeSunGravitationalMoment+ 1);
            Eigen::MatrixXd normalizedCosineCoefficients
                    = Eigen::MatrixXd::Zero(maximumDegreeSunGravitationalMoment + 1,
                                            maximumDegreeSunGravitationalMoment + 1);

            if (normalizedCosineCoefficients.col(0).size() != valuesSunGravitationalMoments.size()){
                std::runtime_error("error: vector sizes incompatible");
            }

            normalizedCosineCoefficients.col(0) = valuesSunGravitationalMoments.transpose();

            bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        sunGravitationalParameter, sunRadius,
                        normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );


    //        bodySettings[ "Sun" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
    //                    "ECLIPJ2000", "IAU_Sun",
    //                    spice_interface::computeRotationQuaternionBetweenFrames("ECLIPJ2000", "IAU_Sun", initialSimulationTime ),
    //                    initialSimulationTime, 2.0 * mathematical_constants::PI / solarDay,
    //                    sunAngularMomentum );


    //        // Time varying spherical harmonics coefficients Sun
    //        if (includeTimeVaryingGravitationalMomentsSun){

    //            // J2
    //            std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariationSettings;

    //            relativity::variableJ2Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(2));
    //            relativity::variableJ2Interface->setPeriod(periodsSunGravitationalMomentsVariation(2));
    //            relativity::variableJ2Interface->setPhase(phasesSunGravitationalMomentsVariation(2));

    //            std::function< double() > amplitudeFunctionJ2 =
    //                    std::bind( &relativity::VariableJ2Interface::getAmplitude, relativity::variableJ2Interface );
    //            std::function< double() > periodFunctionJ2 =
    //                    std::bind( &relativity::VariableJ2Interface::getPeriod, relativity::variableJ2Interface );
    //            std::function< double() > phaseFunctionJ2 =
    //                    std::bind( &relativity::VariableJ2Interface::getPhase, relativity::variableJ2Interface );


    //            std::function< std::map< double, Eigen::MatrixXd >( ) > tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ2 =
    //                    std::bind(tabulatedSphericalHarmonicsCoefficientCorrections,
    //                              initialSimulationTime, finalSimulationTime,
    //                              amplitudeFunctionJ2, periodFunctionJ2, phaseFunctionJ2);

    //            std::map< double, Eigen::MatrixXd > sineCoefficientCorrectionsJ2 =
    //                    zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

    //            const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettingsJ2 =
    //                    std::make_shared< TabulatedGravityFieldVariationSettingsWithCosineFunction >(
    //                        tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ2,
    //                        sineCoefficientCorrectionsJ2,
    //                        2, 0, interpolatorSettings,
    //                        initialSimulationTime, finalSimulationTime, 3600.0);

    //            gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettingsJ2 );

    //            //J4
    //            if (maximumDegreeSunGravitationalMomentVariation >= 4){

    //                relativity::variableJ4Interface->setAmplitude(amplitudesSunGravitationalMomentsVariation(4));
    //                relativity::variableJ4Interface->setPeriod(periodsSunGravitationalMomentsVariation(4));
    //                relativity::variableJ4Interface->setPhase(phasesSunGravitationalMomentsVariation(4));

    //                std::function< double() > amplitudeFunctionJ4 =
    //                        std::bind( &relativity::VariableJ4Interface::getAmplitude, relativity::variableJ4Interface );
    //                std::function< double() > periodFunctionJ4 =
    //                        std::bind( &relativity::VariableJ4Interface::getPeriod, relativity::variableJ4Interface );
    //                std::function< double() > phaseFunctionJ4 =
    //                        std::bind( &relativity::VariableJ4Interface::getPhase, relativity::variableJ4Interface );


    //                std::function< std::map< double, Eigen::MatrixXd >( ) > tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ4 =
    //                        std::bind(tabulatedSphericalHarmonicsCoefficientCorrections,
    //                                  initialSimulationTime, finalSimulationTime,
    //                                  amplitudeFunctionJ4, periodFunctionJ4, phaseFunctionJ4);

    //                std::map< double, Eigen::MatrixXd > sineCoefficientCorrectionsJ4 =
    //                        zeroTabulatedSphericalHarmonicsCoefficientCorrections(initialSimulationTime, finalSimulationTime);

    //                const std::shared_ptr< GravityFieldVariationSettings > timeVaryingSphericalHarmonicsSettingsJ4 =
    //                        std::make_shared< TabulatedGravityFieldVariationSettingsWithCosineFunction >(
    //                            tabulatedSphericalHarmonicsCoefficientCorrectionsFunctionJ4,
    //                            sineCoefficientCorrectionsJ4,
    //                            2, 0, interpolatorSettings,
    //                            initialSimulationTime, finalSimulationTime, 3600.0);

    //                gravityFieldVariationSettings.push_back( timeVaryingSphericalHarmonicsSettingsJ4 );
    //            }


    //            bodySettings[ "Sun" ]->gravityFieldVariationSettings = gravityFieldVariationSettings;
    //        }
    //        // Prepare angular momentum vector Sun
    ////        const Eigen::Vector3d sunAngularMomentumVectorInSunFrame(0.0, 0.0, sunAngularMomentum);
    ////        const Eigen::Vector3d sunAngularMomentumVectorPerUnitMassInSunFrame =
    ////                sunAngularMomentumVectorInSunFrame /
    ////                (sunGravitationalParameter/physical_constants::GRAVITATIONAL_CONSTANT);

            std::cout << "creating environment..." << std::endl;

            // Create body map
            NamedBodyMap bodyMap = createBodies( bodySettings );
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
                                        std::make_shared< SphericalHarmonicAccelerationSettings > (maximumDegreeSunGravitationalMoment,0));

    //                        if (calculateSchwarzschildCorrection || calculateLenseThirringCorrection || calculateDeSitterCorrection){
                                currentAccelerations[ bodyNames.at( j ) ].push_back(
                                            std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                true, false, false));
    //                        }

    //                        if (includeSEPViolationAcceleration == true){
    //                            currentAccelerations[ bodyNames.at( j ) ].push_back(
    //                                        std::make_shared< SEPViolationAccelerationSettings >(
    //                                            bodyNames, nordtvedtConstraintTrueOrFalse, ignoreNordtvedtConstraintInEstimation));
    //                        }

    //                        if (includeTVGPAcceleration == true){
    //                        currentAccelerations[ bodyNames.at( j ) ].push_back(
    //                                    std::make_shared< TimeVaryingGravitationalParameterAccelerationSettings >(
    //                                        timeVaryingGravitationalParameter));
    //                        }

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
                if ( !(bodyNames.at( i ) == "Mercury") && !(bodyNames.at( i ) == "Sun")){
                    dependentVariablesList.push_back(
                              std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                    central_gravity, "Mercury" , bodyNames.at( i ) ) );
                }
            }

            dependentVariablesList.push_back(
                        std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                        "Mercury", "Sun", maximumDegreeSunGravitationalMoment, 0 ) );

            if (calculateSchwarzschildCorrection){
                dependentVariablesList.push_back(
                          std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 1 ) );
                dependentVariablesList.push_back(
                          std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 2 ) );
            }

    //        if (calculateLenseThirringCorrection){
    //            dependentVariablesList.push_back(
    //                      std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
    //                            relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 3 ) );
    //        }

    //        if (calculateDeSitterCorrection){
    //            dependentVariablesList.push_back(
    //                      std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
    //                            relativistic_correction_acceleration, "Mercury" , "Sun", 0, -1, 4 ) );
    //        }

    //        if (includeSEPViolationAcceleration){
    //        dependentVariablesList.push_back(
    //                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
    //                        sep_violation_acceleration, "Mercury" , "Sun" ) );
    //        }

    //        if (includeTVGPAcceleration){
    //        dependentVariablesList.push_back(
    //                  std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
    //                        time_varying_gravitational_parameter_acceleration, "Mercury" , "Sun" ) );
    //        }

            // Create object with list of dependent variables
            std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );


            //////////////////////////////
            //// PROPAGATION SETTINGS ////
            //////////////////////////////

            std::cout << "defining propagation settings..." << std::endl;

            // Get initial state of bodies to be propagated
            Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                        bodiesToPropagate, centralBodies, bodyMap, initialSimulationTime );

            // Define propagator settings.
            std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                  std::make_shared< propagators::PropagationTimeTerminationSettings >( finalSimulationTime, true );

            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate,
                      systemInitialState, terminationSettings, cowell, dependentVariablesToSave);


            std::shared_ptr< AdamsBashforthMoultonSettings< double > > integratorSettings =
                    std::make_shared< AdamsBashforthMoultonSettings< double > > (
                        initialSimulationTime, initialTimeStep,
                        minimumStepSize, maximumStepSize,
                        relativeErrorTolerence, absoluteErrorTolerence,
                        minimumOrder, maximumOrder);

//            std::shared_ptr< IntegratorSettings< > > integratorSettings =
//                    std::make_shared< IntegratorSettings< > >
//                    ( rungeKutta4, initialSimulationTime, initialTimeStep );



            ////////////////////////////
            //// DYNAMICS SIMULATOR ////
            ////////////////////////////

            std::cout << "running dynamics simulator..." << std::endl;

            SingleArcDynamicsSimulator <> dynamicsSimulator (bodyMap,integratorSettings,propagatorSettings,true,false,true);

            std::cout << "saving integration result and dependent variables..." << std::endl;

            std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::VectorXd > dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory();

            // Retrieve numerically integrated state for each body.
            std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
            std::vector< std::map< double, Eigen::VectorXd > > spiceStatesAtPropagationTimes;
            allBodiesPropagationHistory.resize( bodiesToPropagate.size() );
//            spiceStatesAtPropagationTimes.resize( bodiesToPropagate.size() );

            for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
                 stateIterator != integrationResult.end( ); stateIterator++ )
            {
                for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
                {
                    allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
//                    spiceStatesAtPropagationTimes[ i ][ stateIterator->first ] =
//                            getBodyCartesianStateAtEpoch(bodiesToPropagate.at( i ),"SSB","ECLIPJ2000","None",stateIterator->first);
                }
            }

            for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
            {
                // Write propagation history to file.
                input_output::writeDataMapToTextFile(
                            onlyEveryXthValueFromDataMap( allBodiesPropagationHistory[ i ], onlyEveryXthValue),
                            "StatePropagationHistory" + bodiesToPropagate.at( i ) + ".dat",
                            outputSubFolder,
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    //            input_output::writeDataMapToTextFile(
    //                        spiceStatesAtPropagationTimes[ i ],
    //                        "spiceStatesAtPropagationTimes" + bodiesToPropagate.at( i ) + ".dat",
    //                        outputSubFolder,
    //                        "",
    //                        std::numeric_limits< double >::digits10,
    //                        std::numeric_limits< double >::digits10,
    //                        "," );
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



    //        //////////////////////////////
    //        //// BACKWARD PROPAGATION ////
    //        //////////////////////////////

    //        std::cout << "performing backward propagation..." << std::endl;

    //        std::shared_ptr< PropagationTimeTerminationSettings > backwardTerminationSettings =
    //              std::make_shared< propagators::PropagationTimeTerminationSettings >( initialSimulationTime, true );

    //        Eigen::VectorXd systemFinalState = allBodiesPropagationHistory[ 0 ].at(finalSimulationTime);

    //        std::shared_ptr< TranslationalStatePropagatorSettings< double > > backwardPropagatorSettings =
    //                std::make_shared< TranslationalStatePropagatorSettings< double > >
    //                ( centralBodies, accelerationModelMap, bodiesToPropagate,
    //                  systemFinalState, backwardTerminationSettings);

    //        std::shared_ptr< AdamsBashforthMoultonSettings< double > > backwardIntegratorSettings =
    //                std::make_shared< AdamsBashforthMoultonSettings< double > > (
    //                    finalSimulationTime, -1.0*initialTimeStep,
    //                    -1.0*minimumStepSize, -1.0*maximumStepSize,
    //                    relativeErrorTolerence, absoluteErrorTolerence,
    //                    minimumOrder, maximumOrder);

    //        std::cout << "running dynamics simulator..." << std::endl;

    //        SingleArcDynamicsSimulator <> backwardDynamicsSimulator (bodyMap,backwardIntegratorSettings,backwardPropagatorSettings,true,false,true);

    //        std::cout << "saving integration result and dependent variables..." << std::endl;

    //        std::map< double, Eigen::VectorXd > backwardIntegrationResult = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    //        // Retrieve numerically integrated state for each body.
    //        std::vector< std::map< double, Eigen::VectorXd > > backwardAllBodiesPropagationHistory;
    //        backwardAllBodiesPropagationHistory.resize( bodiesToPropagate.size() );
    //        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = backwardIntegrationResult.begin( );
    //             stateIterator != backwardIntegrationResult.end( ); stateIterator++ )
    //        {
    //            for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
    //            {
    //                backwardAllBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
    //            }
    //        }

    //        for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
    //        {
    //            // Write propagation history to file.
    //            input_output::writeDataMapToTextFile(
    //                        backwardAllBodiesPropagationHistory[ i ],
    //                        "StatePropagationHistory" + bodiesToPropagate.at( i ) + "Backwards.dat",
    //                        outputSubFolder,
    //                        "",
    //                        std::numeric_limits< double >::digits10,
    //                        std::numeric_limits< double >::digits10,
    //                        "," );
    //        }


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

            std::vector< std::shared_ptr < EstimatableParameterSettings > > parameterNames;
            std::vector<double> varianceVector;

            // Add bodies that will be propagated to the parameters to be estimated
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

            // set gravitational parameters asteroids
            for (unsigned int a=0; a<maxAsteroid; a++){
                // check if key exists in map (i.e. asteroid is included in file)
                if (asteroidsINPOP19a.find(asteroidNumbers(a)) != asteroidsINPOP19a.end()){

                    // set GM of the asteroid as estimatable parameter, get apriori sigma from INPOP19a
                    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                                             (std::to_string(2000000 + asteroidNumbers(a)), gravitational_parameter));
                    double asteroidSigma = asteroidsINPOP19a.at(asteroidNumbers(a)).second * convertAsteroidGMtoSI;
                    varianceVector.push_back(asteroidSigma*asteroidSigma);

                }
            }

    //        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                 ("2000001", gravitational_parameter));
    //        double asteroidSigma = 340.331E18 * convertAUperdayToSI;
    //        varianceVector.push_back(asteroidSigma*asteroidSigma);

    //        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                 ("2000002", gravitational_parameter));
    //        asteroidSigma = 183.269E18 * convertAUperdayToSI;
    //        varianceVector.push_back(asteroidSigma*asteroidSigma);

    //        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                 ("2000003", gravitational_parameter));
    //        asteroidSigma = 126.860E18 * convertAUperdayToSI;
    //        varianceVector.push_back(asteroidSigma*asteroidSigma);

    //        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                 ("2000004", gravitational_parameter));
    //        asteroidSigma = 93.970E18 * convertAUperdayToSI;
    //        varianceVector.push_back(asteroidSigma*asteroidSigma);

    //        bool gammaIsEstimated = false;
    //        bool betaIsEstimated = false;
    //        bool nordtvedtParameterIsEstimated = false;

    //        // relativistic parameters
    //        if (calculateSchwarzschildCorrection == true
    //            || includeSEPViolationAcceleration == true){

    //            if ( gammaIsAConsiderParameter == false){
    //                parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                         ("global_metric", ppn_parameter_gamma ) );
    //                varianceVector.push_back(sigmaGamma*sigmaGamma);
    //                gammaIsEstimated = true;
    //            }


    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("global_metric", ppn_parameter_beta ) );
    //            varianceVector.push_back(sigmaBeta*sigmaBeta);
    //            betaIsEstimated = true;
    //        }

    //        // Nordtvedt parameter
    //        int nordtvedtParameterIndex = -1;
    //        if (includeSEPViolationAcceleration &&
    //                (nordtvedtConstraintTrueOrFalse == false || ignoreNordtvedtConstraintInEstimation == false)){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("global_metric", ppn_nordtvedt_parameter ) );
    //            varianceVector.push_back(sigmaNordtvedt*sigmaNordtvedt);
    //            nordtvedtParameterIsEstimated = true;
    //            nordtvedtParameterIndex = varianceVector.size()-1;
    //        }

    //        if (estimatePPNalphas == true){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("global_metric", ppn_parameter_alpha1 ) );
    //            varianceVector.push_back(sigmaAlpha1*sigmaAlpha1);

    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("global_metric", ppn_parameter_alpha2 ) );
    //            varianceVector.push_back(sigmaAlpha2*sigmaAlpha2);
    //        }

    //        // angular momentum
    //        if (calculateLenseThirringCorrection && estimateSunAngularMomentum){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", angular_momentum));
    //            varianceVector.push_back(sigmaSunAngularMomentum*sigmaSunAngularMomentum);
    //        }

    //        // time varying gravitational parameter
    //        if (includeTVGPAcceleration == true){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("global_metric", time_varying_gravitational_parameter));
    //            varianceVector.push_back(sigmaTVGP*sigmaTVGP);
    //        }

    //        // gravity field Sun (mu and spherical harmonics)
    //        parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                 ("Sun", gravitational_parameter));
    //        varianceVector.push_back(sigmaSunGM*sigmaSunGM);

    //        // time varying J2 Sun
    //        double variableJ2parametersVariance = 1.0; //sigma taken as a percentage of the mean
    //        if (estimateJ2Amplitude){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J2_amplitude));
    //            double sigmaJ2Amplitude = variableJ2parametersVariance * relativity::variableJ2Interface->getAmplitude();
    //            varianceVector.push_back(sigmaJ2Amplitude*sigmaJ2Amplitude);
    //        }
    //        if (estimateJ2Period){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J2_period));
    //            double sigmaJ2Period = variableJ2parametersVariance * relativity::variableJ2Interface->getPeriod();
    //            varianceVector.push_back(sigmaJ2Period*sigmaJ2Period);
    //        }
    //        if (estimateJ2Phase){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J2_phase));
    //            double sigmaJ2Phase = variableJ2parametersVariance * relativity::variableJ2Interface->getPhase();
    //            varianceVector.push_back(sigmaJ2Phase*sigmaJ2Phase);
    //        }

    //        // time varying J4 Sun
    //        double variableJ4parametersVariance = 1.0; //sigma taken as a percentage of the mean
    //        if (estimateJ4Amplitude){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J4_amplitude));
    //            double sigmaJ4Amplitude = variableJ4parametersVariance * relativity::variableJ4Interface->getAmplitude();
    //            varianceVector.push_back(sigmaJ4Amplitude*sigmaJ4Amplitude);
    //        }
    //        if (estimateJ4Period){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J4_period));
    //            double sigmaJ4Period = variableJ4parametersVariance * relativity::variableJ4Interface->getPeriod();
    //            varianceVector.push_back(sigmaJ4Period*sigmaJ4Period);
    //        }
    //        if (estimateJ4Phase){
    //            parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                     ("Sun", variable_J4_phase));
    //            double sigmaJ4Phase = variableJ4parametersVariance * relativity::variableJ4Interface->getPhase();
    //            varianceVector.push_back(sigmaJ4Phase*sigmaJ4Phase);
    //        }


    //        // conventional spherical harmonics (always needs to be last estimatable parameter)
    //        std::vector< std::pair< int, int > > blockIndices;
    //        for (unsigned int d=2; d<=maximumDegreeSunGravitationalMoment; d+=2){
    //            blockIndices.push_back(std::make_pair(d,0));
    //            varianceVector.push_back(sigmaValuesSunGravitationalMoments(d)*sigmaValuesSunGravitationalMoments(d));
    //        }
    //        parameterNames.push_back(std::make_shared<SphericalHarmonicEstimatableParameterSettings>
    //                                 (blockIndices,"Sun",spherical_harmonics_cosine_coefficient_block));


            // Create parameters
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                    createParametersToEstimate( parameterNames, bodyMap );

            // Print identifiers and indices of parameters to terminal
            printEstimatableParameterEntries( parametersToEstimate );

            // If gamma, beta, eta, are estimatable parameters and nordtvedt constraint is set to true, enforce it in the estimation
            bool enforceNordtvedtConstraintInEstimation = false;
    //        if ( (gammaIsEstimated || gammaIsAConsiderParameter )
    //            && betaIsEstimated && nordtvedtParameterIsEstimated
    //            && nordtvedtConstraintTrueOrFalse
    //            && ignoreNordtvedtConstraintInEstimation == false){
    //            enforceNordtvedtConstraintInEstimation = true;
    //        }


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
            std::vector< double > baseTimeList =
                    makeObservationTimeList(secondsAfterJ2000(observationInitialTime),
                                            secondsAfterJ2000(observationFinalTime),
                                            observationTimeStep,
                                            trackingArcDuration,
                                            maximumNumberOfTrackingDays,
                                            unit_conversions::convertDegreesToRadians(maxMSEAngleDeg),
                                            flybyList);
            int baseTimeListSize = baseTimeList.size();

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

    //        // Create noise functions per observable
    //        std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;

    //        std::function< double( const double ) > mercuryOrbiterNoiseFunction;
    //        mercuryOrbiterNoiseFunction = [noiseAtMinAngle, noiseAtMaxAngle, maxMSEAngleDeg](const double time){
    //            return noiseSampleBasedOnMSEangle(time, noiseAtMinAngle, noiseAtMaxAngle, maxMSEAngleDeg);
    //        };

    //        noiseFunctions[ one_way_range ] = mercuryOrbiterNoiseFunction;
    //        noiseFunctions[ n_way_range ] = mercuryOrbiterNoiseFunction;




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


    //        // retrieve observation times and measurements from output of main application
    //        Eigen::VectorXd observationMeasurements = input_output::readMatrixFromFile(inputSubFolder + "/ObservationMeasurements.dat");
    //        Eigen::VectorXd observationTimesEigenVector = input_output::readMatrixFromFile(inputSubFolder + "/ObservationTimes.dat");
    //        std::vector< double > observationTimes;
    //        for (unsigned int i=0; i<observationTimesEigenVector.size(); i++){
    //            observationTimes.push_back(observationTimesEigenVector(i));
    //        }

    //        PodInputDataType observationsAndTimes;
    //        observationsAndTimes[n_way_range][twoWayRangeLinkEnds] =
    //                std::make_pair(observationMeasurements,
    //                               std::make_pair(observationTimes, receiver));

            // Simulate observations
            PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                        measurementSimulationInput,
                        orbitDeterminationManager.getObservationSimulators( ));

    //        // similar container, but the "observation" will be the noise value instead of the actual observation
    //        PodInputDataType observationWeightsAndTimes = observationsAndTimes;

    //        // add spacecraft initial position error to the observations
    //        std::map< double, Eigen::Vector3d > interpolatedErrorMatrix;
    //        if (includeSpacecraftPositionError == true){

    //            std::cout << "Adding satellite estimation initial position error..." << std::endl;
    //    //        Eigen::Vector3d constantSatelliteError; constantSatelliteError << 10.0, 10.0, 10.0;

    //            // get interpolated error maps for list of vehicle observation times
    //            std::string vehicleErrorFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/error_inputs_"+vehicle+".txt";
    //            Eigen::MatrixXd error_input = input_output::readMatrixFromFile(vehicleErrorFilename, ",");

    //            interpolatedErrorMatrix =
    //                    interpolatePositionErrorsBasedOnTrueAnomaly(
    //                        error_input, baseTimeList, vehicle, mercuryGravitationalParameter);

    //            int interpolatedErrorMatrixSize = interpolatedErrorMatrix.size();
    //            std::cout << "Size interpolated vehicle error map: "<< interpolatedErrorMatrixSize << std::endl;
    //            if (baseTimeListSize != interpolatedErrorMatrixSize){
    //                std::runtime_error("observation lists unequal!");
    //            }

    //            std::random_device rd;
    //            std::mt19937 gen(rd());

    //            // get to the location in the map where we can find the range observables
    //            PodInputDataType::iterator podInputIterator = observationsAndTimes.begin();
    //            PodInputDataType::iterator weightsIterator = observationWeightsAndTimes.begin();
    //            while (podInputIterator != observationsAndTimes.end()){

    //                if (podInputIterator->first == one_way_range || podInputIterator->first == n_way_range){

    //                    SingleObservablePodInputType::iterator singleObservableIterator = podInputIterator->second.begin();
    //                    SingleObservablePodInputType::iterator weightsIterator2 = weightsIterator->second.begin();
    //                    while (singleObservableIterator != podInputIterator->second.end()){

    //                        // retrieve observations and their respective times for current observable type
    //                        ObservationVectorType allObservations = singleObservableIterator->second.first;
    //                        std::vector< double > allObservationTimes = singleObservableIterator->second.second.first;

    //                        ObservationVectorType newObservations = Eigen::VectorXd(allObservations.size());
    //                        ObservationVectorType observationWeights = Eigen::VectorXd(allObservations.size());

    //                        // for every observation, retrieve and add the range bias that should be added
    //                        for (unsigned int i=0; i<allObservationTimes.size(); i++){

    //                            double observationTime = allObservationTimes.at( i );

    //                            double rangeNoiseLevel = noiseLevelBasedOnMSEangle(
    //                                        observationTime, noiseAtMinAngle, noiseAtMaxAngle, maxMSEAngleDeg);

    //                            double totalErrorLevel = combinedRangeAndSatelliteErrorLevel(
    //                                        observationTime, interpolatedErrorMatrix.at(observationTime), rangeNoiseLevel);

    //                            std::normal_distribution<double> d(0.0, totalErrorLevel);
    //                            double rangeErrorSample = d(gen);

    //                            if (podInputIterator->first == n_way_range){rangeErrorSample *= 2.0;}

    //                            if (rangeErrorSample > 10000.0){
    //                                std::cout<<"ERROR: unusually large range correction at time "<<observationTime<<std::endl;
    //                                std::cout<<" total error level: "<<totalErrorLevel<<std::endl;
    //                                std::cout<<" range correction (from random sample): "<<rangeErrorSample<<std::endl;
    //                            }

    //                            newObservations(i) = allObservations(i) + rangeErrorSample;
    //                            observationWeights(i) = 1.0/(totalErrorLevel*totalErrorLevel);

    //    //                        std::cout<<observationTime<<" // "<<
    //    //                                   currentSatelliteError.transpose()<<" // "<<
    //    //                                   randomErrorSample.transpose()<<" // "<<
    //    //                                   rangeCorrection<<" // "<<
    //    //                                   noiseLevel<<std::endl;
    //                        }


    ////                        for (unsigned int i=0; i<allObservationTimes.size(); i++){

    ////                            double observationTime = allObservationTimes.at( i );

    ////                            Eigen::Vector3d currentSatelliteError = interpolatedErrorMatrix.at(observationTime);
    ////                            Eigen::Vector3d randomErrorSample;
    ////                            for (unsigned int j=0; j<3; j++){
    ////                                std::normal_distribution<double> d(0.0, currentSatelliteError( j ));
    ////                                randomErrorSample( j ) = d(gen);
    ////                            }

    ////                            Eigen::Vector3d mercuryPositionWrtEarth = -getBodyCartesianStateAtEpoch("Earth","Mercury","IAU_Mercury","None",observationTime).segment(0,3);
    ////                            Eigen::Vector3d rangeUnitVector = mercuryPositionWrtEarth / mercuryPositionWrtEarth.norm( );

    ////                            double rangeCorrection = randomErrorSample.dot(rangeUnitVector);
    ////                            if (podInputIterator->first == n_way_range){rangeCorrection *= 2.0;}

    ////                            if (rangeCorrection > 10000.0){
    ////                                std::cout<<"ERROR: unusually large range correction at time "<<observationTime<<std::endl;
    ////                                std::cout<<" interpolated vehicle error: "<<currentSatelliteError.transpose()<<std::endl;
    ////                                std::cout<<" range correction (from random sample): "<<rangeCorrection<<std::endl;
    ////                            }

    ////                            newObservations(i) = allObservations(i) + rangeCorrection;

    ////                            double noiseLevel = abs(rangeCorrection) + noiseLevelBasedOnMSEangle(
    ////                                        observationTime, noiseAtMinAngle, noiseAtMaxAngle, maxMSEAngleDeg);

    ////                            observationWeights(i) = 1.0/(noiseLevel*noiseLevel);

    ////    //                        std::cout<<observationTime<<" // "<<
    ////    //                                   currentSatelliteError.transpose()<<" // "<<
    ////    //                                   randomErrorSample.transpose()<<" // "<<
    ////    //                                   rangeCorrection<<" // "<<
    ////    //                                   noiseLevel<<std::endl;
    ////                        }


    //                        singleObservableIterator->second.first = newObservations;
    //                        weightsIterator2->second.first = observationWeights;
    //                        singleObservableIterator++; weightsIterator2++;
    //                    }
    //                }
    //                podInputIterator++; weightsIterator++;
    //            }
    //        }



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

    //        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ ){
    //            // perturb body positions by  1 meters
    //            parameterPerturbation.segment(i*6,3) = Eigen::Vector3d::Constant( 1.0 );
    //            // perturb body velocities by 0.001 m/s
    //            parameterPerturbation.segment(i*6+3,3) = Eigen::Vector3d::Constant( 0.001 );
    //        }

            // perturb eta slightly to prevent partials from being 0 due to which the estimation is unable to estimate eta
    //        if (nordtvedtParameterIndex >= 0){
    //            parameterPerturbation(nordtvedtParameterIndex) = 1.0E-10;
    //        }

            initialParameterEstimate += parameterPerturbation;

            std::cout << "True parameter values:" << std::endl;
            std::cout << truthParameters.transpose() << std::endl;
            std::cout << "Parameter perturbations:" << std::endl;
            std::cout << parameterPerturbation.transpose() << std::endl;
            std::cout << "Initial guesses:" << std::endl;
            std::cout << initialParameterEstimate.transpose() << std::endl;


            // Define a priori covariance matrix
            Eigen::MatrixXd aprioriCovariance =
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));


            Eigen::VectorXd sigmaEigenVector(truthParameters.rows());
            for( unsigned int i = 0; i < truthParameters.size(); i++ ){
                aprioriCovariance( i,i ) = varianceVector.at( i );
                sigmaEigenVector(i) = sqrt(varianceVector.at( i ));
            }
            std::cout << "a priori covariance matrix:" << std::endl
                      << aprioriCovariance.diagonal().transpose() << std::endl;
            Eigen::MatrixXd inverseOfAprioriCovariance = aprioriCovariance.inverse();

            // Define estimation input
            std::shared_ptr< PodInput< double, double > > podInput =
                    std::make_shared< PodInput< double, double > >(
                        observationsAndTimes, initialParameterEstimate.rows( ),
                        inverseOfAprioriCovariance,
                        initialParameterEstimate - truthParameters );
            podInput->defineEstimationSettings( true, true, true, true, true, true,
                                                enforceNordtvedtConstraintInEstimation );
    //        podInput->manuallySetObservationWeightMatrix(observationWeightsAndTimes);


    //        // save stuff related to observations
    //        input_output::writeDataMapToTextFile( interpolatedErrorMatrix, "interpolatedErrorMatrix.dat", outputSubFolder );
    //        input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ), "ObservationMeasurements.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
    //                                         "ObservationTimes.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
    //                                         "ObservationLinkEnds.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
    //                                         "ObservationObservableTypes.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
    //                                         "ObservationMeasurements.dat", 16, outputSubFolder );

            // Perform estimation
            std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                        podInput, std::make_shared< EstimationConvergenceChecker >( maximumNumberOfIterations ) );

            // Print true estimation error, limited mostly by numerical error
            Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
            Eigen::VectorXd formalError = podOutput->getFormalErrorVector( );

            std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
            std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;
            std::cout << "True to form estimation error ratio is: " << std::endl <<
                         ( podOutput->getFormalErrorVector( ).cwiseQuotient( estimationError ) ).transpose( ) << std::endl;


            /////////////////////////////
            //// CONSIDER PARAMETERS ////
            /////////////////////////////

            // get other required matrices for the calculation
    //        Eigen::MatrixXd initialCovarianceMatrix = podOutput->getUnnormalizedCovarianceMatrix( ); // P
    //        Eigen::VectorXd observationWeightDiagonal = podOutput->weightsMatrixDiagonal_; // diagonal of W
            Eigen::MatrixXd partialDerivativesWrtAsteroidGravitationalParameters
                    = podOutput->getUnnormalizedPartialDerivatives( );

            // save partial derivatives, seperately for every asteroid
            for (unsigned int a=0; a<maxAsteroid; a++){
                // check if key exists in map (i.e. asteroid is included in file)
                if (asteroidsINPOP19a.find(asteroidNumbers(a)) != asteroidsINPOP19a.end()){
                    // add SPICE number of asteroid to bodies
                    Eigen::VectorXd partialsCurrentAsteroid = partialDerivativesWrtAsteroidGravitationalParameters.col(a+6);
                    input_output::writeMatrixToFile( partialsCurrentAsteroid,
                                                     "partialDerivativesWrtAsteroid"+std::to_string(asteroidNumbers(a))+".dat",
                                                     16, outputSubFolder );

                }
            }


            input_output::writeMatrixToFile( sigmaEigenVector,
                                             "aprioriUncertaintiesAsteroidGravitationalParameters.dat", 16, outputSubFolder );

    //        Eigen::MatrixXd considerCovarianceMatrix;
    //        unsigned int maxcovtype = 1;

    //        if (gammaIsAConsiderParameter || ppnAlphasAreConsiderParameters){

    //            std::cout<< "calculating covariance due to consider parameters..."<< std::endl;

    //            // In order to get the partial derivatives of consider parameters wrt observations, run an estimation of the consider parameters
    //            std::vector< std::shared_ptr < EstimatableParameterSettings > > considerParameterNames;
    //            std::vector<double> considerVarianceVector;

    //            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    //            {
    //                int j = 6*i;
    //                // bodies must be included for the estimation to work
    //                considerParameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
    //                                              bodiesToPropagate[i], systemInitialState.segment(j,6), centralBodies[i] ) );

    //                for( unsigned int i = 0; i < 3; i++ ){
    //                    considerVarianceVector.push_back(sigmaPosition*sigmaPosition);
    //                }
    //                for( unsigned int i = 3; i < 6; i++ ){
    //                    considerVarianceVector.push_back(sigmaVelocity*sigmaVelocity);
    //                }
    //            }

    //            if ( gammaIsAConsiderParameter == true ){
    //                considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                         ("global_metric", ppn_parameter_gamma ) );
    //                considerVarianceVector.push_back(sigmaGamma*sigmaGamma);
    //            }

    //            if ( ppnAlphasAreConsiderParameters){
    //                considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                         ("global_metric", ppn_parameter_alpha1 ) );
    //                considerVarianceVector.push_back(sigmaAlpha1*sigmaAlpha1);

    //                considerParameterNames.push_back(std::make_shared<EstimatableParameterSettings >
    //                                         ("global_metric", ppn_parameter_alpha2 ) );
    //                considerVarianceVector.push_back(sigmaAlpha2*sigmaAlpha2);
    //            }

    //            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > considerParametersToEstimate =
    //                    createParametersToEstimate( considerParameterNames, bodyMap );
    //            printEstimatableParameterEntries( considerParametersToEstimate );

    //            OrbitDeterminationManager< double, double > orbitDeterminationManagerConsiderParameters =
    //                    OrbitDeterminationManager< double, double >(
    //                        bodyMap, considerParametersToEstimate, observationSettingsMap,
    //                        integratorSettings, propagatorSettings );

    //            Eigen::Matrix< double, Eigen::Dynamic, 1 > initialConsiderParameterEstimate =
    //                    considerParametersToEstimate->template getFullParameterValues< double >( );
    //            Eigen::Matrix< double, Eigen::Dynamic, 1 > truthConsiderParameters = initialConsiderParameterEstimate;

    //            Eigen::MatrixXd considerParameterAprioriCovariance = // C
    //                Eigen::MatrixXd::Zero( truthConsiderParameters.rows( ), truthConsiderParameters.rows( ));

    //            for( unsigned int i = 0; i < truthConsiderParameters.size(); i++ ){
    //                considerParameterAprioriCovariance( i,i ) = considerVarianceVector.at( i );
    //            }
    //            std::cout << "consider parameter a priori covariance matrix:" << std::endl
    //                      << considerParameterAprioriCovariance.diagonal().transpose() << std::endl;

    //            std::shared_ptr< PodInput< double, double > > podInputConsiderParameters =
    //                    std::make_shared< PodInput< double, double > >(
    //                        observationsAndTimes, initialConsiderParameterEstimate.rows( ),
    //                        considerParameterAprioriCovariance.inverse(),
    //                        initialConsiderParameterEstimate - truthConsiderParameters );
    //            podInputConsiderParameters->defineEstimationSettings( true, true, true, true );
    //            podInputConsiderParameters->manuallySetObservationWeightMatrix(observationWeightsAndTimes);

    //            std::shared_ptr< PodOutput< double > > podOutputConsiderParameters = orbitDeterminationManagerConsiderParameters.estimateParameters(
    //                        podInputConsiderParameters, std::make_shared< EstimationConvergenceChecker >( 1 ) );

    //            Eigen::MatrixXd partialDerivativesOfConsiderParameters = // H_c
    //                    (podOutputConsiderParameters->getUnnormalizedPartialDerivatives()).rightCols(
    //                        truthConsiderParameters.size()-6*numberOfNumericalBodies);

    //            // calculate covariance matrix including contribution from consider parameters
    //            considerCovarianceMatrix = calculateConsiderCovarianceMatrix(
    //                        initialCovarianceMatrix, observationWeightDiagonal,
    //                        considerParameterAprioriCovariance.block(6,6,considerParameterAprioriCovariance.rows()-6, considerParameterAprioriCovariance.cols()-6),
    //                        unnormalizedPartialDerivatives, partialDerivativesOfConsiderParameters);

    //            // asteroids


    //            // get consider correlation matrix
    //            Eigen::VectorXd formalErrorWithConsiderParameters = considerCovarianceMatrix.diagonal( ).cwiseSqrt( );
    //            Eigen::MatrixXd considerCorrelationMatrix = considerCovarianceMatrix.cwiseQuotient(
    //                        formalErrorWithConsiderParameters * formalErrorWithConsiderParameters.transpose() );

    //            input_output::writeMatrixToFile( considerParameterAprioriCovariance, "ConsiderParameterAprioriCovariance.dat", 16, outputSubFolder );
    //            input_output::writeMatrixToFile( considerCorrelationMatrix, "EstimationConsiderCorrelations.dat", 16, outputSubFolder );
    //            input_output::writeMatrixToFile( formalErrorWithConsiderParameters, "ObservationFormalEstimationErrorWithConsiderParameters.dat", 16, outputSubFolder );

    //            maxcovtype = 2;
    //        }


    //        /////////////////////////////////////////////
    //        //// PROVIDE OUTPUT TO CONSOLE AND FILES ////
    //        /////////////////////////////////////////////

    //        std::cout<< "provide output..." << std::endl;

    //        // propagate according to integration history. earlier result is separated here as the times are needed on their own.
    //        std::vector<double> fullStateHistoryTimes;
    //        std::map<double, Eigen::VectorXd> propagationHistory = allBodiesPropagationHistory.at( 0 );
    //        std::map<double, Eigen::VectorXd>::iterator historyit = propagationHistory.begin();

    //        while (historyit != propagationHistory.end()){
    //            fullStateHistoryTimes.push_back(historyit->first);
    //            historyit++;
    //        }

    //        // Propagate covariance matrix (twice, both with and without the consider parameters included in the covariance)

    //        for (unsigned int covtype = 0; covtype<maxcovtype; covtype++){

    //            Eigen::MatrixXd initialCovariance;
    //            std::string saveString;
    //            if (covtype == 0){
    //                std::cout<<"propagating covariance matrix..."<<std::endl;
    //                initialCovariance = initialCovarianceMatrix;
    //                saveString = "";
    //            } else{
    //                std::cout<<"propagating consider covariance matrix..."<<std::endl;
    //                initialCovariance = considerCovarianceMatrix;
    //                saveString = "Consider";
    //            }

    //            std::map< double, Eigen::MatrixXd > propagatedCovariance;
    //            propagateCovariance(propagatedCovariance, initialCovariance,
    //                                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
    //                                fullStateHistoryTimes);
    //                                //baseTimeList
    //                                //60.0, initialEphemerisTime + 3600.0, finalEphemerisTime - 3600.0 );

    //            std::map<double, double> trueAnomalyMap;
    //            std::map<double, Eigen::Vector6d> propagatedErrorUsingCovMatrix;
    //            std::map<double, Eigen::Vector6d> propagatedRSWErrorUsingCovMatrix;

    //            std::map<double, Eigen::MatrixXd>::iterator it = propagatedCovariance.begin();

    //            while (it != propagatedCovariance.end()){

    //                Eigen::Vector6d currentCartesianState = propagationHistory.at(it->first);

    //                // save true anomaly of MESSENGER around Mercury at this time step
    //                Eigen::Vector6d currentKeplerianState = convertCartesianToKeplerianElements(
    //                            currentCartesianState, sunGravitationalParameter);
    //                double currentTrueAnomaly = currentKeplerianState(5);

    //                trueAnomalyMap.insert(std::make_pair( it->first, currentTrueAnomaly ) );

    //                // save propagated error (using covariance) in Cartesian and RSW frame
    //                Eigen::MatrixXd currentCovariance = it->second;
    //                Eigen::Matrix3d currentTransformationToRSW =
    //                        reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(currentCartesianState);

    //                Eigen::Vector6d currentSigmaVector = currentCovariance.diagonal( ).cwiseSqrt( );

    //                Eigen::Vector6d currentSigmaVectorRSW;
    //                currentSigmaVectorRSW.segment(0,3) =
    //                        ( currentTransformationToRSW * currentCovariance.block(0,0,3,3) * currentTransformationToRSW.transpose( ) ).diagonal( ).cwiseSqrt( );
    //                currentSigmaVectorRSW.segment(3,3) =
    //                        ( currentTransformationToRSW * currentCovariance.block(3,3,3,3) * currentTransformationToRSW.transpose( ) ).diagonal( ).cwiseSqrt( );

    //                propagatedErrorUsingCovMatrix.insert(std::make_pair( it->first, currentSigmaVector ) );
    //                propagatedRSWErrorUsingCovMatrix.insert(std::make_pair( it->first, currentSigmaVectorRSW ) );

    //                it++;
    //            }

    //            input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedCovariance, onlyEveryXthValue), "Propagated"+saveString+"Covariance.dat", outputSubFolder );
    //            input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedErrorUsingCovMatrix, onlyEveryXthValue), "propagatedErrorUsing"+saveString+"CovMatrix.dat", outputSubFolder );
    //            input_output::writeDataMapToTextFile( onlyEveryXthValueFromDataMap(propagatedRSWErrorUsingCovMatrix, onlyEveryXthValue), "propagatedRSWErrorUsing"+saveString+"CovMatrix.dat", outputSubFolder );
    //        }


    //        // Save data in files
    //        std::cout<<"writing output to files..."<<std::endl;

    //        // errors
    //        input_output::writeMatrixToFile( truthParameters, "TruthParameters.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( estimationError, "ObservationTrueEstimationError.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( formalError, "ObservationFormalEstimationError.dat", 16, outputSubFolder );

    //        // residuals, parameters
    //        input_output::writeMatrixToFile( podOutput->residuals_, "EstimationResiduals.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ), "ResidualHistory.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ), "ParameterHistory.dat", 16, outputSubFolder );

    //        // covariance and correlation
    //        input_output::writeMatrixToFile( considerCovarianceMatrix, "ConsiderCovarianceMatrix.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( initialCovarianceMatrix, "InitialCovarianceMatrix.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( observationWeightDiagonal, "ObservationWeightDiagonal.dat", 16, outputSubFolder );
    //        input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ), "EstimationCorrelations.dat", 16, outputSubFolder );

    ////        input_output::writeMatrixToFile( unnormalizedPartialDerivatives, "test_unnormalizedPartialDerivatives.dat", 16, outputSubFolder );
    ////        input_output::writeMatrixToFile( partialDerivativesOfConsiderParameters, "test_partialDerivativesOfConsiderParameters.dat", 16, outputSubFolder );

    ////        input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_, "EstimationInformationMatrix.dat", 16, outputSubFolder );
    ////        input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_, "EstimationInformationMatrixNormalization.dat", 16, outputSubFolder );
    ////        input_output::writeMatrixToFile( podInput->getInverseOfAprioriCovariance( ), "InverseAPrioriCovariance.dat", 16, outputSubFolder );
    ////        input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_, "InverseNormalizedCovariance.dat", 16, outputSubFolder );
    ////        input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ), "UnnormalizedCovariance.dat", 16, outputSubFolder );


    //        // dependent variables
    //        input_output::writeDataMapToTextFile(
    //                    onlyEveryXthValueFromDataMap(podOutput->dependentVariableHistoryFinalIteration_.at(0), onlyEveryXthValue), "DependentVariablesHistoryFinalIteration.dat",
    //                    outputSubFolder, "", std::numeric_limits< double >::digits10, std::numeric_limits< double >::digits10, "," );

            std::cout << "done!" << std::endl;
        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;

}
