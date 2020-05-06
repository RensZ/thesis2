   /*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

// custom functions written for the main application are placed here to save space:
#include "tudatApplications/thesis/MyApplications/customFunctions.h"

// Get path for output directory.
namespace tudat_applications
{
    static inline std::string getOutputPath(
            const std::string& extraDirectory = "" )
    {
        std::string outputPath = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Output/" + extraDirectory;
        return outputPath;
    }
}



//! Execute propagation of orbits of Asterix and Obelix around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    using namespace tudat::input_output;


    //body input parameters
    const double sunRadius = 695.7E6; //m, from nasa fact sheet
    const double sunJ2 = 2.20E-7; //from Genova 2018
    const double sunGravitationalParameter = 132712440041.9394E9; //m3/s2, from Genova 2018

    const double mercuryRadius = 2.44E6; //m, from https://pgda.gsfc.nasa.gov/products/71
    const double mercuryGravitationalParameter = (2.2031870798779644e+04)*(1E9); //m3/s2, from https://pgda.gsfc.nasa.gov/products/71

    //a-priori sigma values
    const double sigmaPosition = 1000.0;
    const double sigmaVelocity = 1.0;
    const double sigmaRadiation = 0.5;
    const double sigmaMercuryGM = (8.6276900083604273e-05)*(1E9);
    const double sigmaGamma = 2.3E-5; //Genova 2018
    const double sigmaSunJ2 = 0.25E-7; //Genova 2018
    std::vector<double> varianceVector;

    // integrator settings
    double initialStepSize = 2.0;
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = 10.0;
    double tolerance = 1.0E-15;
//    double initialStepSize = 2.0;
//    double minimumStepSize = 2.0;
//    double maximumStepSize = 2.0;
//    double tolerance = 1.0;

    const unsigned int maxMercuryDegree = 8;
    const unsigned int maxMercuryOrder = 8;

    const unsigned int minMercuryDegree = 2; //only for estimatable parameters, SH field starts at d/o 0/0

    // estimation settings
    const unsigned int numberOfIterations = 1;

    // run simulation for vehicle
    std::string vehicle = "MESSENGER";
//    std::string vehicle = "BepiColombo";

    // define mission dependent variables
    double initialEphemerisTime;
    double finalEphemerisTime;
    double spacecraftMass;
    double referenceAreaRadiation;
    std::string vehicleKernel;
    std::string vehicleName;
    double trackingPeriod;
    double observationInterval;
    double dopplerNoiseUnnormalised;


    if (vehicle == "BepiColombo"){
        initialEphemerisTime = 828273600.0; //April 1st 2026, 16 days after final orbit insertion
        finalEphemerisTime = 891432000.0; //April 1st 2028, 1 month before end of extended misson
        spacecraftMass = 1000.0;
        referenceAreaRadiation = 6.0;
        vehicleKernel = "bc_mpo_mlt_50037_20260314_20280529_v01.bsp";
        vehicleName = "BEPICOLOMBO MPO";
        trackingPeriod = 8.0*60.0*60.0;
        observationInterval = 1000.0;
        dopplerNoiseUnnormalised = 1.5E-6; //Ka tracking only
    }
    else if (vehicle == "MESSENGER"){
        initialEphemerisTime = 354888000.0; // April 1st 2011, 12 days after insertion
        finalEphemerisTime = 478440000.0; // March 1st 2015, 2 months before crash
        spacecraftMass = 1000.0;
        referenceAreaRadiation = 3.0;
        vehicleKernel = "msgr_040803_150430_150430_od431sc_2.bsp";
        vehicleName = "MESSENGER";
        trackingPeriod = 8.0*60.0*60.0;
        observationInterval = 60.0;
        dopplerNoiseUnnormalised = 0.1E-3;
    }
    else{
        std::cout<<"ERROR: Vehicle name not recognized, try a different input"<<std::endl;
    }


    std::string outputSubFolder = vehicle + "/";

    int numberOfSimulationDays = 30.0;
    double arcOverlap = 0.0;
    double observationStartOffset = 1000.0; // or observation generation wil complain
    //    double arcDuration = 1.01 * 86400.0; // integrate for one day
    double arcDuration = trackingPeriod + 2.0*observationStartOffset; // integrate for the duration of the observations
    double arcOffset = (finalEphemerisTime-initialEphemerisTime)/(numberOfSimulationDays-1.0);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

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
    bodyNames[ 7 ] = "Moon";

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 2*86400.0, finalEphemerisTime + 2*86400.0);


    double dayOnMercury = 58.785 * physical_constants::JULIAN_DAY; //https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
    bodySettings[ "Mercury" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Mercury", spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Mercury", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI / dayOnMercury );

    // todo: fix rotation rate above

    // Custom settings Sun
    double sunNormalizedJ2 = sunJ2 / calculateLegendreGeodesyNormalizationFactor(2,0);

    Eigen::MatrixXd normalizedSineCoefficients;
    Eigen::MatrixXd normalizedCosineCoefficients;

    normalizedSineCoefficients   = Eigen::MatrixXd::Zero(3,3);
    normalizedCosineCoefficients = Eigen::MatrixXd::Zero(3,3);

    normalizedCosineCoefficients(0,0) = 1.0;
    normalizedCosineCoefficients(2,0) = sunNormalizedJ2;

    bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                sunGravitationalParameter, sunRadius,
                normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );

    std::cout<< normalizedCosineCoefficients << std::endl;
    std::cout<< normalizedSineCoefficients << std::endl;


    // Mercury SH gravity field from https://pgda.gsfc.nasa.gov/products/71
    // Alternative option from https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/
    Eigen::MatrixXd mercuryCosineCoefficients;
    Eigen::MatrixXd mercurySineCoefficients;

    mercuryCosineCoefficients = Eigen::MatrixXd::Zero(maxMercuryDegree+1,maxMercuryOrder+1);
    mercurySineCoefficients   = Eigen::MatrixXd::Zero(maxMercuryDegree+1,maxMercuryOrder+1);

    mercuryCosineCoefficients(0,0) = 1.0;

    //std::string HgM008File = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/MercurySH/HgM008.txt";
    std::string HgM008File = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/MercurySH/HgM005.txt";

    Eigen::MatrixXd HgM008 =
             tudat::input_output::readMatrixFromFile( HgM008File , ",", "#", 1 );
    Eigen::MatrixXi HgM008i =
             tudat::input_output::readMatrixFromFile( HgM008File , ",", "#", 1 ).cast<int>();
//    std::cout << HgM008 << std::endl;

    for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){

        unsigned int d = HgM008i(i,0);
        unsigned int o = HgM008i(i,1);
        if (d > maxMercuryDegree || o > maxMercuryOrder){
            break;
        }

        std::cout << "d/o " << d << "/" << o << " Cnm: " << HgM008(i,2) << " Snm: " << HgM008(i,3) << std::endl;
        double normalization = 1.0; // calculateLegendreGeodesyNormalizationFactor(d,o);
        mercuryCosineCoefficients(d,o) = HgM008(i,2)/normalization;
        mercurySineCoefficients(d,o) = HgM008(i,3)/normalization;
    }

    std::cout<< mercuryCosineCoefficients << std::endl;
    std::cout<< mercurySineCoefficients << std::endl;

    bodySettings[ "Mercury" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                mercuryGravitationalParameter, mercuryRadius,
                mercuryCosineCoefficients, mercurySineCoefficients, "IAU_Mercury" );


    NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ vehicle ] = std::make_shared< Body >( );
    bodyMap[ vehicle ]->setConstantBodyMass( spacecraftMass );

    // Create radiation pressure settings
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mercury" );
    std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ vehicle ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    vehicleRadiationPressureSettings, vehicle, bodyMap ) );

    bodyMap[ vehicle ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ), "Mercury", "ECLIPJ2000" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;


    //    accelerationsOfVehicle[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
    //                                                   basic_astrodynamics::central_gravity ) );

    accelerationsOfVehicle[ "Mercury" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( maxMercuryDegree, maxMercuryOrder ) );

    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );

    accelerationsOfVehicle[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );

    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                     true, false, false ) );

    //    accelerationsOfVehicle[ "Mercury" ].push_back( std::make_shared< EmpiricalAccelerationSettings >( ) );

    accelerationMap[ vehicle ] = accelerationsOfVehicle;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    bodiesToIntegrate.push_back( vehicle );
    centralBodies.push_back( "Mercury" );


    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::cout<<"creating propagation settings"<<std::endl;

    if (vehicle == "BepiColombo"){
        loadSpiceKernelInTudat(input_output::getSpiceKernelPath() + "bc_sci_v04.tf");
    }
    loadSpiceKernelInTudat(input_output::getSpiceKernelPath() + vehicleKernel);


    // Create propagator settings (including initial state taken from Kepler orbit) for each arc
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;

    std::vector< double > arcStartTimes;
    std::vector< double > arcEndTimes;
    std::vector< Eigen::VectorXd > arcInitialStates;
    double currentTime = initialEphemerisTime;

    while( currentTime <= finalEphemerisTime )
    {

        double currentMSEAngleDegrees = unit_conversions::convertRadiansToDegrees(
                    angleBetween2Bodies(currentTime,"Sun","Mercury","Earth"));

        if (currentMSEAngleDegrees < 145.0){

            arcStartTimes.push_back( currentTime );

            Eigen::Vector6d currentArcInitialState = getBodyCartesianStateAtEpoch(vehicleName,"Mercury","ECLIPJ2000","None",currentTime);

            Eigen::Vector6d currentKeplerianState = convertCartesianToKeplerianElements(currentArcInitialState, mercuryGravitationalParameter);

//            currentKeplerianState( 5 ) = unit_conversions::convertDegreesToRadians(0.0);
//            currentArcInitialState = convertKeplerianToCartesianElements(currentKeplerianState, mercuryGravitationalParameter);

            std::cout<<"a: "<<currentKeplerianState(0);
            std::cout<<" e: "<<currentKeplerianState(1);
            std::cout<<" i: "<<unit_conversions::convertRadiansToDegrees(currentKeplerianState(2));
            std::cout<<" AoP: "<<unit_conversions::convertRadiansToDegrees(currentKeplerianState(3));
            std::cout<<" LoAN: "<<unit_conversions::convertRadiansToDegrees(currentKeplerianState(4));
            std::cout<<" TA: "<<unit_conversions::convertRadiansToDegrees(currentKeplerianState(5));
            std::cout<<" MSE: "<<currentMSEAngleDegrees<<std::endl;

            double currentArcEndTime = currentTime + arcDuration + arcOverlap;
            arcEndTimes.push_back( currentArcEndTime );

            arcInitialStates.push_back( currentArcInitialState );

            std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                  std::make_shared< propagators::PropagationTimeTerminationSettings >( currentArcEndTime, true );

            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                                                  centralBodies, accelerationModelMap, bodiesToIntegrate, currentArcInitialState,
                                                   terminationSettings ) );
        } else{
            std::cout<<"arc at time "<<currentTime<<" not included, MSE: "<<currentMSEAngleDegrees <<std::endl;
        }

        currentTime += arcOffset;
    }

    numberOfSimulationDays = static_cast<double>(arcStartTimes.size());


    // Create propagator settings
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );

    // Create integrator settings
//    std::shared_ptr< AdamsBashforthMoultonSettings< double > > integratorSettings =
//            std::make_shared< AdamsBashforthMoultonSettings< double > > (
//                initialEphemerisTime, 8.0,
//                1.0, 32.0,
//                1E-12, 1E-12,
//                6, 12);


    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                double( initialEphemerisTime ), initialStepSize,
                RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                minimumStepSize, maximumStepSize,
                tolerance, tolerance);

//        std::shared_ptr< IntegratorSettings< > > integratorSettings =
//                std::make_shared< IntegratorSettings< > >
//                ( rungeKutta4, initialEphemerisTime, 1.0 );




    ////////////////////////////
    //// DYNAMICS SIMULATOR ////
    ////////////////////////////

    std::cout << "running dynamics simulator..." << std::endl;

    MultiArcDynamicsSimulator <> dynamicsSimulator (bodyMap,
                                                    integratorSettings,
                                                    propagatorSettings,
                                                    arcStartTimes,
                                                    true,false,true);

    std::cout << "saving integration result..." << std::endl;

    std::vector< std::map< double, Eigen::VectorXd > > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    Eigen::MatrixXd finalStates(arcStartTimes.size(),6);

    // Retrieve numerically integrated state for each body.
    std::map< double, Eigen::VectorXd > propagationHistory;

    for( unsigned int i = 0; i < arcStartTimes.size(); i++ ){

        std::map< double, Eigen::VectorXd > intermediateIntegrationResult = integrationResult.at( i );

        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
             stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
        {
            propagationHistory[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
            if (stateIterator->first == arcEndTimes.at(i)){
                finalStates.row(i) = stateIterator->second.segment( 0, 6 );
            }
        }
    }

    // Write propagation history to file.
    input_output::writeDataMapToTextFile(
                propagationHistory,
                "StatePropagationHistoryVehicle.dat",
                 tudat_applications::getOutputPath( ) + outputSubFolder,
                "",
                std::numeric_limits< double >::digits10,
                std::numeric_limits< double >::digits10,
                "," );



    //////////////////////////////
    //// BACKWARD PROPAGATION ////
    //////////////////////////////

    std::cout << "performing backward propagation..." << std::endl;

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > backwardPropagatorSettingsList;

    for (unsigned i=0; i<arcStartTimes.size(); i++){

        std::shared_ptr< PropagationTimeTerminationSettings > backwardTerminationSettings =
              std::make_shared< propagators::PropagationTimeTerminationSettings >( arcStartTimes.at(i), true );

        backwardPropagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                                              centralBodies, accelerationModelMap, bodiesToIntegrate, finalStates.row(i),
                                               backwardTerminationSettings ) );
    }


    std::shared_ptr< PropagatorSettings< double > > backwardPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( backwardPropagatorSettingsList );

    std::shared_ptr< IntegratorSettings< double > > backwardIntegratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                double( finalEphemerisTime ), -initialStepSize,
                RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                -minimumStepSize, -maximumStepSize,
                tolerance, tolerance);

    MultiArcDynamicsSimulator <> backwardDynamicsSimulator (bodyMap,
                                                    backwardIntegratorSettings,
                                                    backwardPropagatorSettings,
                                                    arcEndTimes,
                                                    true,false,true);


    std::cout << "saving backward integration result..." << std::endl;

    std::vector< std::map< double, Eigen::VectorXd > > backwardIntegrationResult = backwardDynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Retrieve numerically integrated state for each body.
    std::map< double, Eigen::VectorXd > backwardPropagationHistory;
    Eigen::MatrixXd initialStatesFromBackwardPropagation(arcStartTimes.size(),6);
    Eigen::MatrixXd arcInitialStatesMatrix(arcStartTimes.size(),6);


    for( unsigned int i = 0; i < arcStartTimes.size(); i++ ){

        std::map< double, Eigen::VectorXd > intermediateIntegrationResult = backwardIntegrationResult.at( i );

        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
             stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
        {
            backwardPropagationHistory[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
            if (stateIterator->first == arcStartTimes.at(i)){
                initialStatesFromBackwardPropagation.row(i) = stateIterator->second.segment( 0, 6 );
                arcInitialStatesMatrix.row(i) = arcInitialStates.at(i);
            }
        }
    }

    Eigen::MatrixXd initialStatesIntegrationError = initialStatesFromBackwardPropagation - arcInitialStatesMatrix;
    std::cout<<"differences initial states: "<<std::endl<<initialStatesIntegrationError<<std::endl;


    // Write propagation history to file.
    input_output::writeDataMapToTextFile(
                backwardPropagationHistory,
                "StatePropagationHistoryVehicleBackwards.dat",
                 tudat_applications::getOutputPath( ) + outputSubFolder,
                "",
                std::numeric_limits< double >::digits10,
                std::numeric_limits< double >::digits10,
                "," );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;


    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    linkEnds[ receiver ] = std::make_pair( vehicle, "" );
    stationTransmitterLinkEnds.push_back( linkEnds );

    linkEnds.clear( );
    linkEnds[ receiver ] = std::make_pair( "Earth", "" );
    linkEnds[ transmitter ] = std::make_pair( vehicle, "" );
    stationReceiverLinkEnds.push_back( linkEnds );

    LinkEnds twoWayDopplerLinkEnds;
    twoWayDopplerLinkEnds[ transmitter ] = std::make_pair( "Earth", "" );
    twoWayDopplerLinkEnds[ reflector1 ] = std::make_pair( vehicle, "" );
    twoWayDopplerLinkEnds[ receiver ] = std::make_pair( "Earth", "" );


    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "define list of parameters to estimate..." << std::endl;

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;

    // Create concatenated list of arc initial states
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 6 * arcStartTimes.size( ) );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        systemInitialState.segment( i * 6, 6 ) = propagatorSettingsList.at( i )->getInitialStates( );

        varianceVector.push_back(sigmaPosition*sigmaPosition);
        varianceVector.push_back(sigmaPosition*sigmaPosition);
        varianceVector.push_back(sigmaPosition*sigmaPosition);
        varianceVector.push_back(sigmaVelocity*sigmaVelocity);
        varianceVector.push_back(sigmaVelocity*sigmaVelocity);
        varianceVector.push_back(sigmaVelocity*sigmaVelocity);

    }

    parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                  vehicle, systemInitialState, arcStartTimes, "Mercury" ) );



    if (maxMercuryDegree >= minMercuryDegree){
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      minMercuryDegree, 0, maxMercuryDegree, maxMercuryOrder, "Mercury", spherical_harmonics_cosine_coefficient_block ) );
        for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){
            unsigned int d = HgM008i(i,0);
            unsigned int o = HgM008i(i,1);
            if (d >= minMercuryDegree && o >= 0 && d <= maxMercuryDegree && o <= maxMercuryOrder){
                double sigma = HgM008(i,4);
                std::cout<<"Cnm d/o "<<d<<"/"<<o<<" sigma: "<<sigma<<std::endl;
                varianceVector.push_back( sigma*sigma );
            }
        }

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      minMercuryDegree, 1, maxMercuryDegree, maxMercuryOrder, "Mercury", spherical_harmonics_sine_coefficient_block ) );
        for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){
            unsigned int d = HgM008i(i,0);
            unsigned int o = HgM008i(i,1);
            if (d >= minMercuryDegree && o >= 1 && d <= maxMercuryDegree && o <= maxMercuryOrder){
                double sigma = HgM008(i,5);
                std::cout<<"Snm d/o "<<d<<"/"<<o<<" sigma: "<<sigma<<std::endl;
                varianceVector.push_back( sigma*sigma );
            }
        }
    }

    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( vehicle, radiation_pressure_coefficient ) );
    varianceVector.push_back( sigmaRadiation*sigmaRadiation );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap, propagatorSettings );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE OBSERVATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "create observation settings..." << std::endl;

    // Iterate over all observable types and associated link ends, and creating settings for observation
    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define bias and light-time correction settings
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections;
            if( currentObservable == one_way_range )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Zero( ), true );
            }
            observationSettingsMap.insert(
                        std::make_pair( currentLinkEndsList.at( i ),
                                        std::make_shared< ObservationSettings >(
                                            currentObservable, lightTimeCorrections, biasSettings ) ) );
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          INITIALIZE ORBIT DETERMINATION OBJECT     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "create OD object.." << std::endl;

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////          SIMULATE OBSERVATIONS                     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "observation simulation settings..." << std::endl;

    // Generate observation time settings based on inputs (observations start at initial state + 1000s)
    std::vector< double > baseTimeList;
    for( int i = 0; i < numberOfSimulationDays; i++ )
    {

        double observationTimeStart = arcStartTimes.at(i) + observationStartOffset;
        int numberOfObservations = static_cast<int>(trackingPeriod / observationInterval);

        for( int j = 0; j < numberOfObservations; j++ )
        {
            baseTimeList.push_back( observationTimeStart +
                                    static_cast< double >( j ) * observationInterval );
        }

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
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >( receiver, baseTimeList );
        }
    }


    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Define noise levels
    double rangeNoise = 1.0;
    double angularPositionNoise = 1.0E-7;
    double twoWayDopplerNoise = (dopplerNoiseUnnormalised)/physical_constants::SPEED_OF_LIGHT; // Mazarico et al. 2014: 0.1mm/s at 60s integration time
    double oneWayDopplerNoise = twoWayDopplerNoise/sqrt(2.0);
    std::cout << "one way doppler noise: " << oneWayDopplerNoise << std::endl;
    std::cout << "one way doppler noise: " << twoWayDopplerNoise << std::endl;


    // Create noise functions per observable
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;
    noiseFunctions[ one_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, oneWayDopplerNoise }, 0.0 ), std::placeholders::_1 );

    noiseFunctions[ two_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, twoWayDopplerNoise }, 0.0 ), std::placeholders::_1 );


    noiseFunctions[ one_way_range ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, rangeNoise }, 0.0 ), std::placeholders::_1 );


    std::cout<< "simulating observations..." << std::endl;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservationsWithNoise< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), noiseFunctions);



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "estimate parameters..." << std::endl;

    // Perturb parameter estimate
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    bool reintegrateOnFirstIteration = false;
    if (numberOfIterations > 1){
        for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
        {
            parameterPerturbation.segment( 6*i, 3 ) = Eigen::Vector3d::Constant( 1.0 );
            parameterPerturbation.segment( 6*i+3, 3 ) = Eigen::Vector3d::Constant( 1.0E-3 );
        }
        initialParameterEstimate += parameterPerturbation;
        reintegrateOnFirstIteration = true;
    }


    std::cout << "truth parameters:" << truthParameters.transpose() << std::endl;
    std::cout << "added as initial guess:" << parameterPerturbation.transpose() << std::endl;

    std::cout<<"number of parameters: "<<truthParameters.size()
             <<", number of variance entries: "<<varianceVector.size()<<std::endl;

    // Define a priori covariance matrix
    if (truthParameters.size() != static_cast<long>(varianceVector.size())){
        std::cout<<("WARNING: parameters and variance vectors are not the same size")<<std::endl;
    }

    Eigen::MatrixXd aprioriCovariance =
        Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));
    Eigen::MatrixXd inverseOfAprioriCovariance;

    std::cout << "a priori variances:" << std::endl;
    for( unsigned int i = 0; i < truthParameters.size(); i++ ){
        std::cout << varianceVector.at( i ) << " ";
        aprioriCovariance( i,i ) = varianceVector.at( i );
    }
    std::cout << std::endl;

    inverseOfAprioriCovariance = aprioriCovariance.inverse();



    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes,
                initialParameterEstimate.rows( ),
                inverseOfAprioriCovariance,
                initialParameterEstimate - truthParameters );


    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( rangeNoise * rangeNoise );
    weightPerObservable[ angular_position ] = 1.0 / ( angularPositionNoise * angularPositionNoise );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( oneWayDopplerNoise * oneWayDopplerNoise );
    weightPerObservable[ two_way_doppler ] = 1.0 / ( twoWayDopplerNoise * twoWayDopplerNoise );
    weightPerObservable[ one_way_differenced_range ] = 1.0 / ( twoWayDopplerNoise * twoWayDopplerNoise );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->defineEstimationSettings( reintegrateOnFirstIteration, //reintegrate on first iteration
                                        true, //reintegrate variational equations
                                        true, //save information matrix
                                        true, //print output
                                        true ); //save state history


    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "provide output..." << std::endl;

    // propagate according to integration history. earlier result is separated here as the times are needed on their own.
    std::vector<double> fullStateHistoryTimes;
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
                    currentCartesianState, mercuryGravitationalParameter);
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

    // print initial state errors in the RSW frame
    Eigen::Vector6d averageTrueError = Eigen::Vector6d::Zero();
    Eigen::Vector6d averageFormalError = Eigen::Vector6d::Zero();
    Eigen::Vector6d averageTrueErrorRSW = Eigen::Vector6d::Zero();
    Eigen::Vector6d averageFormalErrorRSW = Eigen::Vector6d::Zero();

    Eigen::VectorXd trueEstimationErrorRSW = Eigen::VectorXd::Zero(6*numberOfSimulationDays);
    Eigen::VectorXd formalEstimationErrorRSW = Eigen::VectorXd::Zero(6*numberOfSimulationDays);

    for (int i=0; i<numberOfSimulationDays; i++){

        Eigen::Matrix3d transformationToRSW =
                    reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                        propagatorSettingsList.at( i )->getInitialStates( ));

        trueEstimationErrorRSW.segment  ( 6*i  , 3 ) = transformationToRSW * estimationError.segment( 6*i  , 3 );
        trueEstimationErrorRSW.segment  ( 6*i+3, 3 ) = transformationToRSW * estimationError.segment( 6*i+3, 3 );

        formalEstimationErrorRSW.segment( 6*i  , 3 ) = transformationToRSW * formalError.segment( 6*i  , 3 );
        formalEstimationErrorRSW.segment( 6*i+3, 3 ) = transformationToRSW * formalError.segment( 6*i+3, 3 );

        averageTrueError      += estimationError.segment( 6*i, 6 ).cwiseAbs()/static_cast< double >(numberOfSimulationDays);
        averageFormalError    += formalError.segment( 6*i, 6 ).cwiseAbs()/static_cast< double >(numberOfSimulationDays);
        averageTrueErrorRSW   += trueEstimationErrorRSW.segment( 6*i, 6 ).cwiseAbs()/static_cast< double >(numberOfSimulationDays);
        averageFormalErrorRSW += formalEstimationErrorRSW.segment( 6*i, 6 ).cwiseAbs()/static_cast< double >(numberOfSimulationDays);
    }

    std::cout<<"average true errors of all initial states in Cartesian frame (absolute values):"<<std::endl
             <<averageTrueError.transpose()<<std::endl;

    std::cout<<"average formal errors of all initial states in Cartesian frame (absolute values):"<<std::endl
             <<averageFormalError.transpose()<<std::endl;

    std::cout<<"true to formal ratio:"<<std::endl
             <<averageTrueError.cwiseQuotient(averageFormalError).transpose()<<std::endl;

    std::cout<<"average true errors of all initial states in RSW frame (absolute values):"<<std::endl
             <<averageTrueErrorRSW.transpose()<<std::endl;

    std::cout<<"average formal errors of all initial states in RSW frame (absolute values):"<<std::endl
             <<averageFormalErrorRSW.transpose()<<std::endl;

    std::cout<<"true to formal ratio:"<<std::endl
             <<averageTrueErrorRSW.cwiseQuotient(averageFormalErrorRSW).transpose()<<std::endl;


    // Save data in files
    std::cout<<"writing output to files..."<<std::endl;

    input_output::writeDataMapToTextFile( trueAnomalyMap,
                                          "TrueAnomaly.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedCovariance,
                                          "PropagatedCovariance.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedErrorUsingCovMatrix,
                                          "propagatedErrorUsingCovMatrix.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedRSWErrorUsingCovMatrix,
                                          "propagatedRSWErrorUsingCovMatrix.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "EstimationInformationMatrix.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "EstimationInformationMatrixNormalization.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "EstimationWeightsDiagonal.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "EstimationResiduals.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "EstimationCorrelations.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "ResidualHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "ParameterHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "ObservationMeasurements.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationTimes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "ObservationLinkEnds.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationObservableTypes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "ObservationMeasurements.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( estimationError,
                                     "ObservationTrueEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( formalError,
                                     "ObservationFormalEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( trueEstimationErrorRSW,
                                     "ObservationTrueEstimationErrorRSW.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( formalEstimationErrorRSW,
                                     "ObservationFormalEstimationErrorRSW.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podInput->getInverseOfAprioriCovariance( ),
                                     "InverseAPrioriCovariance.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->inverseNormalizedCovarianceMatrix_,
                                     "InverseNormalizedCovariance.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ),
                                     "UnnormalizedCovariance.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    Eigen::MatrixXd partialDerivatives = podOutput->getUnnormalizedPartialDerivatives();
    input_output::writeMatrixToFile( partialDerivatives,
                                     "UnnormalizedPartialDerivatives.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );


    std::cout<<"done!"<<std::endl;
    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
