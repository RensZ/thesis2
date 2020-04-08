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

//#include <SatellitePropagatorExamples/applicationOutput.h>

// Get path for output directory.
namespace tudat_applications
{
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


    //a-priori sigma values
    const double sigmaPosition = 1000.0;
    const double sigmaVelocity = 1.0;
    const double sigmaRadiation = 0.5;
    const double sigmaGamma = 2.3E-5; //Genova 2018
    const double sigmaSunJ2 = 0.03E-7; //Genova 2018
    std::vector<double> varianceVector;


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

    // Specify initial and final time
    double initialEphemerisTime = 410000000.0; //september 2012
    int numberOfSimulationDays = 5.0;
    double finalEphemerisTime = initialEphemerisTime + numberOfSimulationDays * 86400.0;

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 2*86400.0, finalEphemerisTime + 2*86400.0);

//    bodySettings[ "Sun" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
//                "ECLIPJ2000", "IAU_Sun", spice_interface::computeRotationQuaternionBetweenFrames(
//                    "ECLIPJ2000", "IAU_Sun", initialEphemerisTime ),
//                initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY );

//    bodySettings[ "Mercury" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
//                "ECLIPJ2000", "IAU_Mercury", spice_interface::computeRotationQuaternionBetweenFrames(
//                    "ECLIPJ2000", "IAU_Mercury", initialEphemerisTime ),
//                initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY );

    // todo: fix rotation rate above

    // Custom settings Sun
    const double sunRadius = 695.7E6; //m, from nasa fact sheet
    const double sunJ2 = 2.20E-7; //from Genova 2018
    const double sunGravitationalParameter = 132712440041.9394E9; //m3/s2, from Genova 2018

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

    const double mercuryRadius =  2.4400000000000000E6; //m
    const double mercuryGravitationalParameter = 2.2031863566000000E13; //m3/s2

    const unsigned int maxMercuryDegree = 6;
    const unsigned int maxMercuryOrder = 6;

    Eigen::MatrixXd mercuryCosineCoefficients;
    Eigen::MatrixXd mercurySineCoefficients;

    mercuryCosineCoefficients = Eigen::MatrixXd::Zero(maxMercuryDegree+1,maxMercuryOrder+1);
    mercurySineCoefficients   = Eigen::MatrixXd::Zero(maxMercuryDegree+1,maxMercuryOrder+1);

    mercuryCosineCoefficients(0,0) = 1.0;

    std::string HgM008File = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/HgM008.txt";
    Eigen::MatrixXd HgM008 =
             tudat::input_output::readMatrixFromFile( HgM008File , ",", "#" );
    Eigen::MatrixXi HgM008i =
             tudat::input_output::readMatrixFromFile( HgM008File , ",", "#" ).cast<int>();

    std::cout << HgM008 << std::endl;


    for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){

        unsigned int d = HgM008i(i,0);
        unsigned int o = HgM008i(i,1);
        if (d > maxMercuryDegree || o > maxMercuryOrder){
            break;
        }

        std::cout << "d: " << d << "    o: " << o << "    Snm: " << HgM008(i,2) << "    Cnm: " << HgM008(i,3) << std::endl;
        double normalization = calculateLegendreGeodesyNormalizationFactor(d,o);
        mercuryCosineCoefficients(d,o) = HgM008(i,2)/normalization;
        mercurySineCoefficients(d,o) = HgM008(i,3)/normalization;
    }

    std::cout<< mercuryCosineCoefficients << std::endl;
    std::cout<< mercurySineCoefficients << std::endl;

    bodySettings[ "Mercury" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                mercuryGravitationalParameter, mercuryRadius,
                mercuryCosineCoefficients, mercurySineCoefficients, "IAU_Mercury" );




    NamedBodyMap bodyMap = createBodies( bodySettings );
    bodyMap[ "Vehicle" ] = std::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 1000.0 );

//    // Create aerodynamic coefficient interface settings.
//    double referenceArea = 4.0;
//    double aerodynamicCoefficient = 1.2;
//    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
//            std::make_shared< ConstantAerodynamicCoefficientSettings >(
//                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.01, 0.1 ).finished( ), 1, 1 );

//    // Create and set aerodynamic coefficients object
//    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
//                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 3.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mercury" );
    std::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Vehicle" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    vehicleRadiationPressureSettings, "Vehicle", bodyMap ) );

    bodyMap[ "Vehicle" ]->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ), "Mercury", "ECLIPJ2000" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS               //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodyMap.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -1.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.8, 4.0 ).finished( ), geodetic_position );

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

//    accelerationsOfVehicle[ "Mercury" ].push_back( std::make_shared< EmpiricalAccelerationSettings >( ) );

    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                     true, false, false ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    std::cout<<"accelerations on Vehicle set"<<std::endl;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;

    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Mercury" );


    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );

    std::cout<<"acceleration models created"<<std::endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//    Eigen::Vector6d vehicleInitialKeplerianState;
//    vehicleInitialKeplerianState( semiMajorAxisIndex ) = (5.172314175394792E-05)*physical_constants::ASTRONOMICAL_UNIT;
//    vehicleInitialKeplerianState( eccentricityIndex ) = 6.330051519045469E-01;
//    vehicleInitialKeplerianState( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 8.811002404818556E+01 );
//    vehicleInitialKeplerianState( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
//    vehicleInitialKeplerianState( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
//    vehicleInitialKeplerianState( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );


//    // MESSENGER state at 2012-9-4 from JPL HORIZON
//    Eigen::Vector6d vehicleInitialStateKM;
//    vehicleInitialStateKM <<  -3.903241060109104E+03, 1.132080847551607E+03,  3.115447461584540E+02
//                              -1.543646767332147E+00, 1.748791694677803E-01, -2.355773564131267E+00;

//    Eigen::Vector6d vehicleInitialState = (10^3)*vehicleInitialStateKM;

//    Eigen::Vector6d vehicleInitialKeplerianState = convertCartesianToKeplerianElements(vehicleInitialState,mercuryGravitationalParameter);
//    std::cout<<vehicleInitialKeplerianState<<std::endl;


    loadSpiceKernelInTudat(input_output::getSpiceKernelPath() + "msgr_040803_130501_140822_od318sc_0.bsp");

    // Define arc length
    double arcDuration = 1.01 * 86400.0;
    double arcOverlap = 3600.0;



    // Create propagator settings (including initial state taken from Kepler orbit) for each arc
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;

    std::vector< double > arcStartTimes;
    double currentTime = initialEphemerisTime;
    while( currentTime <= finalEphemerisTime )
    {
        arcStartTimes.push_back( currentTime );

        Eigen::Vector6d currentArcInitialState =
                getBodyCartesianStateAtEpoch("MESSENGER","Mercury","ECLIPJ2000","None",currentTime);



//        Eigen::Vector6d currentArcInitialState = convertKeplerianToCartesianElements(
//                    propagateKeplerOrbit( vehicleInitialKeplerianState, currentTime - initialEphemerisTime,
//                                          mercuryGravitationalParameter ), mercuryGravitationalParameter );

        std::cout<<currentArcInitialState<<std::endl;

        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                                              centralBodies, accelerationModelMap, bodiesToIntegrate, currentArcInitialState,
                                              currentTime + arcDuration + arcOverlap ) );

        currentTime += arcDuration;
    }



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
                double( initialEphemerisTime ), 15.0,
                RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                1.0, 30.0, 1.0E-12, 1.0E-12); //tolerances iets lager?

//        std::shared_ptr< IntegratorSettings< > > integratorSettings =
//                std::make_shared< IntegratorSettings< > >
//                ( rungeKutta4, initialEphemerisTime, 2.0 );



    ////////////////////////////
    //// DYNAMICS SIMULATOR ////
    ////////////////////////////

    std::cout << "running dynamics simulator..." << std::endl;

    MultiArcDynamicsSimulator <> dynamicsSimulator (bodyMap,
                                                    integratorSettings,
                                                    propagatorSettings,
                                                    arcStartTimes,
                                                    true,false,true);

    std::cout << "saving integration result and dependent variables..." << std::endl;

    std::vector< std::map< double, Eigen::VectorXd > > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Retrieve numerically integrated state for each body.
    std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
    allBodiesPropagationHistory.resize( bodiesToIntegrate.size() );

    for( unsigned int i = 0; i < arcStartTimes.size(); i++ ){

        std::map< double, Eigen::VectorXd > intermediateIntegrationResult = integrationResult.at(i);

        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = intermediateIntegrationResult.begin( );
             stateIterator != intermediateIntegrationResult.end( ); stateIterator++ )
        {
            for( unsigned int i = 0; i < bodiesToIntegrate.size(); i++ )
            {
                allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
            }
        }

    }

    for( unsigned int i = 0; i < bodiesToIntegrate.size(); i++ )
    {
        // Write propagation history to file.
        input_output::writeDataMapToTextFile(
                    allBodiesPropagationHistory[ i ],
                    "StatePropagationHistory" + bodiesToIntegrate.at( i ) + ".dat",
                     tudat_applications::getOutputPath( ) + "MercuryOrbiterStateEstimation/",
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );
    }




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
//    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    // one-way
    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );

    // two-way
    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 1 ] );

    //    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 2 ] );
    //    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );


//    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
//    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

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
                                  "Vehicle", systemInitialState, arcStartTimes, "Mercury" ) );

//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "global_metric", ppn_parameter_gamma ) );
//    varianceVector.push_back( sigmaGamma*sigmaGamma );

    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
    varianceVector.push_back( sigmaRadiation*sigmaRadiation );

//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", constant_rotation_rate ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station1" ) );
//    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station2" ) );
//    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
//                                  linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, true ) );
//    parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
//                                  linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, true ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 0, "Sun", spherical_harmonics_cosine_coefficient_block ) );
    varianceVector.push_back( sigmaSunJ2*sigmaSunJ2 );

    //    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
    //                                  2, 0, maxMercuryDegree, maxMercuryOrder, "Mercury", spherical_harmonics_cosine_coefficient_block ) );
    //    for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){
    //        unsigned int d = HgM008i(i,0);
    //        unsigned int o = HgM008i(i,1);
    //        if (d > maxMercuryDegree || o > maxMercuryOrder){
    //            break;
    //        }
    //        if (d >= 2){
    //            double normalization = calculateLegendreGeodesyNormalizationFactor(d,o);
    //            double sigma = HgM008(i,4)/normalization;
    //            varianceVector.push_back( sigma );
    //        }
    //    }

    //    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
    //                                  2, 1, maxMercuryDegree, maxMercuryOrder, "Mercury", spherical_harmonics_sine_coefficient_block ) );
    //    for( unsigned int i = 0; i < HgM008.col(0).size(); i++ ){
    //        unsigned int d = HgM008i(i,0);
    //        unsigned int o = HgM008i(i,1);
    //        if (d > maxMercuryDegree || o > maxMercuryOrder){
    //            break;
    //        }
    //        if (d >= 2 && o >= 1){
    //            double normalization = calculateLegendreGeodesyNormalizationFactor(d,o);
    //            double sigma = HgM008(i,5)/normalization;
    //            varianceVector.push_back( sigma );
    //        }
    //    }


//    // Define required settings for arc-wise empirical accelerations
//    std::map< EmpiricalAccelerationComponents, std::vector< EmpiricalAccelerationFunctionalShapes > > empiricalAccelerationComponents;
//    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( cosine_empirical );
//    empiricalAccelerationComponents[ across_track_empirical_acceleration_component ].push_back( sine_empirical );
//    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( cosine_empirical );
//    empiricalAccelerationComponents[ along_track_empirical_acceleration_component ].push_back( sine_empirical );
//    std::vector< double > empiricalAccelerationArcTimes;
//    empiricalAccelerationArcTimes.push_back( initialEphemerisTime );
//    empiricalAccelerationArcTimes.push_back( initialEphemerisTime + ( finalEphemerisTime - initialEphemerisTime ) / 2.0 );
//    parameterNames.push_back( std::make_shared< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >(
//                                  "Vehicle", "Earth", empiricalAccelerationComponents, empiricalAccelerationArcTimes ) );

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

            // Define settings for observable, no light-time corrections, and biases for selected links
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

    // Define time of first observation
    double observationTimeStart = initialEphemerisTime + 1000.0;

    // Define time between two observations
    double  observationInterval = 60.0;

    // Simulate observations for each day in simulation
    std::vector< double > baseTimeList;
    for( int i = 0; i < numberOfSimulationDays; i++ )
    {
        // Simulate 480 observations per day (observationInterval apart)
        for( unsigned int j = 0; j < 480; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
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

    // Create observation viability settings and calculators
//    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
//    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
//                                                minimum_elevation_angle, std::make_pair( "Earth", "" ), "",
//                                                unit_conversions::convertDegreesToRadians( 5.0 ) ) );
//    PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
//                bodyMap, linkEndsPerObservable, observationViabilitySettings );

    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Define noise levels
    double rangeNoise = 0.1;
    double angularPositionNoise = 1.0E-7;
    double dopplerNoise = (0.1E-3)/physical_constants::SPEED_OF_LIGHT; // Mazarico et al. 2014: 0.1mm/s at 60s integration time
    std::cout << "doppler noise: " << dopplerNoise << std::endl;

    // Defaults
//    double rangeNoise = 0.1;
//    double angularPositionNoise = 1.0E-7;
//    double dopplerNoise = 1.0E-12;

    // Create noise functions per observable
    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;
    noiseFunctions[ one_way_range ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, rangeNoise }, 0.0 ), std::placeholders::_1 );

    noiseFunctions[ angular_position ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, angularPositionNoise }, 0.0 ), std::placeholders::_1 );

    noiseFunctions[ one_way_doppler ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, dopplerNoise }, 0.0 ), std::placeholders::_1 );


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
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        parameterPerturbation.segment( 6*i, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 6*i+3, 3 ) = Eigen::Vector3d::Constant( 1.0E-3 );
    }
    initialParameterEstimate += parameterPerturbation;


    // Define a priori covariance matrix
    Eigen::MatrixXd aprioriCovariance =
        Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));

    std::cout << "a priori variances:" << std::endl;
    for( unsigned int i = 0; i < truthParameters.size(); i++ ){
        std::cout << varianceVector.at( i ) << " ";
        aprioriCovariance( i,i ) = varianceVector.at( i );
    }
    std::cout << std::endl;

    Eigen::MatrixXd inverseOfAprioriCovariance = aprioriCovariance.inverse();


    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                inverseOfAprioriCovariance,
                initialParameterEstimate - truthParameters );


    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( rangeNoise * rangeNoise );
    weightPerObservable[ angular_position ] = 1.0 / ( angularPositionNoise * angularPositionNoise );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( dopplerNoise * dopplerNoise );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->defineEstimationSettings( true, false, true, true, true );



    // Perform estimation
    const unsigned int numberOfIterations = 5;
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout<< "provide output..." << std::endl;

    std::string outputSubFolder = "OutputMercuryOrbiter/";

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;

    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;
    std::cout << "True to form estimation error ratio is: " << std::endl <<
                 ( podOutput->getFormalErrorVector( ).cwiseQuotient( estimationError ) ).transpose( ) << std::endl;

    // Propagate errors to total time period
    std::map< double, Eigen::VectorXd > propagatedErrors;
    propagateFormalErrors(
                propagatedErrors, podOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                60.0, initialEphemerisTime + 3600.0, finalEphemerisTime - 3600.0 );
    input_output::writeDataMapToTextFile( propagatedErrors,
                                          "PropagatedErrors.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "EstimationInformationMatrix.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( podOutput->informationMatrixTransformationDiagonal_,
                                     "EstimationInformationMatrixNormalization.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

//    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
//                                     "EstimationWeightsDiagonal.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );
//    input_output::writeMatrixToFile( podOutput->residuals_,
//                                     "EstimationResiduals.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "EstimationCorrelations.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "ResidualHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "ParameterHistory.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

//    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
//                                     "ObservationMeasurements.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationTimes.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
//                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
//                                     "ObservationLinkEnds.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
//                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
//                                     "ObservationObservableTypes.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );
//    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
//                                     "ObservationMeasurements.dat", 16,
//                                     tudat_applications::getOutputPath( ) + outputSubFolder );


    input_output::writeMatrixToFile( estimationError,
                                     "ObservationTrueEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    Eigen::VectorXd formalError = podOutput->getFormalErrorVector( );
    input_output::writeMatrixToFile( formalError,
                                     "ObservationFormalEstimationError.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );

    // get initial state errors in the RSW frame
    Eigen::Vector6d averageTrueError, averageFormalError;
    Eigen::Vector6d averageTrueErrorRSW, averageFormalErrorRSW;

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

//        std::cout<<formalError.segment( 6*i, 6 ).transpose()<<std::endl;
//        std::cout<<formalError.segment( 6*i, 6 ).cwiseAbs().transpose()<<std::endl;
//        std::cout<<(formalError.segment( 6*i, 6 ).cwiseAbs()/numberOfSimulationDays).transpose()<<std::endl;

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


    input_output::writeMatrixToFile( trueEstimationErrorRSW,
                                     "ObservationTrueEstimationErrorRSW.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );
    input_output::writeMatrixToFile( formalEstimationErrorRSW,
                                     "ObservationFormalEstimationErrorRSW.dat", 16,
                                     tudat_applications::getOutputPath( ) + outputSubFolder );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
