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

#include <tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h>

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


// convert calendar date and time to seconds after J2000
double secondsAfterJ2000(Eigen::Vector6i datetime){
    using namespace tudat;
    using namespace tudat::basic_astrodynamics;
    const unsigned int secondsPerDay = 60*60*24;
    const unsigned int julianDayJ2000 = 2451545;
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay(datetime[0],datetime[1],datetime[2],datetime[3],datetime[4],datetime[5]);
    return (julianDay - julianDayJ2000)*secondsPerDay;
}


// Generate a list of observations given some basic inputs
std::vector<double> makeObservationTimeList(double initialTime,
                                            double endTime,
                                            double timeStep,
                                            std::vector<double> flybyTimes){
    std::vector< double > observationTimeList = flybyTimes;
    double currentTime = initialTime;
    while (currentTime <= endTime){
        observationTimeList.push_back(currentTime);
        currentTime += timeStep;
    }
    return observationTimeList;
}




int main( )
{

    // Define namespaces
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

    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );



    /////////////////////
    //// USER INPUTS ////
    /////////////////////

    // Specify start and end time simulation
    Eigen::Vector6i initialTime, finalTime;
    initialTime << 2008, 1, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    finalTime   << 2015, 5, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss

    // Acceleration settings
    const bool calculateSchwarzschildCorrection = true;
    const bool calculateLenseThirringCorrection = false;
    const bool calculateDeSitterCorrection = false;
    const bool estimateJ4 = false;

    // Parameter inputs
    const double sunRadius = 695.7E6; //m, from nasa fact sheet
    const double sunJ2 = 2.20E-7; //from Genova 2018
    const double sunJ4 = -4.34E-9; //from Antia 2008
    const double sunGravitationalParameter = 132712440041.9394E9; //m3/s2, from Genova 2018
    const double sunAngularMomentum = 190.0E39; //kgm2/s, from Pijpers 1998
    double timeVaryingGravitationalParameter = -1E-13; //PLACEHOLDER

    // Parameter apriori values and uncertainties
    const bool useAprioriValues = true;

//    const double aprioriGamma = 1.0;
//    const double aprioriBeta = 1.0;
//    const double aprioriSunGM = sunGravitationalParameter;
//    const double aprioriSunJ2 = sunJ2;
//    const double aprioriSunJ4 = sunJ4;
//    Eigen::VectorXd aprioriParameters(aprioriGamma, aprioriBeta, aprioriSunGM, aprioriSunJ2);

    const double sigmaPosition = 10.0; //educated guess
    const double sigmaVelocity = 10.0E-6; //educated guess
    const double sigmaGamma = 2.3E-5; //Genova 2018
    const double sigmaBeta = 6.9E-5; //Genova 2018
    const double sigmaSunGM = 0.14E9; //Genova 2018
    const double sigmaSunJ2 = 0.03E-7; //Genova 2018
    const double sigmaSunJ4 = 0.1E-9; //PLACEHOLDER
    const double sigmaTVGP = 1E14; //PLACEHOLDER

    // Planet propagation settings
    const bool propogatePlanets = false; // Propogate the other planets besides Mercury (NOTE: need observations for other planets, or LS can't find solutions for other planets)
    const bool useDirectSpice = false; // Use direct SPICE (more accurate than tabulated Spice)

    // Observation settings
    const bool simulateObservationNoise = true;
    double rangeNoise = 2.0;
    double angularPositionNoise = 1.0E-7;
    double dopplerNoise = 1.0E-12;
    double positionObservableNoise = 5.0;

    // MESSENGER observation schedule
    Eigen::Vector6i observationInitialTime, observationFinalTime;
    observationInitialTime << 2011, 3, 18, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    observationFinalTime   << 2015, 4, 28, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    double observationTimeStep = 60.0*60.0*24.0; // seconds

    // MESSENGER flyby's
    Eigen::Vector6i messengerFlyby1, messengerFlyby2, messengerFlyby3;
    messengerFlyby1 << 2008, 1, 14, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    messengerFlyby2 << 2008, 10, 6, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    messengerFlyby3 << 2009, 9, 29, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss

    // ABM integrator settings
    const double initialTimeStep = 3600;
    const double minimumStepSize = 3600/4;
    const double maximumStepSize = 3600*4;
    const double relativeErrorTolerence = 10E-12;
    const double absoluteErrorTolerence = 10E-12;
    const unsigned int minimumOrder = 6;
    const unsigned int maximumOrder = 12;

    // Parameter estimation settings
    const unsigned int maximumNumberOfIterations = 10;

    // Output location
    std::string outputSubFolder = tudat_applications::getOutputPath( ) + "Output/";


    /////////////////////
    //// ENVIRONMENT ////
    /////////////////////

    std::cout << "building environment..." << std::endl;

    // initial and final time to Julian
    double initialSimulationTime = secondsAfterJ2000(initialTime);
    double finalSimulationTime = secondsAfterJ2000(finalTime);

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

    // Default body settings
    double buffer = maximumOrder*maximumStepSize; //see Tudat libraries 1.1.3.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings = getDefaultBodySettings( bodyNames, initialSimulationTime - buffer, finalSimulationTime + buffer );

    // Use direct SPICE ephemerides
    if (useDirectSpice == true){
        std::string frameOrigin = "SSB";
        std::string frameOrientation = "ECLIPJ2000";
        for( unsigned int i = 0; i < totalNumberOfBodies; i++ )
        {
            bodySettings[bodyNames[i]]->ephemerisSettings = std::make_shared< DirectSpiceEphemerisSettings >( frameOrigin, frameOrientation );
        }
    }

    // Custom settings Sun
    double sunNormalizedJ2 = sunJ2 / calculateLegendreGeodesyNormalizationFactor(2,0);
    double sunNormalizedJ4 = sunJ4 / calculateLegendreGeodesyNormalizationFactor(4,0);

    Eigen::MatrixXd normalizedSineCoefficients;
    Eigen::MatrixXd normalizedCosineCoefficients;
    unsigned int maximumDegree;

    if (estimateJ4 == false){
        normalizedSineCoefficients   = Eigen::MatrixXd::Zero(3,3);
        normalizedCosineCoefficients = Eigen::MatrixXd::Zero(3,3);
        maximumDegree = 2;
    } else{
        normalizedSineCoefficients   = Eigen::MatrixXd::Zero(5,5);
        normalizedCosineCoefficients = Eigen::MatrixXd::Zero(5,5);
        normalizedCosineCoefficients(4,0) = sunNormalizedJ4;
        maximumDegree = 4;
    }
    normalizedCosineCoefficients(0,0) = 1.0;
    normalizedCosineCoefficients(2,0) = sunNormalizedJ2;

    bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                sunGravitationalParameter, sunRadius,
                normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );

    // Create body map
    NamedBodyMap bodyMap = createBodies( bodySettings );
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


    // Prepare angular momentum vector Sun
    const Eigen::Vector3d sunAngularMomentumVector(0.0, 0.0, sunAngularMomentum);

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
                                std::make_shared< SphericalHarmonicAccelerationSettings > (maximumDegree,0));
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                    calculateSchwarzschildCorrection,
                                    calculateLenseThirringCorrection,
                                    calculateDeSitterCorrection,
                                    "",
                                    sunAngularMomentumVector));


                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< TimeVaryingGravitationalParameterAccelerationSettings >(
                                    timeVaryingGravitationalParameter));

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




    //////////////////////////////
    //// PROPAGATION SETTINGS ////
    //////////////////////////////

    std::cout << "defining propagation settings..." << std::endl;

    // Get initial state of bodies to be propagated
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, initialSimulationTime );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate,
              systemInitialState, finalSimulationTime);

    // Define numerical integrator settings.
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

    std::cout << "running dynamics simulator...";

    SingleArcDynamicsSimulator <> dynamicsSimulator (bodyMap,integratorSettings,propagatorSettings,true,false,true);

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Retrieve numerically integrated state for each body.
    std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
    allBodiesPropagationHistory.resize( bodiesToPropagate.size() );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        for( unsigned int i = 0; i < bodiesToPropagate.size(); i++ )
        {
            allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
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
    }

    std::cout << " done, state history saved" << std::endl;





    ////////////////////////////////
    //// CREATE GROUND STATIONS ////
    ////////////////////////////////

    std::cout << "creating ground stations..." << std::endl;


    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;

    groundStationNames.push_back( "Station1" );
    createGroundStation( bodyMap.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );

//    groundStationNames.push_back( "Station2" );
//    createGroundStation( bodyMap.at( "Earth" ), "Station2",
//                         ( Eigen::Vector3d( ) << 0.0, -1.55, 2.0 ).finished( ), geodetic_position );

//    groundStationNames.push_back( "Station3" );
//    createGroundStation( bodyMap.at( "Earth" ), "Station3",
//                         ( Eigen::Vector3d( ) << 0.0, 0.8, 4.0 ).finished( ), geodetic_position );


//    // At centre of Mercury
//    createGroundStation( bodyMap.at( "Mercury" ), "MercuryCenter",
//                         ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( ), geodetic_position );




    ///////////////////////////////////////////
    //// DEFINE LINK ENDS FOR OBSERVATIONS ////
    ///////////////////////////////////////////

    std::cout << "defining link ends for observations..." << std::endl;


    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
//    std::vector< LinkEnds > stationReceiverLinkEnds;
//    std::vector< LinkEnds > stationTransmitterLinkEnds;

//    LinkEnds linkEnds;
//    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
//    {

//        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
//        linkEnds[ receiver ] = std::make_pair( "Mercury", "" );
//        stationTransmitterLinkEnds.push_back( linkEnds );
//        linkEnds.clear( );

//        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
//        linkEnds[ transmitter ] = std::make_pair( "Mercury", "" );
//        stationReceiverLinkEnds.push_back( linkEnds );
//        linkEnds.clear( );
//    }


    // Position Observable link ends
    std::vector< LinkEnds> positionObservableLinkEnds;
    LinkEnds linkEnds2;
    linkEnds2[ observed_body ] = std::make_pair( "Mercury", "" );
    positionObservableLinkEnds.push_back( linkEnds2 );


    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;

//    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );

//    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );

//    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
//    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    linkEndsPerObservable[ position_observable ].push_back( positionObservableLinkEnds[ 0 ] );




    /////////////////////////////////////
    //// CREATE OBSERVATION SETTINGS ////
    /////////////////////////////////////

    std::cout << "creating observation settings..." << std::endl;

    observation_models::ObservationSettingsMap observationSettingsMap;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define settings for observable, no light-time corrections, and biases for selected 1-way range links
            observationSettingsMap.insert( std::make_pair( currentLinkEndsList.at( i ),
                                                           std::make_shared< ObservationSettings >( currentObservable ) ) );
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


    // Add additional parameters to be estimated
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_gamma ) );
    varianceVector.push_back(sigmaGamma*sigmaGamma);

    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_beta ) );
    varianceVector.push_back(sigmaBeta*sigmaBeta);

    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("Sun", gravitational_parameter));
    varianceVector.push_back(sigmaSunGM*sigmaSunGM);

    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("Sun", time_varying_gravitational_parameter));
    varianceVector.push_back(sigmaTVGP*sigmaTVGP);

    // Add gravitational moments Sun
    std::vector< std::pair< int, int > > blockIndices;
    blockIndices.push_back(std::make_pair(2,0));
    varianceVector.push_back(sigmaSunJ2*sigmaSunJ2);
    if (estimateJ4 == true){
        blockIndices.push_back(std::make_pair(4,0));
        varianceVector.push_back(sigmaSunJ4*sigmaSunJ4);
    }
    parameterNames.push_back(std::make_shared<SphericalHarmonicEstimatableParameterSettings>
                             (blockIndices,"Sun",spherical_harmonics_cosine_coefficient_block));


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
    std::vector < double > messengerFlybyList;
    messengerFlybyList.push_back(secondsAfterJ2000(messengerFlyby1));
    messengerFlybyList.push_back(secondsAfterJ2000(messengerFlyby2));
    messengerFlybyList.push_back(secondsAfterJ2000(messengerFlyby3));

    std::vector< double > baseTimeList =
            makeObservationTimeList(secondsAfterJ2000(observationInitialTime),
                                    secondsAfterJ2000(observationFinalTime),
                                    observationTimeStep,
                                    messengerFlybyList);


    std::cout << "amount of observation times (before implementing viability settings): " << baseTimeList.size() << std::endl;

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

//    // Create observation viability settings and calculators
//    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
//    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
//                                                body_avoidance_angle, std::make_pair( "Mercury", "" ), "Sun",
//                                                unit_conversions::convertDegreesToRadians( 35.0 ) ) );
//    PerObservableObservationViabilityCalculatorList viabilityCalculators = createObservationViabilityCalculators(
//                bodyMap, linkEndsPerObservable, observationViabilitySettings );


    // Create noise functions per observable
    if (simulateObservationNoise == false){
        rangeNoise = 0.0;
        angularPositionNoise = 0.0;
        dopplerNoise = 0.0;
    }

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

    noiseFunctions[ position_observable ] =
            std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                       createBoostContinuousRandomVariableGeneratorFunction(
                           normal_boost_distribution, { 0.0, positionObservableNoise }, 0.0 ), std::placeholders::_1 );


    // MESSENGER noise based on SPE
//    noiseFunction[ one_way_range ] =
//            std::bind(  )


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
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ),
                noiseFunctions);




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
        // perturb body positions by random value between -10 and 10 meters
        parameterPerturbation.segment(i*6,3) = Eigen::Vector3d( (rand()%200-100.0)/10.0,
                                                                (rand()%200-100.0)/10.0,
                                                                (rand()%200-100.0)/10.0 );
        // perturb body velocities by random value between -0.01 and 0.01 m/s
        parameterPerturbation.segment(i*6,3) = Eigen::Vector3d( (rand()%20-10.0)/100.0,
                                                                (rand()%20-10.0)/100.0,
                                                                (rand()%20-10.0)/100.0 );
    }

    // perturb parameters with a value between -10ppm and +10ppm
    for( unsigned int i = 0; i < truthParameters.size()-6*numberOfNumericalBodies; i++ ){
        unsigned int j = i + 6*numberOfNumericalBodies;
        double randomPercentage = (rand()%20-10.0)/(1.0E6);
        std::cout << randomPercentage << std::endl;
        parameterPerturbation( j ) = initialParameterEstimate[j]*randomPercentage;
    }

    initialParameterEstimate += parameterPerturbation;


    std::cout << "True parameter values:" << std::endl;
    std::cout << truthParameters << std::endl;

    std::cout << "Initial guesses:" << std::endl;
    std::cout << initialParameterEstimate << std::endl;


    // Define a priori covariance matrix
    Eigen::MatrixXd aprioriCovariance =
        Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));

    Eigen::MatrixXd inverseOfAprioriCovariance;

    if (useAprioriValues == true){
        for( unsigned int i = 0; i < truthParameters.size(); i++ ){
            aprioriCovariance( i,i ) = varianceVector.at( i );
        }
        std::cout << "a priori covariance matrix:" << std::endl << aprioriCovariance << std::endl;
        Eigen::MatrixXd inverseOfAprioriCovariance = aprioriCovariance.inverse();
    } else{
        Eigen::MatrixXd inverseOfAprioriCovariance = Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ));
    }


    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                inverseOfAprioriCovariance,
                initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, true, false, true );

    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 );
    weightPerObservable[ position_observable ] = 1.0 / ( 1.0 * 1.0 );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( maximumNumberOfIterations ) );



    /////////////////////////////////////////////
    //// PROVIDE OUTPUT TO CONSOLE AND FILES ////
    /////////////////////////////////////////////

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd parameterOutcome = podOutput->parameterEstimate_;
    std::cout << "Outcome estimation:" << std::endl << ( parameterOutcome ).transpose( ) << std::endl;

//    std::cout << podOutput->getParameterHistoryMatrix( ) << std::endl;

    Eigen::VectorXd estimationError = parameterOutcome - truthParameters;
    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;

    std::cout << "Writing output to files..." << std::endl;

    std::map< double, Eigen::VectorXd > propagatedErrors;
    propagateFormalErrors(
                propagatedErrors, podOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                60.0, initialSimulationTime + 3600.0, finalSimulationTime - 3600.0 );
    input_output::writeDataMapToTextFile( propagatedErrors,
                                          "EstimationPropagatedErrors.dat",
                                          outputSubFolder );
    input_output::writeMatrixToFile( truthParameters,
                                     "TruthParameters.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getUnnormalizedCovarianceMatrix( ),
                                     "EstimationCovariance.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->normalizedInformationMatrix_,
                                     "EstimationInformationMatrix.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->weightsMatrixDiagonal_,
                                     "EstimationWeightsDiagonal.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->residuals_,
                                     "EstimationResiduals.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getCorrelationMatrix( ),
                                     "EstimationCorrelations.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getResidualHistoryMatrix( ),
                                     "ResidualHistory.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getParameterHistoryMatrix( ),
                                     "ParameterHistory.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( getConcatenatedMeasurementVector( podInput->getObservationsAndTimes( ) ),
                                     "ObservationMeasurements.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedTimeVector( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationTimes.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedGroundStationIndex( podInput->getObservationsAndTimes( ) ).first ),
                                     "ObservationLinkEnds.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector(
                                         getConcatenatedObservableTypes( podInput->getObservationsAndTimes( ) ) ),
                                     "ObservationObservableTypes.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( estimationError,
                                     "ObservationTrueEstimationError.dat", 16, outputSubFolder );
    input_output::writeMatrixToFile( podOutput->getFormalErrorVector( ),
                                     "ObservationFormalEstimationError.dat", 16, outputSubFolder );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    std::cout << "done!" << std::endl;
    return EXIT_SUCCESS;

}
