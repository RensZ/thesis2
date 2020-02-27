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
    int secondsPerDay = 60*60*24;
    int julianDayJ2000 = 2451545;
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay(datetime[0],datetime[1],datetime[2],datetime[3],datetime[4],datetime[5]);
    return (julianDay - julianDayJ2000)*secondsPerDay;
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

    // Specify start and end time
    Eigen::Vector6i initialTime, finalTime;
    initialTime << 2020, 1, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    finalTime   << 2020, 7, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss

    // Planet propagation choices
    bool propogatePlanets = false; // Propogate the other planets besides Mercury (NOTE: need observations as well of the other planets in order to do a viable state estimation)
    bool useDirectSpice = false; // Use direct SPICE (more accurate than tabulated Spice)

    // Observation settings
    bool simulateObservationNoise = false;

    // ABM integrator settings
    double initialTimeStep = 3600;
    double minimumStepSize = 3600/2;
    double maximumStepSize = 3600*2;
    double relativeErrorTolerence = 10E-12;
    double absoluteErrorTolerence = 10E-12;
    int minimumOrder = 6;
    int maximumOrder = 12;

    // Parameter estimation settings
    int maximumNumberOfIterations = 10;

    // Parameter inputs
    double sunRadius = 695.7E6; //m, from nasa fact sheet
    double sunJ2 = 2.0E-7; //from Hamid
    double sunGravitationalParameter = 132712E15; //m3/s2, from nasa fact sheet


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

    Eigen::Matrix3d normalizedCosineCoefficients;
    normalizedCosineCoefficients << 1.0,              0.0, 0.0,
                                    0.0,              0.0, 0.0,
                                    sunNormalizedJ2 , 0.0, 0.0 ;
    Eigen::Matrix3d normalizedSineCoefficients = Eigen::Matrix3d::Zero( );

    bodySettings[ "Sun" ] -> gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                sunGravitationalParameter, sunRadius,
                normalizedCosineCoefficients, normalizedSineCoefficients, "IAU_Sun" );

    // Create body map
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );



    ////////////////////////////////
    //// CREATE GROUND STATIONS ////
    ////////////////////////////////

    std::cout << "creating ground stations..." << std::endl;


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



//    // At centre of Earth
//    std::vector< std::string > EarthStationNames;
//    EarthStationNames.push_back( "Earth1" );
//    createGroundStation( bodyMap.at( "Earth" ), "Earth1",
//                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );

    // At centre of Mercury
    createGroundStation( bodyMap.at( "Mercury" ), "Mercury0",
                         ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( ), spherical_position );





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
                                std::make_shared< SphericalHarmonicAccelerationSettings > (2,0));
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




    ///////////////////////////////////////////
    //// DEFINE LINK ENDS FOR OBSERVATIONS ////
    ///////////////////////////////////////////

    std::cout << "defining link ends for observations..." << std::endl;


    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Mercury", "Mercury0" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Mercury", "Mercury0" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );



//    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
//    std::vector< LinkEnds > stationReceiverLinkEnds;
//    std::vector< LinkEnds > stationTransmitterLinkEnds;

//    for( unsigned int i = 0; i < EarthStationNames.size( ); i++ )
//    {
//        LinkEnds linkEnds;
//        linkEnds[ transmitter ] = std::make_pair( "Earth", EarthStationNames.at( i ) );
//        linkEnds[ receiver ] = std::make_pair( "Mercury", "Mercury0" );
//        stationTransmitterLinkEnds.push_back( linkEnds );

//        linkEnds.clear( );
//        linkEnds[ receiver ] = std::make_pair( "Earth", EarthStationNames.at( i ) );
//        linkEnds[ transmitter ] = std::make_pair( "Mercury", "Mercury0" );
//        stationReceiverLinkEnds.push_back( linkEnds );
//    }


//    // Define perfect Mercury observables
//    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
//    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );





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

    // Add bodies that will be propagated to the parameters to be estimated
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        int j = 6*i;

        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                      bodiesToPropagate[i], systemInitialState.segment(j,6), centralBodies[i] ) );
    }



    // Add additional parameters to be estimated
//    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
//                             ("global_metric", ppn_parameter_gamma ) );
//    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
//                             ("global_metric", ppn_parameter_beta ) );
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings >
                             ("Sun", gravitational_parameter));
    parameterNames.push_back(std::make_shared<SphericalHarmonicEstimatableParameterSettings>
                             (2,0,2,0,"Sun",spherical_harmonics_cosine_coefficient_block));

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap );

    // Print identifiers and indices of parameters to terminal
    printEstimatableParameterEntries( parametersToEstimate );




    ///////////////////////////////////////////////
    //// INITIALIZE ORBIT DETERMINATION OBJECT ////
    ///////////////////////////////////////////////

    std::cout << "Running OD manager..." << std::endl;

    // Create orbit determination object (propagate orbit, create observation models)
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodyMap, parametersToEstimate, observationSettingsMap,
                integratorSettings, propagatorSettings );




    ///////////////////////////////
    //// SIMULATE OBSERVATIONS ////
    ///////////////////////////////

    std::cout << "Simulating observations..." << std::endl;

    // Define start and end of observations
    double observationTimeStart = initialSimulationTime;
    double observationTimeEnd = finalSimulationTime;

    // Define time between two observations
    double observationInterval = 60.0*60.0*24.0;

    // Generate list of observation times
    std::vector< double > baseTimeList;
    double obsTime = observationTimeStart;
    while (obsTime <= observationTimeEnd){
        baseTimeList.push_back(obsTime);
        obsTime += observationInterval;
    }


    // Create measurement simulation input
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
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
                    std::make_pair( baseTimeList, receiver );
        }
    }

    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate observations
    PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ) );


    /////////////////////////////
    //// ESTIMATE PARAMETERS ////
    /////////////////////////////

    std::cout << "Estimating parameters..." << std::endl;


    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;


    if (simulateObservationNoise == false){

        // If no observation noise is simulated, we need to perturb initial parameter estimates to create an offset with the true solution
        Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

        for( unsigned int i = 0; i < truthParameters.size(); i++ ){
            double randomPercentage = ((rand()%200) - 100.0) / 10000.0; //random percentage between -1% and +1% to add to initial guess
            parameterPerturbation( i ) = initialParameterEstimate[i]*randomPercentage;
        }
        initialParameterEstimate += parameterPerturbation;
    }

    std::cout << "True parameter values:" << std::endl;
    std::cout << truthParameters << std::endl;

    std::cout << "Initial guesses:" << std::endl;
    std::cout << initialParameterEstimate << std::endl;

    // Define estimation input
    std::shared_ptr< PodInput< double, double > > podInput =
            std::make_shared< PodInput< double, double > >(
                observationsAndTimes, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );
    podInput->defineEstimationSettings( true, true, false, true );

    // Define observation weights (constant per observable type)
    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 );
    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    // Perform estimation
    std::shared_ptr< PodOutput< double > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( maximumNumberOfIterations ) );



    /////////////////////////////////////////////
    //// PROVIDE OUTPUT TO CONSOLE AND FILES ////
    /////////////////////////////////////////////

    std::cout << "Writing output to files..." << std::endl;

    std::string outputSubFolder = tudat_applications::getOutputPath( ) + "Output/";

    // Print true estimation error, limited mostly by numerical error
    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    std::cout << "True estimation error is:   " << std::endl << ( estimationError ).transpose( ) << std::endl;
    std::cout << "Formal estimation error is: " << std::endl << podOutput->getFormalErrorVector( ).transpose( ) << std::endl;

    std::map< double, Eigen::VectorXd > propagatedErrors;

    propagateFormalErrors(
                propagatedErrors, podOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                60.0, initialSimulationTime + 3600.0, finalSimulationTime - 3600.0 );
    input_output::writeDataMapToTextFile( propagatedErrors,
                                          "EstimationPropagatedErrors.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder );

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
