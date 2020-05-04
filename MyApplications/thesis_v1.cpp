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
    using namespace tudat::physical_constants;


    ////////////////////////////
    //// APPLICATION INPUTS ////
    ////////////////////////////

    // Acceleration settings
    const bool calculateLenseThirringCorrection = false; //need to implement partials
    const bool calculateDeSitterCorrection = false; //need to implement partials
    const bool estimateJ4 = false;

    // Parameter apriori values and uncertainties
    const bool useAprioriValues = true;

    // Planet propagation settings
    const bool propogatePlanets = false; // Propogate the other planets besides Mercury (NOTE: need observations for other planets, or LS can't find solutions for other planets)
    const bool useDirectSpice = false; // Use direct SPICE (more accurate than tabulated Spice)

    // Parameter estimation settings
    const unsigned int maximumNumberOfIterations = 5;

    // ABM integrator settings (if RK4 is used instead, initialstepsize is taken)
    const double initialTimeStep = 3600;
    const double minimumStepSize = 3600/8;
    const double maximumStepSize = 3600*8;
    const double relativeErrorTolerence = 10E-12;
    const double absoluteErrorTolerence = 10E-12;
    const unsigned int minimumOrder = 6;
    const unsigned int maximumOrder = 12;

    // Observation settings
    const bool simulateObservationNoise = true;
    const bool includeSpacecraftPositionError = true;

    // Use constraint nordtvedt=4*beta-gamma-3
    const bool useNordtvedtConstraint = true;


    ////////////////////////
    //// MISSION INPUTS ////
    ////////////////////////


    // Load json input
    std::string input_filename = "inputs_Genova2018.json"; // Messenger simulation done in Genova et al 2018, Nature Communications
//    std::string input_filename = "inputs_Genova2018_test.json"; // 1-year version of the simulation above for quicker test
//    std::string input_filename = "inputs_Schettino2015.json"; // BepiColombo simulation done in Schettino et al 2015, IEEE

    using json = nlohmann::json;
    std::string json_directory = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/";
    std::ifstream json_file(json_directory + input_filename);
    json json_input;
    json_file >> json_input;

    std::cout<<"inputs imported from json file: "<<json_input<<std::endl;


    // Retrieve simulation start and end time
    std::vector<int> json1 = json_input["initialTime"];
    Eigen::Vector6i initialTime(json1.data()); json1.clear();
    std::vector<int> json2 = json_input["finalTime"];
    Eigen::Vector6i finalTime(json2.data()); json2.clear();

    std::string vehicle = json_input["vehicle"];

    // Acceleration settings
    const bool calculateSchwarzschildCorrection = json_input["calculateSchwarzschildCorrection"];
    const bool includeSEPViolationAcceleration = json_input["includeSEPViolationAcceleration"];
    const bool includeTVGPAcceleration = json_input["includeTVGPAcceleration"];

    // Retrieve input parameters including uncertainties and apriori values
    const double sunJ2 = json_input["sunJ2"];
    const double sunAngularMomentum = json_input["sunAngularMomentum"];
    const double sunGravitationalParameter = json_input["sunGravitationalParameter"];
    const double nordtvedtParameter = json_input["nordtvedtParameter"];
    const double timeVaryingGravitationalParameter = json_input["timeVaryingGravitationalParameter"];

    const double sigmaGamma = json_input["sigmaGamma"];
    const double sigmaBeta = json_input["sigmaBeta"];
    const double sigmaNordtvedt = json_input["sigmaNordtvedt"];
    const double sigmaSunGM = json_input["sigmaSunGM"];
    const double sigmaSunJ2 = json_input["sigmaSunJ2"];
    const double sigmaTVGP = json_input["sigmaTVGP"];

    //    const double aprioriGamma = 1.0;
    //    const double aprioriBeta = 1.0;
    //    const double aprioriSunGM = sunGravitationalParameter;
    //    const double aprioriSunJ2 = sunJ2;
    //    const double aprioriSunJ4 = sunJ4;
    //    Eigen::VectorXd aprioriParameters(aprioriGamma, aprioriBeta, aprioriSunGM, aprioriSunJ2);


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
    double rangeNoise = json_input["rangeNoise"];
    double dopplerNoise = json_input["dopplerNoise"];
    double angularPositionNoise = json_input["angularPositionNoise"];

    const double noiseAtMinAngle = json_input["noiseAtMinAngle"];
    const double noiseAtMaxAngle = json_input["noiseAtMaxAngle"];

    // Location of simulation output
    std::string outputSubFolderName = json_input["outputSubFolderName"];
    std::string outputSubFolder = getOutputPath( ) + outputSubFolderName;

    // Other parameters, currently not included in json
    const double sunRadius = 695.7E6; //m, from nasa fact sheet
    const double sunJ4 = -4.34E-9; //from Antia 2008
    const double sigmaSunJ4 = 0.1E-9; //PLACEHOLDER

    const double sigmaPosition = 1000.0; //educated guess
    const double sigmaVelocity = 1.0; //educated guess

    double positionObservableNoise = 2.0;



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

    // Convert initial and final time to Julian
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
//    bodyNames[ 7 ] = "Uranus";
//    bodyNames[ 8 ] = "Neptune";
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

    // Process Nordtvedt constraint


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
                "Mercury", "Sun", maximumDegree, 0 ) );

    if (calculateSchwarzschildCorrection == true){
    dependentVariablesList.push_back(
              std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    relativistic_correction_acceleration, "Mercury" , "Sun" ) );
    }

    if (includeSEPViolationAcceleration == true){
    dependentVariablesList.push_back(
              std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    sep_violation_acceleration, "Mercury" , "Sun" ) );
    }

    if (includeTVGPAcceleration == true){
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

    // Get initial state of bodies to be propagated
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, initialSimulationTime );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate,
              systemInitialState, finalSimulationTime, cowell, dependentVariablesToSave);



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

    SingleArcDynamicsSimulator <> dynamicsSimulator (bodyMap,integratorSettings,propagatorSettings,true,false,true);

    std::cout << "saving integration result and dependent variables..." << std::endl;

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariablesHistory = dynamicsSimulator.getDependentVariableHistory();

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
    std::vector< double > baseTimeList =
            makeObservationTimeList(secondsAfterJ2000(observationInitialTime),
                                    secondsAfterJ2000(observationFinalTime),
                                    observationTimeStep,
                                    trackingArcDuration,
                                    maximumNumberOfTrackingDays,
                                    unit_conversions::convertDegreesToRadians(180.0-35.0),
                                    flybyList);

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
    if (simulateObservationNoise == false){
        rangeNoise = 0.0;
        angularPositionNoise = 0.0;
        dopplerNoise = 0.0;
    }

    std::map< ObservableType, std::function< double( const double ) > > noiseFunctions;



    std::function< double( const double ) > mercuryOrbiterNoiseFunction;
    mercuryOrbiterNoiseFunction = [noiseAtMinAngle, noiseAtMaxAngle](const double time){
        return noiseSampleBasedOnMSEangle(time, noiseAtMinAngle, noiseAtMaxAngle);
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
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ),
                noiseFunctions);

    // similar container, but the "observation" will be the noise value instead of the actual observation
    PodInputDataType observationWeightsAndTimes = observationsAndTimes;

    // add spacecraft initial position error to the observations
//    PodInputDataType oldObservationsAndTimes = observationsAndTimes; //to verify
    if (includeSpacecraftPositionError == true){

        std::cout << "Adding satellite estimation initial position error" << std::endl;
//        Eigen::Vector3d constantSatelliteError; constantSatelliteError << 10.0, 10.0, 10.0;

        std::string vehicleErrorFilename = "/home/rens/tudatBundle/tudatApplications/thesis/MyApplications/Input/error_inputs_"+vehicle+".txt";
        Eigen::MatrixXd error_input = input_output::readMatrixFromFile(vehicleErrorFilename, ",");

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

                    Eigen::MatrixXd interpolatedErrorMatrix = interpolatePositionErrors(error_input, allObservationTimes);

                    ObservationVectorType newObservations = Eigen::VectorXd(allObservations.size());
                    ObservationVectorType observationWeights = Eigen::VectorXd(allObservations.size());

                    // for every observation, retrieve and add the range bias that should be added
                    for (unsigned int i=0; i<allObservationTimes.size(); i++){
                        double observationTime = allObservationTimes.at( i );

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

                        double noiseLevel = abs(rangeCorrection)
                                + noiseLevelBasedOnMSEangle(observationTime, noiseAtMinAngle, noiseAtMaxAngle);
                        observationWeights(i) = 1.0/(noiseLevel*noiseLevel);

                        std::cout<<observationTime<<" // "<<
                                   currentSatelliteError.transpose()<<" // "<<
                                   randomErrorSample.transpose()<<" // "<<
                                   rangeCorrection<<" // "<<
                                   noiseLevel<<std::endl;
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
    podInput->defineEstimationSettings( true, true, false, true );
//    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->manuallySetObservationWeightMatrix(observationWeightsAndTimes);

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



    // Propagate errors to total time period
    std::cout<<"propagating errors..."<<std::endl;

    Eigen::MatrixXd fullStateHistory = input_output::readMatrixFromFile(
                outputSubFolder + "/StatePropagationHistoryMercury.dat");


    Eigen::MatrixXd initialCovarianceMatrix = podOutput->getUnnormalizedCovarianceMatrix( );
    std::map< double, Eigen::VectorXd > propagatedErrors;
    propagateFormalErrors(
                propagatedErrors, initialCovarianceMatrix,
                orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                //baseTimeList);
                3600.0, initialSimulationTime + 3600.0, finalSimulationTime - 3600.0 );

    if (fullStateHistory.col(0).size() == static_cast<long>(propagatedErrors.size())){
        std::cout<<"check!"<<std::endl;
    } else{
        std::cout<<fullStateHistory.col(0).size()<<propagatedErrors.size()<<std::endl;
    }

    std::map<double, Eigen::Vector6d> propagatedErrorsRSW;
    Eigen::Vector6d currentPropagatedError;
    std::map<double, Eigen::VectorXd>::iterator it = propagatedErrors.begin();
    unsigned int i = 0;

    while (it != propagatedErrors.end()){

        Eigen::Matrix3d currentTransformationToRSW =
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(fullStateHistory.block(i,1,1,3));

        currentPropagatedError.segment(0,3) = currentTransformationToRSW * it->second.segment(0,3);
        currentPropagatedError.segment(3,3) = currentTransformationToRSW * it->second.segment(3,3);

        propagatedErrorsRSW.insert(std::make_pair( it->first, currentPropagatedError ) );

        i++; it++;
    }

    // Propagate covariance matrix
    std::cout<<"propagating covariance matrix..."<<std::endl;

    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(propagatedCovariance, initialCovarianceMatrix,
                        orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( ),
                        //baseTimeList);
                        3600.0, initialSimulationTime + 3600.0, finalSimulationTime - 3600.0 );

    std::map<double, Eigen::VectorXd> propagatedErrorUsingCovMatrix;
    std::map<double, Eigen::VectorXd> propagatedRSWErrorUsingCovMatrix;

    unsigned int j = 0;
    std::map<double, Eigen::MatrixXd>::iterator covit = propagatedCovariance.begin();
    Eigen::Vector6d currentSigmaVectorRSW;
    while (covit != propagatedCovariance.end()){

        Eigen::VectorXd currentSigmaVector = covit->second.diagonal().cwiseSqrt();
        Eigen::Matrix3d currentTransformationToRSW =
                reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(fullStateHistory.block(j,1,1,3));

        currentSigmaVectorRSW.segment(0,3) = currentTransformationToRSW * currentSigmaVector.segment(0,3);
        currentSigmaVectorRSW.segment(3,3) = currentTransformationToRSW * currentSigmaVector.segment(3,3);

        propagatedErrorUsingCovMatrix.insert(std::make_pair( covit->first, currentSigmaVector ) );
        propagatedRSWErrorUsingCovMatrix.insert(std::make_pair( covit->first, currentSigmaVectorRSW ) );

        j++; covit++;
    }

    // Save data in files
    std::cout<<"writing output to files..."<<std::endl;



    input_output::writeDataMapToTextFile( propagatedErrors,
                                          "PropagatedErrors.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedErrorsRSW,
                                          "PropagatedErrorsRSW.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedCovariance,
                                          "PropagatedCovariance.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedErrorUsingCovMatrix,
                                          "propagatedErrorUsingCovMatrix.dat",
                                          outputSubFolder );

    input_output::writeDataMapToTextFile( propagatedRSWErrorUsingCovMatrix,
                                          "propagatedRSWErrorUsingCovMatrix.dat",
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


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    std::cout << "done!" << std::endl;
    return EXIT_SUCCESS;

}
