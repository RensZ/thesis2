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

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <chrono>

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


    /////////////////////
    //// USER INPUTS ////
    /////////////////////

    // Specify start and end time
    Eigen::Vector6i initialTime, finalTime;
    initialTime << 2020, 1, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss
    finalTime   << 2020, 2, 1, 0, 0, 0; // YYYY, MM, DD, hh, mm, ss

    // Planet propagation choices
    bool propogatePlanets = false; // Propogate the other planets besides Mercury
    bool useDirectSpice = true; // Use direct SPICE (more accurate than tbulated Spice)

    // Sun inputs
    double sunRadius = 695.7E6; //m, from nasa fact sheet
    double sunJ2 = 2.0E-7; //from Hamid
    double sunGravitationalParameter = 132712E6; //km3/s2, from nasa fact sheet

    // ABM integrator settings
    double initialTimeStep = 3600;
    double minimumStepSize = 3600/8;
    double maximumStepSize = 3600*8;
    double relativeErrorTolerence = 10E-12;
    double absoluteErrorTolerence = 10E-12;
    int minimumOrder = 6;
    int maximumOrder = 12;




    /////////////////////
    //// ENVIRONMENT ////
    /////////////////////

    std::cout << "building environment..." << std::endl;

    // initial and final time to Julian
    double initialEphemerisTime = secondsAfterJ2000(initialTime);
    double finalEphemerisTime = secondsAfterJ2000(finalTime);

    // Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
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

    // Default body settings
    double buffer = maximumOrder*maximumStepSize; //see Tudat libraries 1.1.3.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

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

    // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
    SelectedAccelerationMap accelerationMap;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerations;
        for( unsigned int j = 0; j < bodiesToPropagate.size( ); j++ )
        {
            // Create central gravity acceleration between each 2 bodies.
            if( i != j )
            {
                currentAccelerations[ bodiesToPropagate.at( j ) ].push_back(
                            std::make_shared< AccelerationSettings >( central_gravity ) );
            }
        }
        accelerationMap[ bodiesToPropagate.at( i ) ] = currentAccelerations;
    }

    // Define central bodies to use in propagation.
    std::vector< std::string > centralBodies;
    centralBodies.resize( numberOfNumericalBodies );

    // Set central body as Solar System Barycenter for Sun and each planet, for Moon use the Earth
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        if (bodiesToPropagate[i] == "moon"){
            centralBodies[i] = "Earth";
        } else{
            centralBodies[i] = "SSB";
        }
    }

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );




    ////////////////////////////////
    //// ESTIMATABLE PARAMETERS ////
    ////////////////////////////////

    std::vector< boost::shared_ptr < EstimatableParameterSettings > > parameterNames;

    parameterNames.push_back(boost::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_gamma ) );
    parameterNames.push_back(boost::make_shared<EstimatableParameterSettings >
                             ("global_metric", ppn_parameter_beta ) );
    parameterNames.push_back(boost::make_shared<EstimatableParameterSettings >
                             ("Sun", gravitational_parameter));
    parameterNames.push_back(boost::make_shared<SphericalHarmonicEstimatableParameterSettings>
                             (2,0,2,0,"Sun",spherical_harmonics_cosine_coefficient_block));

    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
         createParametersToEstimate( parameterNames, bodyMap, );





    //////////////////////////////
    //// PROPAGATION SETTINGS ////
    //////////////////////////////

    std::cout << "defining propagation settings..." << std::endl;

    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodyMap, initialEphemerisTime );

//    // Define dependent variables to save
//    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
//    dependentVariablesList.push_back(
//        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
//        accelerationModelType, bodyUndergoingAcceleration, bodyExertingAcceleration, useNorm, componentIndex );
//    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
//        std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );
//    https://tudat.tudelft.nl/tutorials/tudatFeatures/propagationSetup/propagatorSettingsDependentVariables.html


    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate,
              systemInitialState, finalEphemerisTime);

    // Define numerical integrator settings.
    std::shared_ptr< AdamsBashforthMoultonSettings< double > > integratorSettings =
            std::make_shared< AdamsBashforthMoultonSettings< double > > (
                initialEphemerisTime, initialTimeStep,
                minimumStepSize, maximumStepSize,
                relativeErrorTolerence, absoluteErrorTolerence,
                minimumOrder, maximumOrder);




    //////////////////////////
    //// PROPAGATE ORBITS ////
    //////////////////////////

    std::cout << "propagating orbits..." << std::endl;

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings );

    std::map< double, Eigen::VectorXd > integrationResult =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    // Retrieve numerically integrated state for each body.
    std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
    allBodiesPropagationHistory.resize( numberOfNumericalBodies );
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
         stateIterator != integrationResult.end( ); stateIterator++ )
    {
        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
            allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
        }
    }




    ///////////////////////////////
    //// WRITE OUTPUT TO FILES ////
    ///////////////////////////////

    std::cout << "writing output files..." << std::endl;

    std::string outputSubFolder = "Output/";

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        // Write propagation history to file.
        input_output::writeDataMapToTextFile(
                    allBodiesPropagationHistory[ i ],
                    "innerSolarSystemPropagationHistory" + bodiesToPropagate.at( i ) + ".dat",
                    tudat_applications::getOutputPath( ) + outputSubFolder,
                    "",
                    std::numeric_limits< double >::digits10,
                    std::numeric_limits< double >::digits10,
                    "," );
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    std::cout << "done!" << std::endl;
    return EXIT_SUCCESS;

}
