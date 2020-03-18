#ifndef TIMEVARYINGGRAVITATIONALPARAMETERACCELERATION_H
#define TIMEVARYINGGRAVITATIONALPARAMETERACCELERATION_H

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <memory>
#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{


//! Function to compute the additional acceleration term due to a dynamic GM of the central body.
//Eigen::Vector3d calculateTimeVaryingGravitationalParameterAcceleration(
//        const double centralBodyGravitationalParameter,
//        const Eigen::Vector3d relativePosition,
//        const double timeVaryingGravitationalParameter,
//        const double timeSinceJ2000 );

// moved everything to .h file


class TimeVaryingGravitationalParameterAcceleration:
        public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:


    //! Constructor
    TimeVaryingGravitationalParameterAcceleration(
            std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody,
            std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody,
            std::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            std::function< double( ) > timeVaryingGravitationalParameterFunction
            ):
        AccelerationModel< Eigen::Vector3d >( ),
        stateFunctionOfAcceleratedBody_( stateFunctionOfAcceleratedBody ),
        stateFunctionOfCentralBody_( stateFunctionOfCentralBody ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        timeVaryingGravitationalParameterFunction_( timeVaryingGravitationalParameterFunction )
    {
        this->updateMembers( );
    }


    //! Destructor
    ~TimeVaryingGravitationalParameterAcceleration( ){ }

    //! Function to return the current acceleration
    /*!
     * Returns the acceleration. Value is computed by updateMembers function
     * \return Acceleration.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }


    //! Update member variables used by the relativistic correction acceleration model.
    /*!
     * Updates member variables used by the relativistic correction acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentTime_ = currentTime;

            // Update common variables

            gravitationalParameterOfCentralBody_ =
                    gravitationalParameterFunctionOfCentralBody_( );

            stateOfAcceleratedBody_ = stateFunctionOfAcceleratedBody_( );
            stateOfCentralBody_ = stateFunctionOfCentralBody_( );

            stateOfAcceleratedBodyWrtCentralBody_ = stateOfAcceleratedBody_ - stateOfCentralBody_;

            positionOfAcceleratedBodyWrtCentralBody_ = stateOfAcceleratedBodyWrtCentralBody_.segment(0,3);

            double distance = positionOfAcceleratedBodyWrtCentralBody_.norm();

//            relativePositionArray_ = positionOfAcceleratedBodyWrtCentralBody_;


            // compute acceleration (equation 11 of Genova et al 2018, Nature communications)
            currentAcceleration_ =
                    gravitationalParameterOfCentralBody_ *
                    timeVaryingGravitationalParameter_ *
                    ( currentTime_/(tudat::physical_constants::JULIAN_YEAR) ) *
                    positionOfAcceleratedBodyWrtCentralBody_ /
                    (distance*distance*distance);


            timeVaryingGravitationalParameter_ = timeVaryingGravitationalParameterFunction_( );

//            currentAcceleration_ = calculateTimeVaryingGravitationalParameterAcceleration(
//                    gravitationalParameterOfCentralBody_,
//                    positionOfAcceleratedBodyWrtCentralBody_,
//                    timeVaryingGravitationalParameter_,
//                    currentTime_
//                    );

        }

     }

    //! Function to return the current state of the body undergoing acceleration
    /*!
     * Function to return the current state of the body undergoing acceleration
     * \return Current state of the body undergoing acceleration
     */
    std::function< Eigen::Vector6d( ) > getStateFunctionOfAcceleratedBody( )
    { return stateFunctionOfAcceleratedBody_; }

    //! Function to return the current state of the main body exerting acceleration
    /*!
     * Function to return the current state of the main body exerting acceleration
     * \return Current state of the main body exerting acceleration
     */
    std::function< Eigen::Vector6d( ) > getStateFunctionOfCentralBody( )
    { return stateFunctionOfCentralBody_; }

    //! Function to return the current gravitational parameter of central body
    /*!
     * Function to return the current gravitational parameter of central body
     * \return Current gravitational parameter of central body
     */
    std::function< double( ) > getGravitationalParameterFunctionOfCentralBody( )
    { return gravitationalParameterFunctionOfCentralBody_; }


    //! Function to return the current gravitational parameter of central body
    std::function< double( ) > getTimeVaryingGravitationalParameterFunction( )
    { return timeVaryingGravitationalParameterFunction_; }


private:

    // Functions

    //! State function of vehicle undergoing acceleration
    std::function< Eigen::Vector6d( ) > stateFunctionOfAcceleratedBody_;

    //! State function of main body exerting acceleration (e.g. Earth for an Earth-orbiting satellite).
    std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody_;

    //! Function returning the gravitational parameter of the central body
    std::function< double( ) > gravitationalParameterFunctionOfCentralBody_;

    //! Function returning the time varying gravitational parameter
    std::function< double( ) > timeVaryingGravitationalParameterFunction_;




    // Variables

    //! Current state of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfAcceleratedBody_;

    //! Current state of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfCentralBody_;

    //! Current state of the main body exerting acceleration, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfCentralBodyWrtPrimaryBody_;

    //! Current state of the primary body, as computed by last call to updateMembers function.
    Eigen::Vector6d stateOfAcceleratedBodyWrtCentralBody_;

    //! Current state of the primary body, as computed by last call to updateMembers function.
    Eigen::Vector3d positionOfAcceleratedBodyWrtCentralBody_;

    //! Current time varying gravitational parameter
    double gravitationalParameterOfCentralBody_;

    //! Current gravitational parameter of central body
    double timeVaryingGravitationalParameter_;

    //! acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d currentAcceleration_;

    // declare variables needed in the cpp file?
    Eigen::Array3d relativePositionArray_;


};

}




#endif // TIMEVARYINGGRAVITATIONALPARAMETERACCELERATION_H



