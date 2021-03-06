#ifndef SEPVIOLATIONACCELERATION_H
#define SEPVIOLATIONACCELERATION_H

#include <functional>
#include <iostream>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"


namespace tudat
{
namespace relativity
{

//! Class for SEP violation acceleration
class SEPViolationAcceleration:
        public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor
    SEPViolationAcceleration(
            std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratedBody,
            std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody,
            std::function< Eigen::Vector3d( ) > sepPositionCorrectionFunction,
            std::function< double( ) > gravitationalParameterFunctionOfCentralBody,
            std::function< Eigen::Vector3d( ) > nordtvedtPartialFunction,
            std::function< bool( ) > useNordtvedtConstraintFunction,
            std::function< double( ) > nordtvedtParameterFunction
            ):
        AccelerationModel< Eigen::Vector3d >( ),
        positionFunctionOfAcceleratedBody_( positionFunctionOfAcceleratedBody ),
        positionFunctionOfCentralBody_( positionFunctionOfCentralBody ),
        sepPositionCorrectionFunction_( sepPositionCorrectionFunction ),
        gravitationalParameterFunctionOfCentralBody_( gravitationalParameterFunctionOfCentralBody ),
        nordtvedtPartialFunction_(nordtvedtPartialFunction),
        useNordtvedtConstraintFunction_(useNordtvedtConstraintFunction),
        nordtvedtParameterFunction_( nordtvedtParameterFunction )
    {
        this->updateMembers( );
    }


    //! Destructor
    ~ SEPViolationAcceleration( ){ }

    //! Function to return the current acceleration
    Eigen::Vector3d getAcceleration( )
    {
        return sepRelativeAcceleration_;
    }

    //! Function to return the current acceleration
    Eigen::Vector3d getSEPPositionCorrection( )
    {
        return sepPositionCorrection_.cast<double>();
    }

    //! Update member variables used by the relativistic correction acceleration model.
    /*!
     * Updates member variables used by the relativistic correction acceleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * acceleration is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current position.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {

            currentTime_ = currentTime;

            // SEP corrections can be small enough for the calculatiosn to be sensitive to floating point precision
            // therefore everything is converted to long double and back for the calculations
            positionOfAcceleratedBody_ = positionFunctionOfAcceleratedBody_( ).cast<long double>();
            positionOfCentralBody_ = positionFunctionOfCentralBody_( ).cast<long double>();
            sepPositionCorrection_ = sepPositionCorrectionFunction_( ).cast<long double>();
            gravitationalParameterOfCentralBody_ =
                    static_cast< long double>(gravitationalParameterFunctionOfCentralBody_( ));

            sepCorrectedPositionOfCentralBody_ = positionOfCentralBody_ + sepPositionCorrection_;
            nordtvedtPartial_ = nordtvedtPartialFunction_( );

            centralGravityAcceleration_ = gravitation::computeGravitationalAccelerationLong(
                        positionOfAcceleratedBody_,
                        gravitationalParameterOfCentralBody_,
                        positionOfCentralBody_);

            sepCorrectedCentralGravityAcceleration_ = gravitation::computeGravitationalAccelerationLong(
                        positionOfAcceleratedBody_,
                        gravitationalParameterOfCentralBody_,
                        sepCorrectedPositionOfCentralBody_);

            nordtvedtParameter_ = nordtvedtParameterFunction_( );
            if (nordtvedtParameter_ == 0.0){
                sepRelativeAcceleration_ = Eigen::Vector3d::Zero();
            } else{
                sepRelativeAcceleration_ =
                        (sepCorrectedCentralGravityAcceleration_ - centralGravityAcceleration_).cast<double>();
            }

//            std::cout<<sepPositionCorrection_.transpose()<<" / ";
//            std::cout<<nordtvedtParameter_<<" / ";
//            std::cout<<sepRelativeAcceleration_.transpose()<<std::endl;
        }
    }

    //! Function to return the current position of the body undergoing acceleration
    /*!
     * Function to return the current position of the body undergoing acceleration
     * \return Current position of the body undergoing acceleration
     */
    std::function< Eigen::Vector3d( ) > getPositionFunctionOfAcceleratedBody( )
    { return positionFunctionOfAcceleratedBody_; }

    //! Function to return the current position of the main body exerting acceleration
    /*!
     * Function to return the current position of the main body exerting acceleration
     * \return Current position of the main body exerting acceleration
     */
    std::function< Eigen::Vector3d( ) > getPositionFunctionOfCentralBody( )
    { return positionFunctionOfCentralBody_; }


    std::function< Eigen::Vector3d( ) > getSEPPositionCorrectionFunction( )
    { return sepPositionCorrectionFunction_; }

    std::function< Eigen::Vector3d( ) > getNordtvedtPartialFunction( )
    { return nordtvedtPartialFunction_; }

    std::function< bool( ) > getUseNordtvedtConstraintFunction( )
    { return useNordtvedtConstraintFunction_; }




    //! Function to return the current gravitational parameter of central body
    std::function< double( ) > getGravitationalParameterFunctionOfCentralBody( )
    { return gravitationalParameterFunctionOfCentralBody_; }


    //! Function to return the time varying gravitational parameter
    std::function< double( ) > getNordtvedtParameterFunction( )
    { return nordtvedtParameterFunction_; }


private:

    // Functions

    //! Function returning the gravitational parameter of the central body
    std::function< std::string( ) > acceleratedBodyNameFunction_;

    //! Function returning the gravitational parameter of the central body
    std::function< std::string( ) > centralBodyNameFunction_;


    //! position function of vehicle undergoing acceleration
    std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratedBody_;

    //! position function of main body exerting acceleration (e.g. Earth for an Earth-orbiting satellite).
    std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody_;

    //! position function of main body exerting acceleration (e.g. Earth for an Earth-orbiting satellite).
    std::function< Eigen::Vector3d( ) > sepPositionCorrectionFunction_;



    //! Function returning the gravitational parameter of the accelerated body
    std::function< double( ) > gravitationalParameterFunctionOfAcceleratedBody_;

    //! Function returning the gravitational parameter of the central body
    std::function< double( ) > gravitationalParameterFunctionOfCentralBody_;

    //! position function of main body exerting acceleration (e.g. Earth for an Earth-orbiting satellite).
    std::function< Eigen::Vector3d( ) > nordtvedtPartialFunction_;

    std::function< bool( ) > useNordtvedtConstraintFunction_;


    //! Function returning the time varying gravitational parameter
    std::function< double( ) > nordtvedtParameterFunction_;



    // Variables

    //! Current position of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Matrix<long double, 3, 1> positionOfAcceleratedBody_;

    //! Current position of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Matrix<long double, 3, 1> positionOfCentralBody_;


    //! Current position of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Matrix<long double, 3, 1> sepPositionCorrection_;

    //! Current position of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Matrix<long double, 3, 1> sepCorrectedPositionOfCentralBody_;

    //! Current position of the body undergoing acceleration, as computed by last call to updateMembers function.
    Eigen::Vector3d nordtvedtPartial_;



    //! Current time varying gravitational parameter
    long double gravitationalParameterOfCentralBody_;

    //! Current gravitational parameter of central body
    double nordtvedtParameter_;


    //! acceleration, as computed by last call to updateMembers function
    Eigen::Matrix<long double, 3, 1> centralGravityAcceleration_;

    //! acceleration, as computed by last call to updateMembers function
    Eigen::Matrix<long double, 3, 1> sepCorrectedCentralGravityAcceleration_;

    //! acceleration, as computed by last call to updateMembers function
    Eigen::Vector3d sepRelativeAcceleration_;



};


} // namespace relativity
} // namespace tudat

#endif // SEPVIOLATIONACCELERATION_H
