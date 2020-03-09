#ifndef TIMEVARYINGGRAVITATIONALPARAMETER_H
#define TIMEVARYINGGRAVITATIONALPARAMETER_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include <tudatApplications/thesis/MyApplications/TVGPInterface.h>

namespace tudat
{

namespace estimatable_parameters{


class TimeVaryingGravitationalParameter: public EstimatableParameter< double >
{

public:

    //! Constructor
    TimeVaryingGravitationalParameter(
            const std::shared_ptr < TVGPInterface > tvgpInterface,
            std::string& associatedBody ):
        EstimatableParameter< double >( time_varying_gravitational_parameter, associatedBody),
        tvgpInterface_( tvgpInterface )
    { }

    //! Destructor
    ~TimeVaryingGravitationalParameter(){}

    //! Function to get the current value.

    double getParameterValue( )
    {
        return tvgpInterface_->getTimeVaryingGravitationalParameter( );
    }

    //! Function to reset the value of the PPN parameter gamma.
    void setParameterValue( double parameterValue )
    {
        tvgpInterface_->setTimeVaryingGravitationalParameter( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }


protected:

private:

//! Object containing the radiation pressure coefficient to be estimated.
std::shared_ptr< TVGPInterface > tvgpInterface_;


};


};


}



#endif // TIMEVARYINGGRAVITATIONALPARAMETER_H
