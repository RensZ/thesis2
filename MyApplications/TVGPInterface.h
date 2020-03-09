/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_METRIC_H
#define TUDAT_METRIC_H

#include <memory>

#include <Eigen/Core>

namespace tudat
{


//! Class that stores the PPN parameters, typically used as a 'global' environment property stored in ppnParameterSet variable
class TVGPInterface
{
public:

    //! Constructor
    TVGPInterface( const double timeVaryingGravitationalParameter):
        timeVaryingGravitationalParameter_( timeVaryingGravitationalParameter )
    { }

    //! Destructor
    ~TVGPInterface( ){ }

    //! Function to retrieve value
    double getTimeVaryingGravitationalParameter( )
    {
        return timeVaryingGravitationalParameter_;
    }


    //! Function to reset value
    void setTimeVaryingGravitationalParameter( const double timeVaryingGravitationalParameter )
    {
        timeVaryingGravitationalParameter_ = timeVaryingGravitationalParameter;
    }


protected:


    //! Value
    double timeVaryingGravitationalParameter_;

};

////! Global PPN parameter set, initialized upon compilation (with values equal to GR).
//extern std::shared_ptr< PPNParameterSet > ppnParameterSet;

////! Global parameter denoting EP violation in proper time rate, initialized to GR value of 0 upon compilation.
//extern double equivalencePrincipleLpiViolationParameter;

}

#endif // TVGPINTERFACE_H
