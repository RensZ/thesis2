/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TVGPINTERFACE_H
#define TVGPINTERFACE_H

#include <memory>

#include <Eigen/Core>

namespace tudat
{


class TVGPInterface
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     * \param gravitationalParameter Gravitational parameter associated with gravity field
     * \param updateInertiaTensor Function that is to be called to update the inertia tensor (typicaly in Body class; default none)
     */
    TVGPInterface( const double TimeVaryingGravitationalParameter):
        TimeVaryingGravitationalParameter_( TimeVaryingGravitationalParameter )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~TVGPInterface( ) { }

    //! Set the gravitational parameter.
    /*!
     * Define the gravitational parameter in meter^3 per second^2.
     * \param TimeVaryingGravitationalParameter New gravitational parameter associated with gravity field.
     */
    void resetTimeVaryingGravitationalParameter( const double TimeVaryingGravitationalParameter )
    {
        TimeVaryingGravitationalParameter_ = TimeVaryingGravitationalParameter;
    }

    //! Get the gravitational parameter.
    /*!
     * Return the gravitational parameter in meter^3 per second^2.
     * \return Gravitational parameter.
     */
    double getTimeVaryingGravitationalParameter( )
    {
        return TimeVaryingGravitationalParameter_;
    }


protected:

    //! Gravitational parameter.
    /*!
     * The gravitational parameter in meter^3 per second^2.
     */
    double TimeVaryingGravitationalParameter_;

private:
};

}

#endif // TVGPINTERFACE_H
