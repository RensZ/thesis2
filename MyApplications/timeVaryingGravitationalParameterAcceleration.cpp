#include <Eigen/Geometry>
#include <tudatApplications/thesis/MyApplications/timeVaryingGravitationalParameterAcceleration.h>

namespace tudat{

Eigen::Vector3d calculateTimeVaryingGravitationalParameterAcceleration(
        const double centralBodyGravitationalParameter,
        const Eigen::Vector3d& relativePosition,
        const double timeVaryingGravitationalParameter,
        const double timeSinceJ2000 )
{
    // convert to arrays for easy element-wise multiplication
    Eigen::Array3d relativePositionArray = relativePosition;
    Eigen::Array3d relativePositionArrayCubed =
            relativePositionArray * relativePositionArray * relativePositionArray;

    // formula source: Equation 11 of Genova et al 2018, Nature communications
    Eigen::Vector3d acceleration =
            centralBodyGravitationalParameter *
            timeVaryingGravitationalParameter *
            timeSinceJ2000 *
            relativePositionArray /
            relativePositionArrayCubed;
    return acceleration;
}

}


//timeVaryingGravitationalParameterAcceleration::timeVaryingGravitationalParameterAcceleration()
//{

//}
