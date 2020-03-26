#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

static inline std::string getOutputPath(
        const std::string& extraDirectory = "" );

double secondsAfterJ2000(Eigen::Vector6i datetime);

double mercurySunEarthAngle(const double time);

std::vector<double> makeObservationTimeList(const double initialTime,
                                            const double endTime,
                                            const double timeStep,
                                            const unsigned int maximumNumberOfObservations, //if not applicable, set to a very large number.
                                            const double maxMSEangle, //radians. if not applicable, set to value greater than pi.
                                            const std::vector<double> flybyTimes);

double noiseBasedOnMSEangle(const double time);


#endif // CUSTOMFUNCTIONS_H
