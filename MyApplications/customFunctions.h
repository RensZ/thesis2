#ifndef CUSTOMFUNCTIONS_H
#define CUSTOMFUNCTIONS_H

#include <Tudat/SimulationSetup/tudatEstimationHeader.h>

double secondsAfterJ2000(Eigen::Vector6i datetime);
double secondsAfterJ2000(std::vector<int> datetime);

//double mercurySunEarthAngle(const double time);

double angleBetween2Bodies(const double time,
                           const std::string centralBody,
                           const std::string body1,
                           const std::string body2);

std::vector<double> makeObservationTimeList(const double initialTime,
                                            const double endTime,
                                            const double timeStep,
                                            const unsigned int maximumNumberOfObservations, //if not applicable, set to a very large number.
                                            const double maxMSEangle, //radians. if not applicable, set to value greater than pi.
                                            const std::vector<double> flybyTimes);

double noiseBasedOnMSEangle(const double time,
                            const double noiseAtMinAngle,
                            const double noiseAtMaxAngle);

double averageOfDoubleVector(std::vector<double> input);

#endif // CUSTOMFUNCTIONS_H
