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
                                            const double arcDuration,
                                            const unsigned int maximumNumberOfObservations, //if not applicable, set to a very large number.
                                            const double maxMSEangle, //radians. if not applicable, set to value greater than pi.
                                            const std::vector<double> flybyTimes);

double noiseLevelBasedOnMSEangle(const double time,
                            const double noiseAtMinAngle,
                            const double noiseAtMaxAngle,
                            const double maxAngleDeg);

double noiseSampleBasedOnMSEangle(const double time,
                            const double noiseAtMinAngle,
                            const double noiseAtMaxAngle,
                            const double maxAngleDeg);

double noiseSampleBasedOnMSEangleForMultipleMissions(const double time,
                            std::vector< double > noiseAtMinAngleVector,
                            std::vector< double > noiseAtMaxAngleVector,
                            std::vector< double > maxAngleDegVector,
                            std::vector< std::vector< double > > seperateBaseTimeLists);


double averageOfDoubleVector(std::vector<double> input);

Eigen::MatrixXd interpolatePositionErrors(const Eigen::MatrixXd errorMatrix,
                                          const std::vector<double> timesAtWhichToInterpolate);

Eigen::MatrixXd interpolatePositionErrorsBasedOnTrueAnomaly(const Eigen::MatrixXd errorMatrix,
                                                            const std::vector<double> timesAtWhichToInterpolate,
                                                            const std::string vehicle,
                                                            const double mercuryGravitationalParameter);

Eigen::Matrix3d transformAngularMomentumFromLocalToGlobalFrame(
        const Eigen::Matrix3d angularMomentumInLocalFrame,
        const std::string localFrame,
        const std::string globalFrame,
        const double currentTime);

double simpleSine(
        const double amplitude,
        const double period,
        const double phase,
        const double time);

std::map< double, Eigen::MatrixXd > tabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime,
        const double amplitude,
        const double period,
        const double phase);

std::map< double, Eigen::MatrixXd > zeroTabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime);

#endif // CUSTOMFUNCTIONS_H
