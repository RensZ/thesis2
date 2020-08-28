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

double noiseLevelBasedOnMSEangleForMultipleMissions(const double time,
                            std::vector< double > noiseAtMinAngleVector,
                            std::vector< double > noiseAtMaxAngleVector,
                            std::vector< double > maxAngleDegVector, //in degrees
                            std::vector< std::vector< double > > seperateBaseTimeLists);

double noiseSampleBasedOnMSEangleForMultipleMissions(const double time,
                            std::vector< double > noiseAtMinAngleVector,
                            std::vector< double > noiseAtMaxAngleVector,
                            std::vector< double > maxAngleDegVector,
                            std::vector< std::vector< double > > seperateBaseTimeLists);


double averageOfDoubleVector(std::vector<double> input);

Eigen::MatrixXd interpolatePositionErrors(const Eigen::MatrixXd errorMatrix,
                                          const std::vector<double> timesAtWhichToInterpolate);

std::map< double, Eigen::Vector3d > interpolatePositionErrorsBasedOnTrueAnomaly(
        const Eigen::MatrixXd errorMatrix,
        const std::vector<double> timesAtWhichToInterpolate,
        const std::string vehicle,
        const double mercuryGravitationalParameter);

Eigen::Matrix3d transformAngularMomentumFromLocalToGlobalFrame(
        const Eigen::Matrix3d angularMomentumInLocalFrame,
        const std::string localFrame,
        const std::string globalFrame,
        const double currentTime);

double phaseAccordingToSolarMinimum(
        const double solarMinimumEpoch,
        const double period);

double simpleSine(
        const double amplitude,
        const double period,
        const double phase,
        const double time);

std::map< double, Eigen::MatrixXd > tabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime,
        std::function< double() > amplitudeFunction,
        std::function< double() > periodFunction,
        std::function< double() > phaseFunction);

std::map< double, Eigen::MatrixXd > zeroTabulatedSphericalHarmonicsCoefficientCorrections(
        const double initialTime,
        const double finalTime);

Eigen::MatrixXd calculateConsiderCovarianceMatrix(
        const Eigen::MatrixXd P,
        const Eigen::VectorXd W_diagonal,
        const Eigen::MatrixXd C,
        const Eigen::MatrixXd Hx,
        const Eigen::MatrixXd Hc,
        std::string outputSubFolder);

Eigen::MatrixXd calculateConsiderCovarianceOfAsteroids(
        const Eigen::MatrixXd P,
        const Eigen::VectorXd W_diagonal,
        const Eigen::MatrixXd Hx,
        std::string inputFolderAsteroids,
        std::string inputFolderPartials,
        std::string outputSubFolder);

std::map< double, Eigen::MatrixXd > onlyEveryXthValueFromDataMap(
        std::map< double, Eigen::MatrixXd > inputMap,
        const int reductionFactor);

std::map< double, Eigen::Vector6d > onlyEveryXthValueFromDataMap(
        std::map< double, Eigen::Vector6d > inputMap,
        const int reductionFactor);

std::map< double, Eigen::VectorXd > onlyEveryXthValueFromDataMap(
        std::map< double, Eigen::VectorXd > inputMap,
        const int reductionFactor);

std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > onlyEveryXthValueFromDataMap(
        std::map< double, Eigen::Matrix<long double, Eigen::Dynamic, 1> > inputMap,
        const int reductionFactor);

std::string printScenario(const int scenario);

double combinedRangeAndSatelliteErrorLevel( const double observationTime,
                                            const Eigen::Vector3d satellitePositionErrorLevel,
                                            const double rangeNoiseLevel);

double satelliteErrorLevel( const double observationTime,
                            const Eigen::Vector3d satellitePositionErrorLevel);

std::map< unsigned int, std::pair< double, double > > readAsteroidsFile (
        std::string filename, std::string delimiter);

#endif // CUSTOMFUNCTIONS_H
