/**
 * @file random-generator.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "random-generator.h"

/**
 * @brief Random sampling of position, according to https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf 
 * 
 * @return true When the point is valid, i. e. is in limits (does not check nn) 
 */
template<>
bool RandomGenerator<Point2D<double>>::RandomPointInDistance(const Point2D<double>& center, Point2D<double>& point, const double distance) {
  // rotation in polar
  double phi{uniDistAngle(rndEng)};
  
  point.SetPosition(center[0] + cos(phi) * distance, center[1] + sin(phi) * distance);
  return isInLimits(point);
}

template<>
bool RandomGenerator<Point2DDubins<double>>::RandomPointInDistance(const Point2DDubins<double>& center, Point2DDubins<double>& point, const double distance) {
  // rotation in polar
  double phi{uniDistAngle(rndEng)};
  
  point.SetPosition(center[0] + cos(phi) * distance, center[1] + sin(phi) * distance);
  point.SetAngle(uniDistAngle(rndEng));
  return isInLimits(point);
}

template<>
bool RandomGenerator<Point3D<double>>::RandomPointInDistance(const Point3D<double>& center, Point3D<double>& point, const double distance) {
  Point3D<double> temp;

  double phi{uniDistAngle(rndEng)};
  double theta{uniDistAngle(rndEng)};
  temp.SetPosition(center[0] + cos(theta) * sin(phi) * distance, center[1] + sin(theta) * sin(phi)  * distance, center[2] + cos(phi) * distance);

  // rotation
  double s{RandomProbability()};
  double sigOne{sqrt(1-s)};
  double sigTwo{sqrt(s)};
  double thetaOne{2*M_PI*RandomProbability()};
  double thetaTwo{2*M_PI*RandomProbability()};

  Quaternion<double> rotation{cos(thetaTwo) * sigTwo, sin(thetaOne) * sigOne, cos(thetaOne) * sigOne, sin(thetaTwo) * sigTwo};
  temp.SetRotation(rotation);

  // the distance is not exactly the expected one (might be much greater), therefore a point in the same direction with the correct distance is needed
  point = center.GetStateInDistance(temp, distance);

  return isInLimits(point);
}

/**
 * @brief Random space uniformly sampled in the configuration space of particular problem, according to 
 * https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf 
 * 
 * @param point Output = uniformly sampled point in the configuration space
 */
void RandomGenerator<Point2D<double>>::RandomPointInSpace(Point2D<double>& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng));
}

void RandomGenerator<Point2DDubins<double>>::RandomPointInSpace(Point2DDubins<double>& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng));
  point.SetAngle(uniDistAngle(rndEng));
}

void RandomGenerator<Point3D<double>>::RandomPointInSpace(Point3D<double>& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  // rotation
  double s{RandomProbability()};
  double sigOne{sqrt(1-s)};
  double sigTwo{sqrt(s)};
  double thetaOne{2*M_PI*RandomProbability()};
  double thetaTwo{2*M_PI*RandomProbability()};

  Quaternion<double> rotation{cos(thetaTwo) * sigTwo, sin(thetaOne) * sigOne, cos(thetaOne) * sigOne, sin(thetaTwo) * sigTwo};
  point.SetRotation(rotation);
}

/**
 * @brief Private function to check, whether the point lies in the workspace/environment
 * 
 * @param p Point to check
 * @return true Point DOES lie in the workspace
 * @return false Point DOES NOT lie in the workspace
 */
bool RandomGenerator<Point2D<double>>::isInLimits(Point2D<double>& p) {
  bool valid{true};
  valid &= p[0] >= limits.minX;
  valid &= p[0] <= limits.maxX;
  valid &= p[1] >= limits.minY;
  valid &= p[1] <= limits.maxY;

  return valid;
}

bool RandomGenerator<Point2DDubins<double>>::isInLimits(Point2DDubins<double>& p) {
  bool valid{true};
  valid &= p[0] >= limits.minX;
  valid &= p[0] <= limits.maxX;
  valid &= p[1] >= limits.minY;
  valid &= p[1] <= limits.maxY;

  return valid;
}

bool RandomGenerator<Point3D<double>>::isInLimits(Point3D<double>& p) {
  bool valid{true};
  valid &= p[0] >= limits.minX;
  valid &= p[0] <= limits.maxX;
  valid &= p[1] >= limits.minY;
  valid &= p[1] <= limits.maxY;
  valid &= p[2] >= limits.minZ;
  valid &= p[2] <= limits.maxZ;

  return valid;
}
