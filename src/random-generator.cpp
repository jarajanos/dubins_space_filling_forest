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
bool RandomGenerator<Point2D>::RandomPointInDistance(const Point2D& center, Point2D& point, const double distance) {
  // rotation in polar
  double phi{uniDistAngle(rndEng)};
  
  point.SetPosition(center[0] + cos(phi) * distance, center[1] + sin(phi) * distance);
  return limits.IsInLimits(point);
}

template<>
bool RandomGenerator<Point2DDubins>::RandomPointInDistance(const Point2DDubins& center, Point2DDubins& point, const double distance) {
  // rotation in polar
  double phi{uniDistAngle(rndEng)};
  
  point.SetPosition(center[0] + cos(phi) * distance, center[1] + sin(phi) * distance);
  point.SetHeading(uniDistAngle(rndEng));

  opendubins::State a{center[0], center[1], center.GetHeading()};
  opendubins::State b{point[0], point[1], point.GetHeading()};
  opendubins::Dubins dubPath{a, b, Point2DDubins::DubinsRadius};

  // get point in exact distance
  point = Point2DDubins(dubPath.getState(distance));
  return limits.IsInLimits(point);
}

template<>
bool RandomGenerator<Point3D>::RandomPointInDistance(const Point3D& center, Point3D& point, const double distance) {
  Point3D temp;

  double phi{uniDistAngle(rndEng)};
  double theta{uniDistAngle(rndEng)};
  temp.SetPosition(center[0] + cos(theta) * sin(phi) * distance, center[1] + sin(theta) * sin(phi)  * distance, center[2] + cos(phi) * distance);

  // rotation
  bool inLimits{false};
  double s, sigOne, sigTwo, thetaOne, thetaTwo;
  while (!inLimits) {
    s = RandomProbability();
    sigOne = sqrt(1-s);
    sigTwo = sqrt(s);
    thetaOne = 2*M_PI*RandomProbability();
    thetaTwo = 2*M_PI*RandomProbability();  
  
    Quaternion rotation{cos(thetaTwo) * sigTwo, sin(thetaOne) * sigOne, cos(thetaOne) * sigOne, sin(thetaTwo) * sigTwo};
    temp.SetRotation(rotation);

    // the distance is not exactly the expected one (might be much greater), therefore a point in the same direction with the correct distance is needed
    point = center.GetStateInDistance(temp, distance);
    inLimits = limits.IsInLimits(point);
  }

  return inLimits;
}

template<>
bool RandomGenerator<Point3DDubins>::RandomPointInDistance(const Point3DDubins& center, Point3DDubins& point, const double distance) {
  Point3DDubins temp;
  
  double phi{uniDistAngle(rndEng)};
  double theta{uniDistAngle(rndEng)};
  temp.SetPosition(center[0] + cos(theta) * sin(phi) * distance, center[1] + sin(theta) * sin(phi)  * distance, center[2] + cos(phi) * distance);

  // rotation
  temp.SetHeading(uniDistAngle(rndEng));
  temp.SetPitch(uniDistAngle(rndEng));

  opendubins::State3D a{center[0], center[1], center[2], center.GetHeading(), center.GetPitch()};
  opendubins::State3D b{temp[0], temp[1], temp[2], temp.GetHeading(), temp.GetPitch()};
  opendubins::Dubins3D dubPath{a, b, Point3DDubins::DubinsRadius, -Point3DDubins::MaxPitch, Point3DDubins::MaxPitch};

  // get point in exact distance
  point = Point3DDubins(dubPath.getState(distance));
  return limits.IsInLimits(point);
}

/**
 * @brief Random space uniformly sampled in the configuration space of particular problem, according to 
 * https://ri.cmu.edu/pub_files/pub4/kuffner_james_2004_1/kuffner_james_2004_1.pdf 
 * 
 * @param point Output = uniformly sampled point in the configuration space
 */
template<>
void RandomGenerator<Point2D>::RandomPointInSpace(Point2D& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng));
}

template<>
void RandomGenerator<Point2DDubins>::RandomPointInSpace(Point2DDubins& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng));
  point.SetHeading(uniDistAngle(rndEng));
}

template<>
void RandomGenerator<Point3D>::RandomPointInSpace(Point3D& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  // rotation
  bool inLimits{false};
  double s, sigOne, sigTwo, thetaOne, thetaTwo;
  while (!inLimits) {
    s = RandomProbability();
    sigOne = sqrt(1-s);
    sigTwo = sqrt(s);
    thetaOne = 2*M_PI*RandomProbability();
    thetaTwo = 2*M_PI*RandomProbability();  
  
    Quaternion rotation{cos(thetaTwo) * sigTwo, sin(thetaOne) * sigOne, cos(thetaOne) * sigOne, sin(thetaTwo) * sigTwo};
    point.SetRotation(rotation);
    inLimits = limits.IsInLimits(point);
  }
}

template<>
void RandomGenerator<Point3DDubins>::RandomPointInSpace(Point3DDubins& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  // rotation
  point.SetHeading(uniDistAngle(rndEng));
  point.SetPitch(uniDistPitch(rndEng));
}
