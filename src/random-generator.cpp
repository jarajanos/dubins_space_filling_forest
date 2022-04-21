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

// Marsaglia, G. (1972), "Choosing a Point from the Surface of a Sphere", 
template<>
bool RandomGenerator<Point3D>::RandomPointInDistance(const Point3D& center, Point3D& point, const double distance) {
  Point3D temp;

  // rotation
  bool inLimits{false};
  double s, sigOne, sigTwo, thetaOne, thetaTwo;
  int iter{0};
  while (!inLimits && iter < maxIter) {
    double u{normProb(rndEng)};
    double v{normProb(rndEng)};
    double w{normProb(rndEng)};

    double norm{sqrt(u*u + v*v + w*w)};
    temp.SetPosition(center[0] + u * distance / norm, center[1] + v * distance / norm, center[2] + w * distance / norm);

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
    ++iter;
  }

  return inLimits;
}

// Marsaglia, G. (1972), "Choosing a Point from the Surface of a Sphere", 
template<>
bool RandomGenerator<Point3DDubins>::RandomPointInDistance(const Point3DDubins& center, Point3DDubins& point, const double distance) {
  Point3DDubins temp;

  double dist{std::numeric_limits<double>::max()};
  opendubins::Dubins3D dubPath;
  opendubins::State3D a{center[0], center[1], center[2], center.GetHeading(), center.GetPitch()};
  opendubins::State3D b;
  int iter{0};
  while (dist == std::numeric_limits<double>::max() && iter < maxIter) {
    double u{normProb(rndEng)};
    double v{normProb(rndEng)};
    double w{normProb(rndEng)};

    double norm{sqrt(u*u + v*v + w*w)};
    temp.SetPosition(center[0] + u * distance / norm, center[1] + v * distance / norm, center[2] + w * distance / norm);

    // rotation
    temp.SetHeading(uniDistAngle(rndEng));
    temp.SetPitch(uniDistAngle(rndEng));
    
    b = opendubins::State3D{temp[0], temp[1], temp[2], temp.GetHeading(), temp.GetPitch()};
    dubPath = opendubins::Dubins3D{a, b, Point3DDubins::DubinsRadius, -Point3DDubins::MaxPitch, Point3DDubins::MaxPitch};
    dist = dubPath.getLength();
    ++iter;
  }

  if (iter != maxIter) {
    // get point in exact distance
    point = Point3DDubins(dubPath.getState(distance));
    return limits.IsInLimits(point);
  }

  return false;
}

template<>
bool RandomGenerator<Point3DPolynom>::RandomPointInDistance(const Point3DPolynom& center, Point3DPolynom& point, const double distance) {
  Vec3 pos0{center[0], center[1], center[2]};
  Vec3 vel0{center[3], center[4], center[5]};
  Vec3 acc0{center[6], center[7], center[8]};
  Vec3 gravity{0,0,-this->problem.Gravity};

  double dist{std::numeric_limits<double>::max()};
  Point3DPolynom p;
  int iter{0};
  while(dist == std::numeric_limits<double>::max() && iter < maxIter) {
    ++iter;

    double u{normProb(rndEng)};
    double v{normProb(rndEng)};
    double w{normProb(rndEng)};

    double norm{sqrt(u*u + v*v + w*w)};
    Vec3 position{center[0] + u * distance / norm, center[1] + v * distance / norm, center[2] + w * distance / norm};
    p = Point3DPolynom(position);
    if (!limits.IsInLimits(p)) {
      continue;
    }

    RapidTrajectoryGenerator traj(pos0, vel0, acc0, gravity);
    traj.SetGoalPosition(position);
    if (!this->problem.FreeSampling) {
      std::normal_distribution<double> velDist{problem.AvgVelocity, problem.AvgVelocity * AVG_VEL_STD_MULT};
      std::uniform_int_distribution<int> signDist{0, 1};
      Vec3 velocity{TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng)};
      p.SetVelocity(velocity.x, velocity.y, velocity.z);
      traj.SetGoalVelocity(velocity);

      std::uniform_real_distribution<double> accDist{problem.MinThrust, problem.MaxThrust};
      Vec3 acceleration{accDist(rndEng), accDist(rndEng), accDist(rndEng)};
      p.SetAcceleration(acceleration.x, acceleration.y, acceleration.z);
      traj.SetGoalAcceleration(acceleration);
    }

    traj.Generate(this->problem.SegmentTime);
    if (this->problem.FreeSampling) {
      Vec3 finalVel{traj.GetVelocity(this->problem.SegmentTime)};
      Vec3 finalAcc{traj.GetAcceleration(this->problem.SegmentTime)};

      p.SetVelocity(finalVel.x, finalVel.y, finalVel.z);
      p.SetAcceleration(finalAcc.x, finalAcc.y, finalAcc.z);
    }

    dist = traj.GetCost();
  }

  if (iter != maxIter) {
    point = p;
    return true;
  }

  return false;
}

template<>
bool RandomGenerator<Point2DPolynom>::RandomPointInDistance(const Point2DPolynom& center, Point2DPolynom& point, const double distance) {
  Vec2 pos0{center[0], center[1]};
  Vec2 vel0{center[2], center[3]};
  Vec2 acc0{center[4], center[5]};
  Vec2 gravity{0,-this->problem.Gravity};

  double dist{std::numeric_limits<double>::max()};
  Point2DPolynom p;
  int iter{0};
  while(dist == std::numeric_limits<double>::max() && iter < maxIter) {
    ++iter;

    double u{normProb(rndEng)};
    double v{normProb(rndEng)};

    double norm{sqrt(u*u + v*v)};
    Vec2 position{center[0] + u * distance / norm, center[1] + v * distance / norm};
    p = Point2DPolynom(position);
    if (!limits.IsInLimits(p)) {
      continue;
    }

    RapidTrajectoryGenerator2D traj(pos0, vel0, acc0, gravity);
    traj.SetGoalPosition(position);
    if (!this->problem.FreeSampling) {
      std::normal_distribution<double> velDist{problem.AvgVelocity, problem.AvgVelocity * AVG_VEL_STD_MULT};
      std::uniform_int_distribution<int> signDist{0, 1};
      Vec2 velocity{TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng)};
      p.SetVelocity(velocity.x, velocity.y);
      traj.SetGoalVelocity(velocity);

      std::uniform_real_distribution<double> accDist{problem.MinThrust, problem.MaxThrust};
      Vec2 acceleration{accDist(rndEng), accDist(rndEng)};
      p.SetAcceleration(acceleration.x, acceleration.y);
      traj.SetGoalAcceleration(acceleration);
    }

    traj.Generate(this->problem.SegmentTime);
    if (this->problem.FreeSampling) {
      Vec2 finalVel{traj.GetVelocity(this->problem.SegmentTime)};
      Vec2 finalAcc{traj.GetAcceleration(this->problem.SegmentTime)};

      p.SetVelocity(finalVel.x, finalVel.y);
      p.SetAcceleration(finalAcc.x, finalAcc.y);
    }

    dist = traj.GetCost();
  }

  if (iter != maxIter) {
    point = p;
    return true;
  }

  return false;
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
void RandomGenerator<Point2DPolynom>::RandomPointInSpace(Point2DPolynom &point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng));

  if (!this->problem.FreeSampling) {
    std::normal_distribution<double> velDist{problem.AvgVelocity, problem.AvgVelocity * AVG_VEL_STD_MULT};
    std::uniform_int_distribution<int> signDist{0, 1};
    point.SetVelocity(TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng));

    std::uniform_real_distribution<double> accDist{problem.MinThrust, problem.MaxThrust};
    point.SetAcceleration(accDist(rndEng), accDist(rndEng));
  }
}

template<>
void RandomGenerator<Point3D>::RandomPointInSpace(Point3D& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  // rotation
  bool inLimits{false};
  double s, sigOne, sigTwo, thetaOne, thetaTwo;
  int iter{0};
  while (!inLimits && iter < maxIter) {
    s = RandomProbability();
    sigOne = sqrt(1-s);
    sigTwo = sqrt(s);
    thetaOne = 2*M_PI*RandomProbability();
    thetaTwo = 2*M_PI*RandomProbability();  
  
    Quaternion rotation{cos(thetaTwo) * sigTwo, sin(thetaOne) * sigOne, cos(thetaOne) * sigOne, sin(thetaTwo) * sigTwo};
    point.SetRotation(rotation);
    inLimits = limits.IsInLimits(point);
    ++iter;
  }
}

template<>
void RandomGenerator<Point3DDubins>::RandomPointInSpace(Point3DDubins& point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  // rotation
  point.SetHeading(uniDistAngle(rndEng));
  point.SetPitch(uniDistPitch(rndEng));
}

template<>
void RandomGenerator<Point3DPolynom>::RandomPointInSpace(Point3DPolynom &point) {
  point.SetPosition(uniSpaceX(rndEng), uniSpaceY(rndEng), uniSpaceZ(rndEng));

  if (!this->problem.FreeSampling) {
    std::normal_distribution<double> velDist{problem.AvgVelocity, problem.AvgVelocity * AVG_VEL_STD_MULT};
    std::uniform_int_distribution<int> signDist{0, 1};
    point.SetVelocity(TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng), TO_SIGN(signDist(rndEng)) * velDist(rndEng));

    std::uniform_real_distribution<double> accDist{problem.MinThrust, problem.MaxThrust};
    point.SetAcceleration(accDist(rndEng), accDist(rndEng), accDist(rndEng));
  }
}
