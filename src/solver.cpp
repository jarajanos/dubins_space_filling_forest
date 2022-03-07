/**
 * @file solver.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "solver.h"


/**
 * @brief Optimized versions of local planners (might be also solved using GetStateInDistance which however results in more operations)
 * 
 * @tparam  Type of point/problem
 * @param start Starting point for the local planner
 * @param finish End point for the local planner
 * @return true When the path between start and finish is clear for the robot
 * @return false When the path between start and finish is blocked for the robot
 */
template<>
bool SolverBase<Point2D>::isPathFree(const Point2D start, const Point2D finish) {
  double distance{start.Distance(finish)};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};
  PointVector2D direction{start, finish};

  Point2D position;
  for (int index{1}; index < parts && isFree; ++index) {
    position.SetPosition(start[0] + direction[0] * index / parts, start[1] + direction[1] * index / parts);

    isFree &= !problem.Env.Collide(position);
  }

  return isFree;
}

template<>
bool SolverBase<Point2DDubins>::isPathFree(const Point2DDubins start, const Point2DDubins finish) {
  opendubins::State startDub{start[0], start[1], start.GetHeading()};
  opendubins::State finishDub{finish[0], finish[1], finish.GetHeading()};
  opendubins::Dubins pathDub{startDub, finishDub, this->problem.DubinsRadius};
  double distance{pathDub.length};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State temp{pathDub.getState(index * distance / parts)};
    Point2DDubins position{temp};
    
    isFree &= !problem.Env.Collide(position);
    isFree &= problem.Env.Limits.IsInLimits(position);
  }

  return isFree;
}

template<>
bool SolverBase<Point3D>::isPathFree(const Point3D start, const Point3D finish) {
  double distance{start.Distance(finish)};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};
  PointVector3D direction{start, finish};

  Quaternion firstAngle{start.GetRotation()};
  Quaternion secondAngle{finish.GetRotation()};
  double product{firstAngle | secondAngle}; // inner product of quaternions
  if (product < 0) {
    secondAngle = secondAngle.Inverse();
    secondAngle.Set(0, -secondAngle[0]);
  }

  Point3D position;
  for (int index{1}; index < parts && isFree; ++index) {
    position.SetPosition(start[0] + direction[0] * index / parts, start[1] + direction[1] * index / parts, start[2] + direction[2] * index / parts);
    
    Quaternion rot{firstAngle[0] + index * (secondAngle[0] - firstAngle[0]) / parts, firstAngle[1] + index * (secondAngle[1] - firstAngle[1]) / parts, firstAngle[2] + index * (secondAngle[2] - firstAngle[2]) / parts, firstAngle[3] + index * (secondAngle[3] - firstAngle[3]) / parts};
    rot.Normalize();
    position.SetRotation(rot);

    isFree &= !problem.Env.Collide(position);
  }

  return isFree;
}

template<>
bool SolverBase<Point3DDubins>::isPathFree(const Point3DDubins start, const Point3DDubins finish) {
  opendubins::State3D startDub{start[0], start[1], start[2], start.GetHeading(), start.GetPitch()};
  opendubins::State3D finishDub{finish[0], finish[1], finish[2], finish.GetHeading(), finish.GetPitch()};
  opendubins::Dubins3D pathDub{startDub, finishDub, this->problem.DubinsRadius, this->problem.Env.Limits.mins[3], this->problem.Env.Limits.maxs[3]};
  double distance{pathDub.length};
  if (distance == std::numeric_limits<double>::max()) {
    return false;
  }

  double parts{distance / problem.CollisionDist};
  bool isFree{true};

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State3D temp{pathDub.getState(index * distance / parts)};
    Point3DDubins position{temp};
    position.SetPitch(0); // fixes problem with opendubins
    
    isFree &= !problem.Env.Collide(position);
    isFree &= problem.Env.Limits.IsInLimits(position);
  }

  return isFree;
}

template<>
bool SolverBase<Point3DPolynom>::isPathFree(const Point3DPolynom start, const Point3DPolynom finish) {
  Vec3 startPos{start[0], start[1], start[2]};
  Vec3 startVel{start[3], start[4], start[5]};
  Vec3 startAcc{start[6], start[7], start[8]};

  Vec3 finishPos{finish[0], finish[1], finish[2]};
  Vec3 finishVel{finish[3], finish[4], finish[5]};
  Vec3 finishAcc{finish[6], finish[7], finish[8]};

  Vec3 gravity{0,0,-GRAVITY};

  RapidTrajectoryGenerator traj(startPos, startVel, startAcc, gravity);
  traj.SetGoalPosition(finishPos);
  traj.SetGoalVelocity(finishVel);
  traj.SetGoalAcceleration(finishAcc);

  traj.Generate(this->problem.SegmentTime);

  double parts{this->problem.SegmentTime / this->problem.CtrlInterval};
  bool isFree{true};
  for (int index{0}; index < parts && isFree; ++index) {
    Vec3 pos{traj.GetPosition(index * this->problem.CtrlInterval)};
    Point3DPolynom pos3D{pos};

    isFree &= !problem.Env.Collide(pos3D);
  }
  traj.GetPosition(this->problem.CtrlInterval);

  return isFree;
}

// approximation for getAllPaths
template <>
double Solver<Point2DDubins>::computeDistance(std::deque<Point2DDubins> &plan) {
  Point2D actual; 
  Point2D previous;
  bool first{true};
  double distance{0};
  for (Point2DDubins pos : plan) {
    actual = Point2D(pos[0], pos[1]);
    if (!first) {
      distance += previous.Distance(actual);
    } else {
      first = false;
    }
    
    previous = actual;
  }

  return distance;
}

// approximation for getAllPaths
template <>
double Solver<Point3DDubins>::computeDistance(std::deque<Point3DDubins> &plan) {
  Point3D actual; 
  Point3D previous;
  bool first{true};
  double distance{0};
  for (Point3DDubins pos : plan) {
    actual = Point3D(pos[0], pos[1], pos[2], pos.GetHeading(), pos.GetPitch(), 0);
    if (!first) {
      distance += previous.Distance(actual);
    } else {
      first = false;
    }
    
    previous = actual;
  }

  return distance;
}
