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
bool Solver<Point2D>::isPathFree(Point2D &start, Point2D &finish) {
  double distance{start.Distance(finish)};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};
  PointVector2D direction{start, finish};

  Point2D position;
  for (int index{1}; index < parts && isFree; ++index) {
    position.SetPosition(start[0] + direction[0] * index / parts, start[1] + direction[1] * index / parts);

    isFree &= problem.Env.Collide(position);
  }

  return isFree;
}

template<>
bool Solver<Point2DDubins>::isPathFree(Point2DDubins &start, Point2DDubins &finish) {
  opendubins::State startDub{start[0], start[1], start.GetAngle()};
  opendubins::State finishDub{finish[0], finish[1], finish.GetAngle()};
  opendubins::Dubins pathDub{startDub, finishDub, this->problem.DubinsRadius};
  double distance{pathDub.length};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State temp{pathDub.getState(index * distance / parts)};
    Point2DDubins position{temp};
    
    isFree &= problem.Env.Collide(position);
  }

  // TODO checking both ways should be perhaps replaced by two point states
  opendubins::Dubins pathDubBack{finishDub, startDub, this->problem.DubinsRadius};
  distance = pathDubBack.length;
  parts = distance / problem.CollisionDist;

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State temp{pathDubBack.getState(index * distance / parts)};
    Point2DDubins position{temp};
    
    isFree &= problem.Env.Collide(position);
  }

  return isFree;
}

template<>
bool Solver<Point3D>::isPathFree(Point3D &start, Point3D &finish) {
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

    isFree &= problem.Env.Collide(position);
  }

  return isFree;
}
