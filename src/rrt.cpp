/**
 * @file rrt.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "rrt.h"

template<> 
double RapidExpTreeBase<Point3DPolynom>::getApproxDistance(Node<Point3DPolynom> &start, Node<Point3DPolynom> &goal) {
  return start.Position.EuclideanDistance(goal.Position);
}

template<> 
double RapidExpTreeBase<Point2DPolynom>::getApproxDistance(Node<Point2DPolynom> &start, Node<Point2DPolynom> &goal) {
  return start.Position.EuclideanDistance(goal.Position);
}
