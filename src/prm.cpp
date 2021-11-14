/**
 * @file prm.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 1.0
 * @date 12. 11. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "prm.h"

template <>
void ProbRoadMaps<Point2DDubins>::getPaths() {
  //TODO
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Dijkstra<Point2DDubins> dijkstra;
    std::vector<int> goals;
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (i == j) {
        continue;
      } 
      goals.push_back(j);
    }
    std::deque<DistanceHolder<Point2DDubins>> plans{dijkstra.findPath(i, goals, allPoints)};
    for (auto &holder : plans) {
      int id1{holder.Node1->ID};
      int id2{holder.Node2->ID};
      this->neighboringMatrix.AddLink(holder, id1, id2, 0, 0);
    }
  }
}

template <> 
void ProbRoadMaps<Point2DDubins>::getConnected() {
  //TODO
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{i}; j < this->problem.GetNumRoots(); ++j) {
      if (!this->neighboringMatrix.Exists(i, j, 0, 0)) {
        continue;
      }

      if (!this->connected[i]) {
        this->connected[i] = true;
      }

      if (!this->connected[j]) {
        this->connected[j] = true;
      }
    }
  }
}