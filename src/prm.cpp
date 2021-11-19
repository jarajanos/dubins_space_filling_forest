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
      if (!holder.Exists()) {
        continue;
      }
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

template <> 
void ProbRoadMaps<Point2DDubins>::saveTrees(const FileStruct file) {
  INFO("Saving trees");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    if (file.type == Obj) {
      ERROR("Not implemented");
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 
          opendubins::State startDub{node.Position[0], node.Position[1], node.Position.GetAngle()};
          opendubins::State finishDub{neigh->Position[0], neigh->Position[1], neigh->Position.GetAngle()};
          opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
          
          Point2DDubins lastPoint{startDub};
          Point2DDubins actPoint;
          double length{pathFromClosest.length};
          double parts{length / this->problem.CollisionDist};
          for (int index{1}; index < parts; ++index) {
            actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
            fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.ID << DELIMITER_OUT << node.GetAge() << "\n";
            lastPoint = actPoint;
          }
          actPoint = Point2DDubins(finishDub);
          fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.ID << DELIMITER_OUT << node.GetAge() << "\n";
        }
      }
    } else {
      throw std::string("Unimplemented file type");
    }

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}
