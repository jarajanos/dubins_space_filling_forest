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

template<>
void Solver<Point2DDubins>::initNeighboringMatrix() {
  this->neighboringMatrix = DistanceMatrix<DistanceHolder<Point2DDubins>>(problem.GetNumRoots(), problem.DubinsResolution);
}

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
bool Solver<Point2D>::isPathFree(const Point2D start, const Point2D finish) {
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
bool Solver<Point2DDubins>::isPathFree(const Point2DDubins start, const Point2DDubins finish) {
  opendubins::State startDub{start[0], start[1], start.GetAngle()};
  opendubins::State finishDub{finish[0], finish[1], finish.GetAngle()};
  opendubins::Dubins pathDub{startDub, finishDub, this->problem.DubinsRadius};
  double distance{pathDub.length};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State temp{pathDub.getState(index * distance / parts)};
    Point2DDubins position{temp};
    
    isFree &= !problem.Env.Collide(position);
  }

  return isFree;
}

template<>
bool Solver<Point3D>::isPathFree(const Point3D start, const Point3D finish) {
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
void Solver<Point2DDubins>::saveTrees(const FileStruct file) {
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
      ERROR("Not implemented yet");
      // fileStream << "o Trees\n";
      // for (int i{0}; i < this->allNodes.size(); ++i) {
      //   Point2DDubins temp{this->allNodes[i]->Position / problem.Env.ScaleFactor};
      //   fileStream << "v" << DELIMITER_OUT;
      //   temp.PrintPosition(fileStream);
      //   fileStream << "\n";
      // }

      // for (int i{0}; i < this->trees.size(); ++i) {
      //   for (Node<R> &node : this->trees[i].Leaves) {
      //     if (node.DistanceToRoot != 0) {
      //       fileStream << "l" << DELIMITER_OUT << node.ID + 1 << DELIMITER_OUT << node.Closest->ID + 1 << "\n";
      //     }
      //   }
      // }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<Point2DDubins> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            opendubins::State finishDub{node.Position[0], node.Position[1], node.Position.GetAngle()};
            opendubins::State startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position.GetAngle()};
            opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
            
            Point2DDubins lastPoint{startDub};
            Point2DDubins actPoint;
            double length{pathFromClosest.length};
            double parts{length / this->problem.CollisionDist};
            for (int index{1}; index < parts; ++index) {
              actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
              fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
              lastPoint = actPoint;
            }
            actPoint = Point2DDubins(finishDub);
            fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
          }
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

// template<>
// void Solver<Point2DDubins>::savePaths(const FileStruct file) {
//   INFO("Saving paths");
//   std::ofstream fileStream{file.fileName.c_str()};
//   if (!fileStream.good()) {
//     std::stringstream message;
//     message << "Cannot create file at: " << file.fileName;

//     WARN(message.str());
//     return;
//   }

//   if (fileStream.is_open()) {
//     int numRoots{this->problem.GetNumRoots()};
//     if (file.type == Obj) {
//       ERROR("Not implemented yet");
//       // fileStream << "o Paths\n";
//       // for (int i{0}; i < this->allNodes.size(); ++i) {
//       //   R temp{this->allNodes[i]->Position / problem.Env.ScaleFactor};
//       //   fileStream << "v" << DELIMITER_OUT;
//       //   temp.PrintPosition(fileStream);
//       //   fileStream << "\n";
//       // }
      
//       // for (int i{0}; i < numRoots; ++i) {
//       //   for (int j{i + 1}; j < numRoots; ++j) {
//       //     DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
//       //     if (holder.Node1 == NULL) {
//       //       continue;
//       //     }

//       //     std::deque<Node<R> *> &plan{holder.Plan};
//       //     for (int k{0}; k < plan.size() - 1; ++k) {
//       //       fileStream << "l" << DELIMITER_OUT << plan[k]->ID + 1 << DELIMITER_OUT << plan[k+1]->ID + 1 << "\n";
//       //     }
//       //   }
//       // }
//     } else if (file.type == Map) {
//       fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
//       for (int i{0}; i < numRoots; ++i) {
//         for (int j{0}; j < numRoots; ++j) {
//           DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(i, j)};
//           if (holder.Node1 == NULL) {
//             continue;
//           }

//           std::deque<Node<Point2DDubins> *> &plan{holder.Plan};
//           for (int k{0}; k < plan.size() - 1; ++k) {
//             opendubins::State finishDub{plan[k]->Position[0], plan[k]->Position[1], plan[k]->Position.GetAngle()};
//             opendubins::State startDub{plan[k+1]->Position[0], plan[k+1]->Position[1], plan[k+1]->Position.GetAngle()};
//             opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
            
//             Point2DDubins lastPoint{startDub};
//             Point2DDubins actPoint;
//             double length{pathFromClosest.length};
//             double parts{length / this->problem.CollisionDist};
//             for (int index{1}; index < parts; ++index) {
//               actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
//               fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
//               lastPoint = actPoint;
//             }
//             actPoint = Point2DDubins(finishDub);
//             fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
//           }
//           fileStream << "\n";
//         }
//       }
//     } else {
//       throw std::string("Unimplemented file type");
//     }

//     fileStream.flush();
//     fileStream.close();
//   } else {
//     std::stringstream message;
//     message << "Cannot open file at: " << file.fileName;
//     WARN(message.str());
//   }  
// }
