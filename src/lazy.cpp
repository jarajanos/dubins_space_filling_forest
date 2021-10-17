/**
 * @file lazy.cpp
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "lazy.h"

// template<>
// void LazyTSP<Point2DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) {
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
      
//       for (auto &pair : selectedPaths) {
//         int first, second;
//         std::tie(first, second) = pair;
//         DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(first, second)};

//         std::deque<Node<Point2DDubins> *> &plan{holder.Plan};
//         for (int k{0}; k < plan.size() - 1; ++k) {
//           opendubins::State finishDub{plan[k]->Position[0], plan[k]->Position[1], plan[k]->Position.GetAngle()};
//           opendubins::State startDub{plan[k+1]->Position[0], plan[k+1]->Position[1], plan[k+1]->Position.GetAngle()};
//           opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
          
//           Point2DDubins lastPoint{startDub};
//           Point2DDubins actPoint;
//           double length{pathFromClosest.length};
//           double parts{length / this->problem.CollisionDist};
//           for (int index{1}; index < parts; ++index) {
//             actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
//             fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
//             lastPoint = actPoint;
//           }
//           actPoint = Point2DDubins(finishDub);
//           fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
//         }
//         fileStream << "\n";
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
