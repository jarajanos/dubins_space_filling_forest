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
void Solver<Point2DDubins>::getAllPaths() {
  int numRoots{this->problem.GetNumRoots()};
  int numAngles{this->problem.DubinsResolution};
  for (int k{0}; k < this->connectedTrees.size(); ++k) {
    int id3{connectedTrees[k]->Root->ID};
    for (int angle3{0}; angle3 < numAngles; ++angle3) {
      for (int i{0}; i < this->connectedTrees.size(); ++i) {
        int id1{connectedTrees[i]->Root->ID};
        if (i == k) {
          continue;
        }

        for (int angle1{0}; angle1 < numAngles; ++angle1) {
          if (!this->neighboringMatrix.Exists(id1, id3, angle1, angle3)) {
            continue;
          }

          DistanceHolder<Point2DDubins> &holder1{this->neighboringMatrix(id1, id3, angle1, angle3)};
          for (int j{0}; j < this->connectedTrees.size(); ++j) {
            int id2{connectedTrees[j]->Root->ID};
            if (i == j) {
              continue;
            }
            for (int angle2{0}; angle2 < numAngles; ++angle2) {
              if (!this->neighboringMatrix.Exists(id3, id2, angle3, angle2)) {
                continue;
              }

              DistanceHolder<Point2DDubins> &holder2{this->neighboringMatrix(id3, id2, angle3, angle2)};

              Node<Point2DDubins> *node1;
              Node<Point2DDubins> *node2;
              std::deque<Point2DDubins> plan1;
              std::deque<Point2DDubins> plan2;
              std::deque<Point2DDubins> finalPlan;

              bool reversed1{false};
              bool reversed2{false};
              plan1 = holder1.Plan;
              if (holder1.Node1->SourceTree->Root->ID == id1) {
                node1 = holder1.Node1;
              } else {
                node1 = holder1.Node2;
                std::reverse(plan1.begin(), plan1.end());
                reversed1 = true;
              }

              plan2 = holder2.Plan;
              if (holder2.Node2->SourceTree->Root->ID == id2) {
                node2 = holder2.Node2;
              } else {
                node2 = holder2.Node1;
                std::reverse(plan2.begin(), plan2.end());
                reversed2 = true;
              }

              Point2DDubins last{plan1.back()};   // root of id3
              plan1.pop_back();
              plan2.pop_front();
              while (!plan1.empty() && !plan2.empty() && plan1.back() == plan2.front().GetInvertedPoint()) {
                last = plan1.back();
                plan1.pop_back();
                plan2.pop_front();
              }

              if (!this->isPathFree(last, plan2.front())) {
                continue; // path is blocked, shortcut not possible TODO: let it sink and think ? maybe ?
              }

              while (plan1.size() > 0) {
                finalPlan.push_back(plan1.front());
                plan1.pop_front();
              }
              finalPlan.push_back(last);
              while (plan2.size() > 0) {
                finalPlan.push_back(plan2.front());
                plan2.pop_front();
              }

              double distance{computeDistance(finalPlan)};
              if (!this->neighboringMatrix.Exists(id1, id2, angle1, angle2)) {
                DistanceHolder<Point2DDubins> newLink{node1, node2, distance, finalPlan};
                this->neighboringMatrix.AddLink(newLink, id1, id2, angle1, angle2, true);
              } else {
                DistanceHolder<Point2DDubins> &link{this->neighboringMatrix(id1, id2, angle1, angle2)};
                if (distance < link.Distance - SFF_TOLERANCE) {
                  link.Node1 = node1;
                  link.Node2 = node2;
                  link.Distance = distance;
                  link.Plan = finalPlan;
                }
              }
            }
          }
        }
      }
    }
  }
}

template <> 
void Solver<Point2DDubins>::computeTsp() {
  INFO("Computing TSP");
  TSPMatrix<Point2DDubins> gatsp{this->problem, this->neighboringMatrix};
  TSPMatrix<Point2DDubins> atsp{gatsp.TransformGATSPtoATSP()};

  // converted to ATSP, run desired solver
  TSPOrder tempSol;
  if (this->problem.TspType == Concorde) {
    tempSol = atsp.SolveByConcorde();
  } else if (this->problem.TspType == LKH) {
    tempSol = atsp.SolveByLKH();
  } else {
    ERROR("TSP solver not implemented");
    exit(1);
  }

  if (tempSol.size() == 0) {
    WARN("TSP not solved");
    if (SaveTSPPaths <= this->problem.SaveOpt) {
      this->problem.SaveOpt = this->problem.SaveOpt - SaveTSPPaths;
    }
  }

  gatspSolution = atsp.TransformATSPSolToGATSP(tempSol);

  // convert also to tsp solution
  for (auto &pair : gatspSolution) {
    auto [ nodeId, angleId ] = pair;
    tspSolution.push_back(nodeId);
  }
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
            if (!node.Closest->IsRoot()) {
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
            } else {
              for (auto &angle : node.GetExpandedAngles()) {
                opendubins::State finishDub{node.Position[0], node.Position[1], node.Position.GetAngle()};
                opendubins::State startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position.GetAngle() + (2 * angle * M_PI) / problem.DubinsResolution};
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

template<>
void Solver<Point2DDubins>::savePaths(const FileStruct file) {
  INFO("Saving paths");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    int numRoots{this->problem.GetNumRoots()};
    int numAngles{this->problem.DubinsResolution};
    if (file.type == Obj) {
      ERROR("Not implemented yet");
      // fileStream << "o Paths\n";
      // for (int i{0}; i < this->allNodes.size(); ++i) {
      //   R temp{this->allNodes[i]->Position / problem.Env.ScaleFactor};
      //   fileStream << "v" << DELIMITER_OUT;
      //   temp.PrintPosition(fileStream);
      //   fileStream << "\n";
      // }
      
      // for (int i{0}; i < numRoots; ++i) {
      //   for (int j{i + 1}; j < numRoots; ++j) {
      //     DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
      //     if (holder.Node1 == NULL) {
      //       continue;
      //     }

      //     std::deque<Node<R> *> &plan{holder.Plan};
      //     for (int k{0}; k < plan.size() - 1; ++k) {
      //       fileStream << "l" << DELIMITER_OUT << plan[k]->ID + 1 << DELIMITER_OUT << plan[k+1]->ID + 1 << "\n";
      //     }
      //   }
      // }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{0}; j < numRoots; ++j) {
          for (int k{0}; k < numAngles; ++k) {
            for (int l{0}; l < numAngles; ++l) {
              if (!this->neighboringMatrix.Exists(i, j, k, l)) {
                continue;  
              }
              DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(i, j, k, l)};
              std::deque<Point2DDubins> &plan{holder.Plan};
              for (int m{0}; m < plan.size() - 1; ++m) {
                opendubins::State startDub{plan[m][0], plan[m][1], plan[m].GetAngle()};
                opendubins::State finishDub{plan[m+1][0], plan[m+1][1], plan[m+1].GetAngle()};
                opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
                
                Point2DDubins lastPoint{startDub};
                Point2DDubins actPoint;
                double length{pathFromClosest.length};
                double parts{length / this->problem.CollisionDist};
                for (int index{1}; index < parts; ++index) {
                  actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
                  fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
                  lastPoint = actPoint;
                }
                actPoint = Point2DDubins(finishDub);
                fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
              }
              fileStream << "\n";    
            }
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


template <> 
void Solver<Point2DDubins>::saveTspPaths(const FileStruct file) {
  INFO("Saving TSP paths");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    int numRoots{this->problem.GetNumRoots()};
    int numAngles{this->problem.DubinsResolution};
    if (file.type == Obj) {
      ERROR("Not implemented yet");
    } else if (file.type == Map) {
      fileStream << "#TspPaths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        auto [ actNode, actAngle ] = gatspSolution[i];
        auto [ nextNode, nextAngle ] = gatspSolution[(i + 1) % numRoots];
        if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
        std::deque<Point2DDubins> &plan{holder.Plan};
        for (int m{0}; m < plan.size() - 1; ++m) {
          opendubins::State startDub{plan[m][0], plan[m][1], plan[m].GetAngle()};
          opendubins::State finishDub{plan[m+1][0], plan[m+1][1], plan[m+1].GetAngle()};
          opendubins::Dubins pathFromClosest{startDub, finishDub, this->problem.DubinsRadius};
          
          Point2DDubins lastPoint{startDub};
          Point2DDubins actPoint;
          double length{pathFromClosest.length};
          double parts{length / this->problem.CollisionDist};
          for (int index{1}; index < parts; ++index) {
            actPoint = Point2DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
            fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
            lastPoint = actPoint;
          }
          actPoint = Point2DDubins(finishDub);
          fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
        }
        fileStream << "\n";    
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

template <>
void Solver<Point2DDubins>::saveTsp(const FileStruct file) {
  ERROR("Saving DTSP file is not supported");
}

template<>
void Solver<Point2DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
  INFO("Saving parameters");
  std::ofstream fileStream{file.fileName.c_str(), std::ios_base::openmode::_S_app};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << this->problem.ID << CSV_DELIMITER;
    fileStream << this->problem.Repetition << CSV_DELIMITER;
    fileStream << iterations << CSV_DELIMITER;
    fileStream << (solved ? "solved" : "unsolved") << CSV_DELIMITER;
  
    fileStream << "[";
    for (int i{0}; i < this->connectedTrees.size(); ++i) {
      fileStream << this->connectedTrees[i]->Root->ID;
      if (i + 1 != connectedTrees.size()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
  
    // do not print distance matrix as for non-dubins case -- would be too big
    fileStream << CSV_NO_PATH;

    fileStream << "]" << CSV_DELIMITER;
    fileStream << elapsedTime.count() << "\n";

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}
