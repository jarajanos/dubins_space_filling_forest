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

template<>
void Solver<Point3DDubins>::initNeighboringMatrix() {
  this->neighboringMatrix = DistanceMatrix<DistanceHolder<Point3DDubins>>(problem.GetNumRoots(), problem.DubinsResolution);
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
bool Solver<Point3DDubins>::isPathFree(const Point3DDubins start, const Point3DDubins finish) {
  opendubins::State3D startDub{start[0], start[1], start[2], start.GetHeading(), start.GetPitch()};
  opendubins::State3D finishDub{finish[0], finish[1], finish[2], finish.GetHeading(), finish.GetPitch()};
  opendubins::Dubins3D pathDub{startDub, finishDub, this->problem.DubinsRadius, -this->problem.MaxPitch, this->problem.MaxPitch};
  double distance{pathDub.length};
  double parts{distance / problem.CollisionDist};
  bool isFree{true};

  for (int index{1}; index < parts && isFree; ++index) {
    opendubins::State3D temp{pathDub.getState(index * distance / parts)};
    Point3DDubins position{temp};
    
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

template<>
void Solver<Point3DDubins>::getAllPaths() {
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

          DistanceHolder<Point3DDubins> &holder1{this->neighboringMatrix(id1, id3, angle1, angle3)};
          for (int j{0}; j < this->connectedTrees.size(); ++j) {
            int id2{connectedTrees[j]->Root->ID};
            if (i == j) {
              continue;
            }
            for (int angle2{0}; angle2 < numAngles; ++angle2) {
              if (!this->neighboringMatrix.Exists(id3, id2, angle3, angle2)) {
                continue;
              }

              DistanceHolder<Point3DDubins> &holder2{this->neighboringMatrix(id3, id2, angle3, angle2)};

              Node<Point3DDubins> *node1;
              Node<Point3DDubins> *node2;
              std::deque<Point3DDubins> plan1;
              std::deque<Point3DDubins> plan2;
              std::deque<Point3DDubins> finalPlan;

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

              Point3DDubins last{plan1.back()};   // root of id3
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
                DistanceHolder<Point3DDubins> newLink{node1, node2, distance, finalPlan};
                this->neighboringMatrix.AddLink(newLink, id1, id2, angle1, angle2, true);
              } else {
                DistanceHolder<Point3DDubins> &link{this->neighboringMatrix(id1, id2, angle1, angle2)};
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

template <> 
void Solver<Point3DDubins>::computeTsp() {
  INFO("Computing TSP");
  TSPMatrix<Point3DDubins> gatsp{this->problem, this->neighboringMatrix};
  TSPMatrix<Point3DDubins> atsp{gatsp.TransformGATSPtoATSP()};

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
              opendubins::State finishDub{node.Position[0], node.Position[1], node.Position.GetHeading()};
              opendubins::State startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position.GetHeading()};
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
                opendubins::State finishDub{node.Position[0], node.Position[1], node.Position.GetHeading()};
                opendubins::State startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position.GetHeading() + (2 * angle * M_PI) / problem.DubinsResolution};
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
void Solver<Point3DDubins>::saveTrees(const FileStruct file) {
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
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;
      fileStream << "o Trees\n";
      for (int i{0}; i < this->allNodes.size(); ++i) {
        Node<Point3DDubins> &node{*this->allNodes[i]};
        if (node.IsRoot()) {
          continue;
        }
        
        if (!node.Closest->IsRoot()) {
          opendubins::State3D finishDub{node.Position[0], node.Position[1], node.Position[2], node.Position.GetHeading(), node.Position.GetPitch()};
          opendubins::State3D startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position[2], node.Closest->Position.GetHeading(), node.Closest->Position.GetPitch()};
          opendubins::Dubins3D pathFromClosest{startDub, finishDub, this->problem.DubinsRadius, -this->problem.MaxPitch, this->problem.MaxPitch};
          
          Point3DDubins lastPoint{startDub};
          Point3DDubins actPoint;
          double length{pathFromClosest.length};
          double parts{length / this->problem.CollisionDist};
          unsigned startingInd{vertexInd};
          for (int index{1}; index < parts; ++index) {
            actPoint = Point3DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
            fileStream << "v" << DELIMITER_OUT;
            Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
            temp.PrintPosition(fileStream);
            fileStream << "\n";
            ++vertexInd;
          }
          actPoint = Point3DDubins(finishDub);
          fileStream << "v" << DELIMITER_OUT;
          Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
          temp.PrintPosition(fileStream);
          fileStream << "\n";
          ++vertexInd;

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
        } else {
          for (auto &angle : node.GetExpandedAngles()) {
            opendubins::State3D finishDub{node.Position[0], node.Position[1], node.Position[2], node.Position.GetHeading(), node.Position.GetPitch()};
            opendubins::State3D startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position[2], node.Closest->Position.GetHeading() + (2 * angle * M_PI) / problem.DubinsResolution, node.Closest->Position.GetPitch()};
            opendubins::Dubins3D pathFromClosest{startDub, finishDub, this->problem.DubinsRadius, -this->problem.MaxPitch, this->problem.MaxPitch};
            
            Point3DDubins lastPoint{startDub};
            Point3DDubins actPoint;
            double length{pathFromClosest.length};
            double parts{length / this->problem.CollisionDist};
            unsigned startingInd{vertexInd};
            for (int index{1}; index < parts; ++index) {
              actPoint = Point3DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
              fileStream << "v" << DELIMITER_OUT;
              Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
              temp.PrintPosition(fileStream);
              fileStream << "\n";
              ++vertexInd;
            }
            actPoint = Point3DDubins(finishDub);
            fileStream << "v" << DELIMITER_OUT;
            Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
            temp.PrintPosition(fileStream);
            fileStream << "\n";
            ++vertexInd;

            vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
          }
        }
      }

      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<Point3DDubins> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            if (!node.Closest->IsRoot()) {
              opendubins::State3D finishDub{node.Position[0], node.Position[1], node.Position[2], node.Position.GetHeading(), node.Position.GetPitch()};
              opendubins::State3D startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position[2], node.Closest->Position.GetHeading(), node.Closest->Position.GetPitch()};
              opendubins::Dubins3D pathFromClosest{startDub, finishDub, this->problem.DubinsRadius, -this->problem.MaxPitch, this->problem.MaxPitch};
              
              Point3DDubins lastPoint{startDub};
              Point3DDubins actPoint;
              double length{pathFromClosest.length};
              double parts{length / this->problem.CollisionDist};
              for (int index{1}; index < parts; ++index) {
                actPoint = Point3DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
                fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
                lastPoint = actPoint;
              }
              actPoint = Point3DDubins(finishDub);
              fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
            } else {
              for (auto &angle : node.GetExpandedAngles()) {
                opendubins::State3D finishDub{node.Position[0], node.Position[1], node.Position[2], node.Position.GetHeading(), node.Position.GetPitch()};
                opendubins::State3D startDub{node.Closest->Position[0], node.Closest->Position[1], node.Closest->Position[2], node.Closest->Position.GetHeading()  + (2 * angle * M_PI) / problem.DubinsResolution, node.Closest->Position.GetPitch()};
                opendubins::Dubins3D pathFromClosest{startDub, finishDub, this->problem.DubinsRadius, -this->problem.MaxPitch, this->problem.MaxPitch};
              
                Point3DDubins lastPoint{startDub};
                Point3DDubins actPoint;
                double length{pathFromClosest.length};
                double parts{length / this->problem.CollisionDist};
                for (int index{1}; index < parts; ++index) {
                  actPoint = Point3DDubins(pathFromClosest.getState(index * this->problem.CollisionDist));
                  fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
                  lastPoint = actPoint;
                }
                actPoint = Point3DDubins(finishDub);
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
                Point2DDubins actPoint{plan[m]};
                Point2DDubins lastPoint{plan[m + 1]};
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

template<>
void Solver<Point3DDubins>::savePaths(const FileStruct file) {
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
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;
      fileStream << "o Paths\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{0}; j < numRoots; ++j) {
          for (int k{0}; k < numAngles; ++k) {
            for (int l{0}; l < numAngles; ++l) {
              if (!this->neighboringMatrix.Exists(i, j, k, l)) {
                continue;  
              }
              DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(i, j, k, l)};
              std::deque<Point3DDubins> &plan{holder.Plan};
              unsigned startingInd{vertexInd};
              for (int m{0}; m < plan.size(); ++m) {
                Point3DDubins actPoint{plan[m]};

                fileStream << "v" << DELIMITER_OUT;
                Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
                temp.PrintPosition(fileStream);
                fileStream << "\n";
              }
              vertexInd += plan.size();
              vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
            }
          }
        }
      }

      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{0}; j < numRoots; ++j) {
          for (int k{0}; k < numAngles; ++k) {
            for (int l{0}; l < numAngles; ++l) {
              if (!this->neighboringMatrix.Exists(i, j, k, l)) {
                continue;  
              }
              DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(i, j, k, l)};
              std::deque<Point3DDubins> &plan{holder.Plan};
              for (int m{0}; m < plan.size() - 1; ++m) {
                Point3DDubins actPoint{plan[m]};
                Point3DDubins lastPoint{plan[m + 1]};

                fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";   
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
          Point2DDubins actPoint{plan[m]};
          Point2DDubins lastPoint{plan[m + 1]};

          fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";
          lastPoint = actPoint;
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
void Solver<Point3DDubins>::saveTspPaths(const FileStruct file) {
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
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;
      fileStream << "o Paths\n";
      for (int i{0}; i < numRoots; ++i) {
        auto [ actNode, actAngle ] = gatspSolution[i];
        auto [ nextNode, nextAngle ] = gatspSolution[(i + 1) % numRoots];
        if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
        std::deque<Point3DDubins> &plan{holder.Plan};
        unsigned startingInd{vertexInd};
        for (int m{0}; m < plan.size(); ++m) {
          Point3DDubins actPoint{plan[m]};

          fileStream << "v" << DELIMITER_OUT;
          Point3DDubins temp{actPoint / problem.Env.ScaleFactor}; 
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        }
        vertexInd += plan.size();
        vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));    
      }
      
      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#TspPaths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        auto [ actNode, actAngle ] = gatspSolution[i];
        auto [ nextNode, nextAngle ] = gatspSolution[(i + 1) % numRoots];
        if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
        std::deque<Point3DDubins> &plan{holder.Plan};
        for (int m{0}; m < plan.size() - 1; ++m) {
          Point3DDubins actPoint{plan[m]};
          Point3DDubins lastPoint{plan[m + 1]};

          fileStream << actPoint / problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / problem.Env.ScaleFactor << "\n";   
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
void Solver<Point2DDubins>::saveTsp(const FileStruct file) {
  ERROR("Saving DTSP file is not supported");
}

template <>
void Solver<Point3DDubins>::saveTsp(const FileStruct file) {
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
    bool first{true};
    for (int i{0}; i < this->connected.size(); ++i) {
        if (!first) {
          fileStream << CSV_DELIMITER_2;
        }
        
        if (this->connected[i]) {
            fileStream << i + 1;
            first = false;
        }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
  
    // do not print distance matrix as for non-dubins case -- would be too big (except single-goal)
    if (!this->problem.HasGoal || !this->neighboringMatrix.Exists(0, 1, 0, 0)) {
      fileStream << CSV_NO_PATH;
    } else {
      double dist{this->neighboringMatrix(0, 1, 0, 0).Distance};
      if (dist == std::numeric_limits<double>::max()) {
        fileStream << CSV_NO_PATH;
      } else {
        fileStream << dist / problem.Env.ScaleFactor;
      }
    }

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

template<>
void Solver<Point3DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
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
    bool first{true};
    for (int i{0}; i < this->connected.size(); ++i) {
        if (!first) {
          fileStream << CSV_DELIMITER_2;
        }
        
        if (this->connected[i]) {
            fileStream << i + 1;
            first = false;
        }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
  
    // do not print distance matrix as for non-dubins case -- would be too big (except single-goal)
    if (!this->problem.HasGoal || !this->neighboringMatrix.Exists(0, 1, 0, 0)) {
      fileStream << CSV_NO_PATH;
    } else {
      double dist{this->neighboringMatrix(0, 1, 0, 0).Distance};
      if (dist == std::numeric_limits<double>::max()) {
        fileStream << CSV_NO_PATH;
      } else {
        fileStream << dist / problem.Env.ScaleFactor;
      }
    }

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
