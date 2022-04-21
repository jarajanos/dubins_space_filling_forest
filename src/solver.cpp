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

  Vec3 gravity{0,0,-this->problem.Gravity};

  double totalDist{start.EuclideanDistance(finish)};
  double totalTime{totalDist / this->problem.AvgVelocity};

  RapidTrajectoryGenerator traj(startPos, startVel, startAcc, gravity);
  traj.SetGoalPosition(finishPos);
  traj.SetGoalVelocity(finishVel);
  traj.SetGoalAcceleration(finishAcc);

  traj.Generate(totalTime);

  double parts{totalTime / this->problem.CtrlInterval};
  bool isFree;
  
  RapidTrajectoryGenerator::InputFeasibilityResult result{traj.CheckInputFeasibility(this->problem.MinThrust, this->problem.MaxThrust, this->problem.MaxRotSpeed, this->problem.CtrlInterval)};
  isFree = (result == traj.InputFeasible || result == traj.InputIndeterminable);
  
  for (int index{0}; index < parts && isFree; ++index) {
    Vec3 pos{traj.GetPosition(index * this->problem.CtrlInterval)};
    Point3DPolynom pos3D{pos};

    isFree &= !problem.Env.Collide(pos3D);
    isFree &= problem.Env.Limits.IsInLimits(pos3D);
  }

  return isFree;
}

template<>
bool SolverBase<Point2DPolynom>::isPathFree(const Point2DPolynom start, const Point2DPolynom finish) {
  Vec2 startPos{start[0], start[1]};
  Vec2 startVel{start[2], start[3]};
  Vec2 startAcc{start[4], start[5]};

  Vec2 finishPos{finish[0], finish[1]};
  Vec2 finishVel{finish[2], finish[3]};
  Vec2 finishAcc{finish[4], finish[5]};

  Vec2 gravity{0,-this->problem.Gravity};

  double totalDist{start.EuclideanDistance(finish)};
  double totalTime{totalDist / this->problem.AvgVelocity};

  RapidTrajectoryGenerator2D traj(startPos, startVel, startAcc, gravity);
  traj.SetGoalPosition(finishPos);
  traj.SetGoalVelocity(finishVel);
  traj.SetGoalAcceleration(finishAcc);

  traj.Generate(totalTime);

  double parts{totalTime / this->problem.CtrlInterval};
  bool isFree;
  
  RapidTrajectoryGenerator2D::InputFeasibilityResult result{traj.CheckInputFeasibility(this->problem.MinThrust, this->problem.MaxThrust, this->problem.MaxRotSpeed, this->problem.CtrlInterval)};
  isFree = (result == traj.InputFeasible || result == traj.InputIndeterminable);
  
  for (int index{0}; index < parts && isFree; ++index) {
    Vec2 pos{traj.GetPosition(index * this->problem.CtrlInterval)};
    Point2DPolynom pos2D{pos};

    isFree &= !problem.Env.Collide(pos2D);
    isFree &= problem.Env.Limits.IsInLimits(pos2D);
  }

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

template<>
void Solver<Point3DPolynom, false>::saveTrees(const FileStruct file) {
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
        Node<Point3DPolynom> &node{*(this->allNodes[i])};
        if (node.IsRoot()) {
          continue;
        }

        unsigned startingInd{vertexInd};
        Point3DPolynom actPoint{node.Position};
        auto path{node.Closest->Position.SampleTrajectory(actPoint, this->problem.CtrlInterval)};

        for (auto &point : path) {
          fileStream << "v" << DELIMITER_OUT;
          Point3DPolynom temp{point / this->problem.Env.ScaleFactor}; 
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        }

        vertexInd += path.size();
        vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
      }

      for (auto &tup : vertexRanges) {
        auto [start, end] = tup;
        end -= 1;
        for (int i = start; i < end; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << "\n";
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<Point3DPolynom> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            fileStream << node.Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.Closest->Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
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
void Solver<Point2DPolynom, false>::saveTrees(const FileStruct file) {
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
        Node<Point2DPolynom> &node{*(this->allNodes[i])};
        if (node.IsRoot()) {
          continue;
        }

        unsigned startingInd{vertexInd};
        Point2DPolynom actPoint{node.Position};
        auto path{node.Closest->Position.SampleTrajectory(actPoint, this->problem.CtrlInterval)};

        for (auto &point : path) {
          fileStream << "v" << DELIMITER_OUT;
          Point2DPolynom temp{point / this->problem.Env.ScaleFactor}; 
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        }

        vertexInd += path.size();
        vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
      }

      for (auto &tup : vertexRanges) {
        auto [start, end] = tup;
        end -= 1;
        for (int i = start; i < end; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << "\n";
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<Point2DPolynom> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            fileStream << node.Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.Closest->Position / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
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
void Solver<Point3DPolynom, false>::savePaths(const FileStruct file) {
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
    if (file.type == Obj) {
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;

      fileStream << "o Paths\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<Point3DPolynom> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<Point3DPolynom> &plan{holder.Plan};
          unsigned startingInd{vertexInd};
          for (int k{0}; k < plan.size() - 1; ++k) {
            Point3DPolynom actPoint{plan[k]};
            Point3DPolynom nextPoint{plan[k + 1]};

            std::deque<Point3DPolynom> trajectory{actPoint.SampleTrajectory(nextPoint, this->problem.CtrlInterval)};
            for (auto &p : trajectory) {
              Point3DPolynom temp{p / this->problem.Env.ScaleFactor};
              fileStream << "v" << DELIMITER_OUT;
              temp.PrintPosition(fileStream);
              fileStream << "\n";
            }
            vertexInd += trajectory.size();
            vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
          } 
        }
      }
      
      for (auto &tup : vertexRanges) {
        auto [start, end] = tup;
        end -= 1;
        for (int i = start; i < end; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << "\n";
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<Point3DPolynom> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<Point3DPolynom> &plan{holder.Plan};
          for (int k{0}; k < plan.size() - 1; ++k) {
            if (plan[k] == plan[k + 1]) {
              continue;
            }
            fileStream << plan[k] / this->problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / this->problem.Env.ScaleFactor << "\n";
          }
          fileStream << "\n";
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
void Solver<Point2DPolynom, false>::savePaths(const FileStruct file) {
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
    if (file.type == Obj) {
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;

      fileStream << "o Paths\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<Point2DPolynom> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<Point2DPolynom> &plan{holder.Plan};
          unsigned startingInd{vertexInd};
          for (int k{0}; k < plan.size() - 1; ++k) {
            Point2DPolynom actPoint{plan[k]};
            Point2DPolynom nextPoint{plan[k + 1]};

            std::deque<Point2DPolynom> trajectory{actPoint.SampleTrajectory(nextPoint, this->problem.CtrlInterval)};
            for (auto &p : trajectory) {
              Point2DPolynom temp{p / this->problem.Env.ScaleFactor};
              fileStream << "v" << DELIMITER_OUT;
              temp.PrintPosition(fileStream);
              fileStream << "\n";
            }
            vertexInd += trajectory.size();
            vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
          } 
        }
      }
      
      for (auto &tup : vertexRanges) {
        auto [start, end] = tup;
        end -= 1;
        for (int i = start; i < end; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << "\n";
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<Point2DPolynom> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<Point2DPolynom> &plan{holder.Plan};
          for (int k{0}; k < plan.size() - 1; ++k) {
            if (plan[k] == plan[k + 1]) {
              continue;
            }
            fileStream << plan[k] / this->problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / this->problem.Env.ScaleFactor << "\n";
          }
          fileStream << "\n";
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
