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
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;
      fileStream << "o Trees\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 

          unsigned startingInd{vertexInd};
          auto path{node.Position.SampleDubinsPathTo(neigh->Position, this->problem.CollisionDist)};
          for (int i{0}; i < path.size(); ++i) {
            fileStream << "v" << DELIMITER_OUT;
            Point2DDubins temp{path[i] / problem.Env.ScaleFactor}; 
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
          vertexInd += path.size();

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));    
        }
      }
      
      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 
          opendubins::State startDub{node.Position[0], node.Position[1], node.Position.GetHeading()};
          opendubins::State finishDub{neigh->Position[0], neigh->Position[1], neigh->Position.GetHeading()};
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

template <> 
void ProbRoadMaps<Point3DDubins>::saveTrees(const FileStruct file) {
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
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 

          unsigned startingInd{vertexInd};
          auto path{node.Position.SampleDubinsPathTo(neigh->Position, this->problem.CollisionDist)};
          for (int i{0}; i < path.size(); ++i) {
            fileStream << "v" << DELIMITER_OUT;
            Point3DDubins temp{path[i] / problem.Env.ScaleFactor}; 
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
          vertexInd += path.size();

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));    
        }
      }
      
      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 
          opendubins::State3D finishDub{node.Position[0], node.Position[1], node.Position[2], node.Position.GetHeading(), node.Position.GetPitch()};
          opendubins::State3D startDub{neigh->Position[0], neigh->Position[1], neigh->Position[2], neigh->Position.GetHeading(), neigh->Position.GetPitch()};
          opendubins::Dubins3D pathFromClosest{startDub, finishDub, this->problem.DubinsRadius, this->problem.Env.Limits.mins[3], this->problem.Env.Limits.maxs[3]};
          
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
void ProbRoadMaps<Point2DPolynom, false>::saveTrees(const FileStruct file) {
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
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [ neigh, dist ] = pair;
          if (neigh == &node) {
            continue;
          } 

          unsigned startingInd{vertexInd};
          auto path{node.Position.SampleTrajectory(neigh->Position, this->problem.CtrlInterval)};
          for (auto &p : path) {
            Point2DPolynom temp{p / this->problem.Env.ScaleFactor};
            fileStream << "v" << DELIMITER_OUT;
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
          vertexInd += path.size();

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
        }
      }

      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 

          auto path{node.Position.SampleTrajectory(neigh->Position, this->problem.CtrlInterval)};
          for (int i{0}; i < path.size() - 1; ++i) {
            fileStream << path[i] / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT
                       << path[i+1] / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT << 0 << DELIMITER_OUT
                       << node.GetAge() << "\n";
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
void ProbRoadMaps<Point3DPolynom, false>::saveTrees(const FileStruct file) {
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
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [ neigh, dist ] = pair;
          if (neigh == &node) {
            continue;
          } 

          unsigned startingInd{vertexInd};
          auto path{node.Position.SampleTrajectory(neigh->Position, this->problem.CtrlInterval)};
          for (auto &p : path) {
            Point3DPolynom temp{p / this->problem.Env.ScaleFactor};
            fileStream << "v" << DELIMITER_OUT;
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
          vertexInd += path.size();

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
        }
      }

      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (neigh == &node) {
            continue;
          } 

          auto path{node.Position.SampleTrajectory(neigh->Position, this->problem.CtrlInterval)};
          for (int i{0}; i < path.size() - 1; ++i) {
            fileStream << path[i] / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT
                       << path[i+1] / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT << 0 << DELIMITER_OUT
                       << node.GetAge() << "\n";
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
