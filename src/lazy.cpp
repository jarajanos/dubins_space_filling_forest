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

template<>
void LazyTSP<Point2DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) {
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
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(first, second, 0, 0)};
        std::deque<Point2DDubins> &plan{holder.Plan};
        unsigned startingInd{vertexInd};
        for (int m{0}; m < plan.size(); ++m) {
          Point2DDubins actPoint{plan[m]};

          fileStream << "v" << DELIMITER_OUT;
          Point2DDubins temp{actPoint / problem.Env.ScaleFactor}; 
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
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<Point2DDubins> &holder{this->neighboringMatrix(first, second, 0, 0)};

        std::deque<Point2DDubins> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          opendubins::State finishDub{plan[k][0], plan[k][1], plan[k].GetHeading()};
          opendubins::State startDub{plan[k+1][0], plan[k+1][1], plan[k+1].GetHeading()};
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

template<>
void LazyTSP<Point3DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) {
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
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(first, second, 0, 0)};
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
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<Point3DDubins> &holder{this->neighboringMatrix(first, second, 0, 0)};
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
