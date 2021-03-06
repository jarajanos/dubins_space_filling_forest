/**
 * @file lazy-sff.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 1.0
 * @date 17. 02. 2022
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "lazy-sff.h"

template<>
void LazySpaceForest<Point2DDubins>::Solve() {

}

template<>
void LazySpaceForest<Point3DDubins>::Solve() {
  // convert solved problem to 3D
  Problem<Point3D> convProblem{this->problem};
  SpaceForest<Point3D> baseSolver{convProblem};

  if (SaveGoals <= this->problem.SaveOpt) {
    this->saveCities(this->problem.FileNames[SaveGoals]);
  }

  int iter{0};
  bool solved{false};
  std::chrono::duration<double> elapsedTime;
  baseSolver.PureSolve(iter, solved, elapsedTime);

  if (solved) {
    // multi-goal
    if (!this->problem.HasGoal) {
      auto normalPaths{baseSolver.GetNeighboringMatrix()};
      TSPOrder baseTsp{baseSolver.GetTSPSolution()};
      int numRoots{this->problem.GetNumRoots()};
      int sumPoints{0};

      for (int i{0}; i < numRoots; ++i) {
        int actNode{baseTsp[i]};
        int nextNode{baseTsp[(i + 1) % numRoots]};

        auto link{normalPaths(actNode, nextNode)};
        sumPoints += link.Plan.size();
      }

      // distance matrix of size sum of all nodes in TSP path - numRoots (duplicates) + 1 (intensional duplicate of first node)
      std::vector<std::vector<std::vector<double>>> dtpMatrix;
      int numNodes{sumPoints - numRoots};
      dtpMatrix.resize(numNodes);
      for (int i{0}; i < numNodes; ++i) {
        dtpMatrix[i].resize(this->problem.DubinsResolution);
        for (int j{0}; j < this->problem.DubinsResolution; ++j) {
          dtpMatrix[i][j].resize(this->problem.DubinsResolution, -1);
        }
      }

      // fill in matrix with dubins distances
      int idx{0};
      Point3DDubins start;
      for (int i{0}; i < numRoots; ++i) {
        int actNode{baseTsp[i]};
        int nextNode{baseTsp[(i + 1) % numRoots]};

        auto link{normalPaths(actNode, nextNode)};
        if (i == 0) {
          start = link.Plan[0];
          pathPoints.push_back(start);
        }
        for (int l{1}; l < link.Plan.size(); ++l) {
          Point3DDubins finish{link.Plan[l]};
          pathPoints.push_back(finish);
          for (int j{0}; j < this->problem.DubinsResolution; ++j) {
            for (int k{0}; k < this->problem.DubinsResolution; ++k) {
              Point3DDubins tempStart{start}, tempFinish{finish};
              tempStart.SetHeading(j, this->problem.DubinsResolution);
              tempFinish.SetHeading(k, this->problem.DubinsResolution);
              if (this->isPathFree(tempStart, tempFinish)) {
                dtpMatrix[idx][j][k] = tempStart.Distance(tempFinish);
              }
            }
          }
          start = finish;
          idx += 1;
        }
      }

      // dijkstra
      std::vector<std::vector<DtpNode>> index;
      index.resize(numNodes + 1);
      for (int i{0}; i < numNodes + 1; ++i) {
        index[i].resize(this->problem.DubinsResolution);
        for (int j{0}; j < this->problem.DubinsResolution; ++j) {
          index[i][j].Position = i;
          index[i][j].ActAngle = j;
        } 
      }

      std::deque<DtpNode> initialNull;
      Heap<DtpNode> heap{initialNull, nullptr, false};
      
      // initialize
      for (int i{0}; i < this->problem.DubinsResolution; ++i) {
        index[0][i].Distance = 0;
        index[0][i].FirstAngle = i;
        index[0][i].LastAngle = i;
        heap.Push(&(index[0][i]));
      }

      DtpNode *finalNode{nullptr};
      while (!heap.Empty()) {
        DtpNode *actNode{heap.Get(0)};
        heap.Pop();

        if (actNode->Position == numNodes) {
          if (actNode->ActAngle == actNode->FirstAngle) {
            // found path
            finalNode = actNode;
            break;
          }

          // invalid connection
          continue;
        }

        int nextPos{actNode->Position + 1};
        for (int i{0}; i < this->problem.DubinsResolution; ++i) {
          double dist{dtpMatrix[actNode->Position][actNode->ActAngle][i]};
          if (index[nextPos][i].FirstAngle == -1 && dist != -1) {
            index[nextPos][i].FirstAngle = actNode->FirstAngle;
            index[nextPos][i].LastAngle = actNode->ActAngle;
            index[nextPos][i].Distance = actNode->Distance + dist;
            heap.Push(&(index[nextPos][i]));
          }
        }
      }

      if (finalNode == nullptr) {
        solved = false;
      } else {
        DtpNode *node{finalNode};
        this->finalDistance = node->Distance;
        int order{node->Position};
        this->gatspSolution.push_front(std::tuple<int, int>(order--, node->LastAngle));
        for (; order >= 0; --order) {
          node = &(index[order][node->LastAngle]);
          this->gatspSolution.push_front(std::tuple<int, int>(order, node->LastAngle));
        }
      }
    // single goal
    } else {
      DistanceHolder<Point3D> &solution{baseSolver.GetNeighboringMatrix()(0,1)};

      // distance matrix of size of all nodes in the path
      std::vector<std::vector<std::vector<double>>> dtpMatrix;
      int numNodes{solution.Plan.size()};
      dtpMatrix.resize(numNodes);
      for (int i{0}; i < numNodes; ++i) {
        dtpMatrix[i].resize(this->problem.DubinsResolution);
        for (int j{0}; j < this->problem.DubinsResolution; ++j) {
          dtpMatrix[i][j].resize(this->problem.DubinsResolution, -1);
        }
      }

      // fill in matrix with dubins distances
      int idx{0};
      Point3DDubins start{solution.Plan[0]};
      pathPoints.push_back(start);
      for (int l{1}; l < numNodes; ++l) {
        Point3DDubins finish{solution.Plan[l]};
        pathPoints.push_back(finish);
        for (int j{0}; j < this->problem.DubinsResolution; ++j) {
          for (int k{0}; k < this->problem.DubinsResolution; ++k) {
            Point3DDubins tempStart{finish}, tempFinish{start};
            tempStart.SetHeading(j, this->problem.DubinsResolution);
            tempFinish.SetHeading(k, this->problem.DubinsResolution);
            if (this->isPathFree(tempStart, tempFinish)) {
              dtpMatrix[idx][j][k] = tempStart.Distance(tempFinish);
            }
          }
        }
        start = finish;
        idx += 1;
      }

      // dijkstra
      std::vector<std::vector<DtpNode>> index;
      index.resize(numNodes);
      for (int i{0}; i < numNodes; ++i) {
        index[i].resize(this->problem.DubinsResolution);
        for (int j{0}; j < this->problem.DubinsResolution; ++j) {
          index[i][j].Position = i;
          index[i][j].ActAngle = j;
        } 
      }

      std::deque<DtpNode> initialNull;
      Heap<DtpNode> heap{initialNull, nullptr, false};
      
      // initialize
      for (int i{0}; i < this->problem.DubinsResolution; ++i) {
        index[0][i].Distance = 0;
        index[0][i].FirstAngle = i;
        index[0][i].LastAngle = i;
        heap.Push(&(index[0][i]));
      }

      DtpNode *finalNode{nullptr};
      while (!heap.Empty()) {
        DtpNode *actNode{heap.Get(0)};
        heap.Pop();

        if (actNode->Position == numNodes - 1) {
          // found path
          finalNode = actNode;
          break;
        }

        int nextPos{actNode->Position + 1};
        for (int i{0}; i < this->problem.DubinsResolution; ++i) {
          double dist{dtpMatrix[actNode->Position][actNode->ActAngle][i]};
          if (index[nextPos][i].FirstAngle == -1 && dist != -1) {
            index[nextPos][i].FirstAngle = actNode->FirstAngle;
            index[nextPos][i].LastAngle = actNode->ActAngle;
            index[nextPos][i].Distance = actNode->Distance + dist;
            heap.Push(&(index[nextPos][i]));
          }
        }
      }

      if (finalNode == nullptr) {
        solved = false;
      } else {
        DtpNode *node{finalNode};
        this->finalDistance = node->Distance;
        int order{node->Position};
        this->gatspSolution.push_front(std::tuple<int, int>(order--, node->LastAngle));
        for (; order >= 0; --order) {
          node = &(index[order][node->LastAngle]);
          this->gatspSolution.push_front(std::tuple<int, int>(order, node->LastAngle));
        }
      }
    }
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams], iter, solved, elapsedTime);
  }

  if (SaveTSPFile <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSPFile]);
  }

  if (solved && SaveTSPPaths <= this->problem.SaveOpt) {
    this->saveTspPaths(this->problem.FileNames[SaveTSPPaths]);
  }
}

template<> 
void LazySpaceForest<Point2DDubins>::saveTspPaths(const FileStruct file) {

}

template<>
void LazySpaceForest<Point3DDubins>::saveTspPaths(const FileStruct file) {
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
      fileStream << "o TspPaths\n";
      unsigned vertexInd{1};
      std::deque<std::tuple<unsigned,unsigned>> vertexRanges;
      int numPoints{static_cast<int>(pathPoints.size()) - 1};
      for (int i{0}; i < numPoints; ++i) {
        auto [ actNode, actAngle ] = this->gatspSolution[i];
        auto [ nextNode, nextAngle ] = this->gatspSolution[i + 1];

        Point3DDubins tempStart{this->pathPoints[nextNode]}, tempFinish{this->pathPoints[actNode]};
        tempStart.SetHeading(nextAngle, this->problem.DubinsResolution);
        tempFinish.SetHeading(actAngle, this->problem.DubinsResolution);
        auto pathSeg{tempStart.SampleDubinsPathTo(tempFinish, this->problem.CollisionDist)};
        unsigned startingInd{vertexInd};
        for (int m{0}; m < pathSeg.size(); ++m) {
          Point3DDubins actPoint{pathSeg[m]};

          fileStream << "v" << DELIMITER_OUT;
          Point3DDubins temp{actPoint / this->problem.Env.ScaleFactor}; 
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        } 
        vertexInd += pathSeg.size();
        vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1)); 
      }
      
      for (auto &pair : vertexRanges) {
        auto [ from, to ] = pair;
        for (unsigned i{from}; i < to; ++i) {
          fileStream << "l" << DELIMITER_OUT << i << DELIMITER_OUT << i + 1 << '\n';
        }
      }
    } else if (file.type == Map) {
      fileStream << "#TspPaths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      // for (int i{0}; i < numRoots; ++i) {
      //   auto [ actNode, actAngle ] = this->gatspSolution[i];
      //   auto [ nextNode, nextAngle ] = this->gatspSolution[(i + 1) % numRoots];
      //   if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
      //     ERROR("Invalid TSP solution");
      //     exit(1);
      //   }

      //   DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
      //   std::deque<R> &plan{holder.Plan};
      //   for (int m{0}; m < plan.size() - 1; ++m) {
      //     R actPoint{plan[m]};
      //     R lastPoint{plan[m + 1]};

      //     fileStream << actPoint / this->problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / this->problem.Env.ScaleFactor << "\n";
      //     lastPoint = actPoint;
      //   }  
      // }
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
void LazySpaceForest<Point3DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
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
  
    fileStream << "[" << CSV_NO_PATH << "]" << CSV_DELIMITER << "[";
  
    // do not print distance matrix
    if (finalDistance != std::numeric_limits<double>::max()) {
      fileStream << finalDistance / this->problem.Env.ScaleFactor << CSV_DELIMITER;
    } else {
      fileStream << CSV_NO_PATH;
    } 
    fileStream << "]" << CSV_DELIMITER;
    
    if (finalDistance != std::numeric_limits<double>::max()) {
      fileStream << finalDistance / this->problem.Env.ScaleFactor << CSV_DELIMITER;
    } else {
      fileStream << CSV_NO_PATH << CSV_DELIMITER;
    }

    fileStream << elapsedTime.count() << "\n";

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}
