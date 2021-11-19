/**
 * @file prm.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief
 * @version 1.0
 * @date 11. 11. 2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef __PRM_CLASS_H__
#define __PRM_CLASS_H__

#include "dijkstra.h"
#include "opendubins/dubins.h"
#include "solver.h"

template <class R> class ProbRoadMaps : public Solver<R> {
public:
  ProbRoadMaps(Problem<R> &problem);

  void Solve() override;

private:
  std::deque<PrmNode<R>> allPoints;

  FlannHolder<Node<R>> flannIndex;

  void getPaths() override;
  void getConnected() override;
  void saveTrees(const FileStruct file) override;
};

template <> void ProbRoadMaps<Point2DDubins>::getPaths();
template <> void ProbRoadMaps<Point2DDubins>::getConnected();
template <> void ProbRoadMaps<Point2DDubins>::saveTrees(const FileStruct file);

template <class R>
ProbRoadMaps<R>::ProbRoadMaps(Problem<R> &problem) : Solver<R>(problem) {
  for (int j{0}; j < this->problem.Roots.size(); ++j) {
    PrmNode<R> &node{
        allPoints.emplace_back(this->problem.Roots[j], nullptr, nullptr, 0)};
    this->allNodes.emplace_back((Node<R> *)&node);
  }

  if (this->problem.HasGoal) {
    PrmNode<R> &node{
        allPoints.emplace_back(this->problem.Goal, nullptr, nullptr, 0)};
    this->allNodes.emplace_back((Node<R> *)&node);
  }
}

template <class R> void ProbRoadMaps<R>::Solve() {
  StopWatch sw;

  if (SaveGoals <= this->problem.SaveOpt) {
    this->saveCities(this->problem.FileNames[SaveGoals]);
  }

  bool solved{true};
  sw.Start();
  // sample random points
  for (int i{0}; i < this->problem.MaxIterations; ++i) {
    R rndPoint;

    this->rnd.RandomPointInSpace(rndPoint);
    if (this->problem.Env.Collide(rndPoint)) { // check limits and collisions
      continue;
    }

    allPoints.emplace_back(rndPoint, nullptr, nullptr, i + 1);
  }

  // calculate edges
  // create queries for knn/add points to FLANN KD-tree
  size_t newPointsCount{allPoints.size()};
  flann::Matrix<float> newPoints{new float[newPointsCount * PROBLEM_DIMENSION],
                                 newPointsCount, PROBLEM_DIMENSION};
  for (int j{0}; j < newPointsCount; ++j) {
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      newPoints[j][i] = allPoints[j].Position[i];
    }
  }
  flannIndex.CreateIndex(newPoints);

  // perform knn
  std::vector<std::vector<int>> indices;
  std::vector<std::vector<float>> dists;
  // (euler * (1 + 1 / config_dimension) + 1) * log(num_points) -> k_PRM_Class*
  // = (euler * (1 + 1 / config_dimension) + 1)
  double nnD{(M_E * (1 + 1 / 3.0) + 1) * log10(allPoints.size())};
  int nn{static_cast<int>(nnD)};
  flannIndex.Index->knnSearch(newPoints, indices, dists, nn,
                              flann::SearchParams(FLANN_NUM_SEARCHES));

  for (int i{0}; i < allPoints.size(); ++i) {
    PrmNode<R> &act{allPoints[i]};
    std::vector<int> &indRow{indices[i]};
    for (int neighInd{0}; neighInd < nn; ++neighInd) {
      int neighID{indRow[neighInd]};
      PrmNode<R> &neighbour{allPoints[neighID]};

      bool isPathBetweenFree{
          this->isPathFree(act.Position, neighbour.Position)};
      if (isPathBetweenFree) {
        double distBetween{act.Position.Distance(neighbour.Position)};
        act.VisibleNodes.emplace(&neighbour, distBetween);
      }
    }
  }

  getPaths();

  sw.Stop();

  getConnected();
  this->getAllPaths();

  for (int i{0}; i < this->connected.size(); ++i) {
    if (!this->connected[i]) {
      solved = false;
      break;
    }
  }

  if (this->problem.ComputeTSP) {
    this->computeTsp();
  }

  if (SaveTree <= this->problem.SaveOpt) {
    this->saveTrees(this->problem.FileNames[SaveTree]);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams],
                     this->problem.MaxIterations, solved, sw.GetElapsed());
  }

  if (SaveRoadmap <= this->problem.SaveOpt) {
    this->savePaths(this->problem.FileNames[SaveRoadmap]);
  }

  if (SaveTSPFile <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSPFile]);
  }

  if (this->problem.ComputeTSP && SaveTSPPaths <= this->problem.SaveOpt) {
    this->saveTspPaths(this->problem.FileNames[SaveTSPPaths]);
  }
}

template <class R> void ProbRoadMaps<R>::getPaths() {
  for (int i{1}; i < this->problem.GetNumRoots(); ++i) {
    Dijkstra<R> dijkstra;
    std::vector<int> goals;
    for (int j{0}; j < i; ++j) {
      goals.push_back(j);
    }
    std::deque<DistanceHolder<R>> plans{dijkstra.findPath(i, goals, allPoints)};
    for (auto &holder : plans) {
      int id1{holder.Node1->ID};
      int id2{holder.Node2->ID};
      this->neighboringMatrix(id1, id2) = holder;
    }
  }
}

template <class R> void ProbRoadMaps<R>::getConnected() {
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{i}; j < this->problem.GetNumRoots(); ++j) {
      if (!this->neighboringMatrix.Exists(i, j) ||
          this->neighboringMatrix(i, j).Distance ==
              std::numeric_limits<double>::max()) {
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

template <class R> void ProbRoadMaps<R>::saveTrees(const FileStruct file) {
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
      fileStream << "o Trees\n";
      for (int i{0}; i < this->allPoints.size(); ++i) {
        R temp{this->allPoints[i].Position / this->problem.Env.ScaleFactor};
        fileStream << "v" << DELIMITER_OUT;
        temp.PrintPosition(fileStream);
        fileStream << "\n";
      }

      for (auto &node : allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (node.ID < neigh->ID) {
            fileStream << "l" << DELIMITER_OUT << node.ID + 1 << DELIMITER_OUT
                       << neigh->ID + 1 << "\n";
          }
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension
                 << "\n";
      for (auto &node : this->allPoints) {
        for (auto &pair : node.VisibleNodes) {
          auto [neigh, dist] = pair;
          if (node.ID < neigh->ID) {
            fileStream << node.Position / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT
                       << neigh->Position / this->problem.Env.ScaleFactor
                       << DELIMITER_OUT << node.ID << DELIMITER_OUT
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

#endif
