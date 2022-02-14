/**
 * @file lazy.h
 * @author janosjar (janosjar@fel.cvut.cz)
 * @brief 
 * @version 0.1
 * @date 2021-08-20
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __LAZY_H__
#define __LAZY_H__

#include <deque>
#include <flann/flann.hpp>
#include "solver.h"
#include "common.h"
#include "tsp-handler.h"

template<class R>
class LazyTSPBase : public Solver<R> {
  public:
    LazyTSPBase(Problem<R> &problem) : Solver<R>(problem) {} 
    ~LazyTSPBase();

    void Solve() override;
  protected:
    int numTrees;
    std::deque<Node<R>> rootNodes;
    std::deque<Tree<R> *> treesToDel;
    
    void runRRT(DistanceHolder<R> *edge, int &iterations);
    void processResults(TSPOrder &solution, std::deque<std::tuple<int,int>> &edgePairs);
    virtual double updateEdge(std::tuple<int, int> selectedEdge, int iter) = 0;
    virtual void rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) = 0;

    void getPaths() override;
    virtual void savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) = 0;
    virtual void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) = 0;
};

template<class R, bool = isDubins<R>::value>
class LazyTSP : public LazyTSPBase<R> {  
};

template<class R>
class LazyTSP<R, false> : public LazyTSPBase<R> {
  public:
    LazyTSP(Problem<R> &problem);

  protected:
    double updateEdge(std::tuple<int, int> selectedEdge, int iter) override;
    void rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) override;

    void savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) override;
    void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) override;
};

template<class R>
class LazyTSP<R, true> : public LazyTSPBase<R> {
  public:
    LazyTSP(Problem<R> &problem);

  protected:
    double updateEdge(std::tuple<int, int> selectedEdge, int iter) override;
    void rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) override;

    void savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) override;
    void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) override;
};

template<> void LazyTSP<Point2DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths);
template<> void LazyTSP<Point3DDubins>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths);

template<class R>
LazyTSP<R, false>::LazyTSP(Problem<R> &problem) : LazyTSPBase<R>(problem) {
  // create adjacency matrix
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Node<R> &node{this->rootNodes.emplace_back(this->problem.Roots[i], nullptr, nullptr, 0, 0)};
    this->allNodes.push_back(&node);
  }
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (this->neighboringMatrix(i, j).Exists()) {
        continue;
      }
      this->neighboringMatrix(i, j) = DistanceHolder<R>(&(this->rootNodes[i]), &(this->rootNodes[j]), this->rootNodes[i].Distance(this->rootNodes[j]));
    }
  }
}

template<class R>
LazyTSP<R, true>::LazyTSP(Problem<R> &problem) : LazyTSPBase<R>(problem) {
  // create adjacency matrix
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Node<R> &node{this->rootNodes.emplace_back(this->problem.Roots[i], nullptr, nullptr, 0, 0)};
    this->allNodes.push_back(&node);
  }
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (this->neighboringMatrix(i, j, 0, 0).Exists()) {
        continue;
      }
      DistanceHolder<R> dist{&(this->rootNodes[i]), &(this->rootNodes[j]), this->rootNodes[i].Distance(this->rootNodes[j])};
      this->neighboringMatrix.AddLink(dist, i, j, 0, 0);
    }
  }
}

template <class R>
LazyTSPBase<R>::~LazyTSPBase() {
  for (auto ptr : treesToDel) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}

template <class R>
void LazyTSPBase<R>::Solve() {
  double prevDist{-1}, newDist{0};
  StopWatch watches;
  std::deque<std::tuple<int, int>> selectedEdges;

  if (SaveGoals <= this->problem.SaveOpt) {
    this->saveCities(this->problem.FileNames[SaveGoals]);
  }

  watches.Start();

  FileStruct tempTsp;
  tempTsp.fileName = TEMP_TSP;
  tempTsp.type = Map;

  bool solved{false};
  int iter{0};
  int maxNumTrees{this->problem.GetNumRoots() * (this->problem.GetNumRoots() - 1) / 2};
  while (!solved && iter != maxNumTrees * this->problem.MaxIterations) {
    prevDist = newDist;

    // run TSP = create file, execute, read output
    TSPMatrix tsp{this->problem, this->neighboringMatrix};
    TSPOrder tspSolution;
    if (this->problem.TspType == Concorde) {
      tspSolution = tsp.SolveByConcorde();
    } else if (this->problem.TspType == LKH) {
      tspSolution = tsp.SolveByLKH();
    } else {
      ERROR("TSP solver not implemented");
      exit(1);
    }

    if (tspSolution.size() == 0) {
      ERROR("Lazy TSP: temporary TSP file error");
      return;
    }

    processResults(tspSolution, selectedEdges);

    // run RRT for selected edges, recompute new distance
    newDist = 0;
    for (auto &pair : selectedEdges) {
      newDist += updateEdge(pair, iter);
    }

    solved = (newDist >= prevDist - SFF_TOLERANCE && newDist <= prevDist + SFF_TOLERANCE);
  }
  watches.Stop();

  if (SaveRoadmap <= this->problem.SaveOpt) {
    this->savePaths(this->problem.FileNames[SaveRoadmap], selectedEdges);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.FileNames[SaveParams], iter, solved, watches.GetElapsed(), selectedEdges);
  }

  if (SaveTSPFile <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSPFile]);
  }
}

template<class R>
double LazyTSP<R, false>::updateEdge(std::tuple<int, int> selectedEdge, int iter) {
  int first, second;
  std::tie(first, second) = selectedEdge;
  DistanceHolder<R> &edge{this->neighboringMatrix(first, second)};
  if (edge.Plan.empty()) {
    this->runRRT(&edge, iter);
  } 
  return edge.Distance;
}

template<class R>
double LazyTSP<R, true>::updateEdge(std::tuple<int, int> selectedEdge, int iter) {
  int first, second;
  std::tie(first, second) = selectedEdge;
  DistanceHolder<R> &edge{this->neighboringMatrix(first, second, 0, 0)};
  if (edge.Plan.empty()) {
    this->runRRT(&edge, iter);
  } 
  return edge.Distance;
}

template <class R>
void LazyTSPBase<R>::getPaths() {

}

template <class R>
void LazyTSPBase<R>::runRRT(DistanceHolder<R> *edge, int &iterations) {
  Node<R> *goal{edge->Node2};
  Tree<R> *rrtTree{new Tree<R>};
  treesToDel.push_back(rrtTree);
  
  Node<R> &start{rrtTree->Leaves.emplace_back(edge->Node1->Position, rrtTree, nullptr, 0, 0)};
  rrtTree->Root = &start;
  this->allNodes.push_back(&start);
  flann::Matrix<float> rootMat{new float[1 * PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
  for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
    rootMat[0][i] = start.Position[i];
  }
  rrtTree->Flann.CreateIndex(rootMat);

  bool solved{false};
  int iter{0};
  while (iter < this->problem.MaxIterations && !solved) {
    ++iter;
    R rndPoint, newPoint;
    Node<R> *newNode;
    this->rnd.RandomPointInSpace(rndPoint);  // NO PRIORITY BIAS!!!

    std::vector<std::vector<int>> indices;
    std::vector<std::vector<float>> dists;
    flann::Matrix<float> rndPointMat{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION}; // just one point to add
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      rndPointMat[0][i] = rndPoint[i];
    }
    rrtTree->Flann.Index->knnSearch(rndPointMat, indices, dists, 1, flann::SearchParams(FLANN_NUM_SEARCHES));
    delete[] rndPointMat.ptr();
    
    Node<R> &neighbor{rrtTree->Leaves[indices[0][0]]};
    Node<R> *nearest{&neighbor};
    
    // get point in this direction, check for collisions
    newPoint = neighbor.Position.GetStateInDistance(rndPoint, this->problem.SamplingDist);
    if (this->problem.Env.Collide(newPoint) || !this->isPathFree(neighbor.Position, newPoint)) {
      continue;    
    }

    // rrt star
    if (this->problem.Optimize) {
      double bestDist{nearest->Distance(newPoint) + nearest->DistanceToRoot()};
      double krrt{2 * M_E * log10(rrtTree->Leaves.size() + 1)};
      indices.clear();
      dists.clear();

      flann::Matrix<float> newPointMat{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
      for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
        newPointMat[0][i] = newPoint[i];
      }
      rrtTree->Flann.Index->knnSearch(newPointMat, indices, dists, krrt, flann::SearchParams(FLANN_NUM_SEARCHES));

      std::vector<int> &indRow{indices[0]};
      for (int &ind : indRow) {
        Node<R> &neighbor{rrtTree->Leaves[ind]};
        double neighDist{neighbor.Distance(newPoint) + neighbor.DistanceToRoot()};
        if (neighDist < bestDist - SFF_TOLERANCE && this->isPathFree(neighbor.Position, newPoint)) {
          bestDist = neighDist;
          nearest = &neighbor;
        }
      }

      newNode = &(rrtTree->Leaves.emplace_back(newPoint, rrtTree, nearest, nearest->Distance(newPoint), iter));
      nearest->Children.push_back(newNode);

      for (int &ind : indRow) {
        Node<R> &neighbor{rrtTree->Leaves[ind]}; // offset goal node
        double newPointDist{neighbor.Distance(newPoint)};
        double proposedDist{bestDist + newPointDist};
        if (proposedDist < neighbor.DistanceToRoot() - SFF_TOLERANCE && this->isPathFree(newPoint, neighbor.Position)) {
          // rewire
          rewireNodes(newNode, neighbor, newPointDist);
        }
      }
    } else {
      newNode = &(rrtTree->Leaves.emplace_back(newPoint, rrtTree, nearest, this->problem.SamplingDist, iter));
      nearest->Children.push_back(newNode);
    }
    this->allNodes.push_back(newNode);

    // add to flann
    flann::Matrix<float> pointToAdd{new float[PROBLEM_DIMENSION], 1, PROBLEM_DIMENSION};
    for (int i{0}; i < PROBLEM_DIMENSION; ++i) {
      pointToAdd[0][i] = newPoint[i];
    }
    rrtTree->Flann.Index->addPoints(pointToAdd);
    rrtTree->Flann.PtrsToDel.push_back(pointToAdd.ptr());

    // check goal
    double goalDistance{goal->Distance(*newNode)};
    if (goalDistance < this->problem.DistTree && this->isPathFree(newNode->Position, goal->Position)) {
      solved = true;
      edge->Distance = goalDistance + newNode->DistanceToRoot();
      
      // fill-in plan
      edge->Plan.push_front(goal->Position);
      Node<R> *nodeToPush{newNode};
      edge->Plan.push_front(nodeToPush->Position);
      while (!nodeToPush->IsRoot()) {
        nodeToPush = nodeToPush->Closest;
        edge->Plan.push_front(nodeToPush->Position);
      }
    }
  }

  if (!solved) {
    edge->Distance = std::numeric_limits<double>::max();
  }

  iterations += iter;
}

template<class R>
void LazyTSP<R, false>::rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) {
  std::deque<Node<R> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring LazyTSP: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest = newDistance;
  newNode->Children.push_back(&neighbor);
}

template<class R> 
void LazyTSP<R, true>::rewireNodes(Node<R> *newNode, Node<R> &neighbor, double newDistance) {
  std::deque<Node<R> *> &children{neighbor.Closest->Children};
  auto iter{std::find(children.begin(), children.end(), &neighbor)};
  if (iter == children.end()) {
    ERROR("Fatal error when rewiring LazyTSP: Node not in children");
    exit(1);
  }
  neighbor.Closest->Children.erase(iter);
  neighbor.Closest = newNode;
  neighbor.DistanceToClosest[0] = newDistance;
  newNode->Children.push_back(&neighbor);
}

template <class R>
void LazyTSPBase<R>::processResults(TSPOrder &solution, std::deque<std::tuple<int,int>> &edgePairs) {
  edgePairs.clear();

  int prevPoint{solution.back()}; // last point of the solution
  for (int i{0}; i < solution.size(); ++i) {
    int actPoint{solution[i]};
    edgePairs.push_back(std::tuple<int, int>(prevPoint, actPoint));
    prevPoint = actPoint;
  }
}

template <class R>
void LazyTSP<R, false>::savePaths(const FileStruct file, const std::deque<std::tuple<int,int>> &selectedPaths) {
  INFO("Saving paths");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;
    ERROR(message.str());
    return;
  }

  if (fileStream.is_open()) {
    int numRoots{this->problem.GetNumRoots()};
    if (file.type == Obj) {
      // fileStream << "o Paths\n";
      // for (int i{0}; i < this->allNodes.size(); ++i) {
      //   R temp{this->allNodes[i]->Position / this->problem.Env.ScaleFactor};
      //   fileStream << "v" << DELIMITER_OUT;
      //   temp.PrintPosition(fileStream);
      //   fileStream << "\n";
      // }
      
      // for (auto &pair : selectedPaths) {
      //   int first, second;
      //   std::tie(first, second) = pair;
      //   DistanceHolder<R> &holder{this->neighboringMatrix(first, second)};

      //   std::deque<R> &plan{holder.Plan};
      //   for (int k{0}; k < plan.size() - 1; ++k) {
      //     fileStream << "l" << DELIMITER_OUT << plan[k]->ID + 1 << DELIMITER_OUT << plan[k+1]->ID + 1 << "\n";
      //   }
      // }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (auto &pair : selectedPaths) {
        int first, second;
        std::tie(first, second) = pair;
        DistanceHolder<R> &holder{this->neighboringMatrix(first, second)};

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          fileStream << plan[k] / this->problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / this->problem.Env.ScaleFactor << "\n";
        }
        fileStream << "\n";
    }
    } else {
      throw std::string("Not implemented");
    }

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    ERROR(message.str());
  }  
}

template<class R>
void LazyTSP<R, false>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) {
  INFO("Saving parameters");
  std::ofstream fileStream{file.fileName.c_str(), std::ios_base::openmode::_S_app};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;
    ERROR(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << this->problem.ID << CSV_DELIMITER;
    fileStream << this->problem.Repetition << CSV_DELIMITER;
    fileStream << iterations << CSV_DELIMITER;
    fileStream << (solved ? "solved" : "unsolved") << CSV_DELIMITER;
    fileStream << "[";
    int iter{0};
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << first;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
    iter = 0;
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << this->neighboringMatrix(first,second).Distance / this->problem.Env.ScaleFactor;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER;
    fileStream << elapsedTime.count() << "\n";
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    ERROR(message.str());
  }
}

template<class R>
void LazyTSP<R, true>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime, const std::deque<std::tuple<int,int>> &selectedEdges) {
  INFO("Saving parameters");
  std::ofstream fileStream{file.fileName.c_str(), std::ios_base::openmode::_S_app};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;
    ERROR(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << this->problem.ID << CSV_DELIMITER;
    fileStream << this->problem.Repetition << CSV_DELIMITER;
    fileStream << iterations << CSV_DELIMITER;
    fileStream << (solved ? "solved" : "unsolved") << CSV_DELIMITER;
    fileStream << "[";
    int iter{0};
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << first;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
    iter = 0;
    for (auto &pair : selectedEdges) {
      ++iter;
      int first, second;
      std::tie(first, second) = pair;
      fileStream << this->neighboringMatrix(first,second, 0, 0).Distance / this->problem.Env.ScaleFactor;
      if (iter != this->problem.GetNumRoots()) {
        fileStream << CSV_DELIMITER_2;
      }
    }
    fileStream << "]" << CSV_DELIMITER;
    fileStream << elapsedTime.count() << "\n";
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    ERROR(message.str());
  }
}

#endif 
