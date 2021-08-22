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

template<class R>
class LazyTSP : public Solver<R> {
  public:
    LazyTSP(Problem<R> &problem); 
    ~LazyTSP();

    void Solve() override;
  private:
    int numTrees;
    inline static std::string resultDelimiter = " , ";
    std::deque<Node<R>> rootNodes;
    std::deque<Tree<Node<R>> *> treesToDel;
    
    void runRRT(DistanceHolder<T, Node<T, R>> *edge, int &iterations);
    void processResults(std::string &line, std::deque<std::tuple<int,int>> &edgePairs, T &pathLength);

    void getPaths() override;
};

template<class R>
LazyTSP<R>::LazyTSP(Problem<R> &problem) : Solver<R>(problem) {
  // create adjacency matrix
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Node<R> &node{rootNodes.emplace_back(this->problem.roots[i], nullptr, nullptr, 0, 0, 0)};
    this->allNodes.push_back(&node);
  }
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (this->neighboringMatrix(i, j).Exists()) {
        continue;
      }
      this->neighboringMatrix(i, j) = DistanceHolder<T, Node<T, R>>(&(rootNodes[i]), &(rootNodes[j]), rootNodes[i].Position.Distance(rootNodes[j].Position));
    }
  }
}

template <class R>
LazyTSP<R>::~LazyTSP() {
  for (auto ptr : treesToDel) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}

template <class R>
void LazyTSP<R>::Solve() {
  double prevDist{-1}, newDist{0};
  StopWatch watches;
  std::string resultLine;
  std::deque<std::tuple<int, int>> selectedEdges;

  watches.Start();

  FileStruct tempTsp;
  tempTsp.fileName = TEMP_TSP;
  tempTsp.type = Map;

  bool solved{false};
  int iter{0};
  int maxNumTrees{this->problem.GetNumRoots() * (this->problem.GetNumRoots() - 1) / 2};
  while (!solved && iter != maxNumTrees * this->problem.MaxIterations) {
    selectedEdges.clear();
    prevDist = newDist;

    // run TSP = create file, execute, read output
    std::string id{"id_" + std::to_string(this->problem.Iteration) + "_"};
    FileStruct runFile{PrefixFileName(tempTsp, id)};
    this->saveTsp(runFile);
    std::string command{this->problem.TspSolver};
    command.append(" --map-type=TSP_FILE --use-path-files-folder=false --use-prm=false --tsp-solver=");
    command.append(this->problem.TspType);
    command.append(" --problem=");
    command.append(runFile.fileName);
    system(command.c_str());

    std::string resultName{TEMP_RESULT};
    std::ifstream resFile{resultName.insert(0, id), std::ios::in};
    if (!resFile.good()) {
      ERROR("Lazy TSP: temporary TSP file error");
      return;
    }

    if (resFile.is_open()) {
      getline(resFile, resultLine);
    } else {
      ERROR("Lazy TSP: temporary TSP file not opened");
      return;
    }

    processResults(resultLine, selectedEdges, newDist);
    newDist = 0;

    // run RRT for selected edges, recompute new distance
    for (auto &pair : selectedEdges) {
      int first, second;
      std::tie(first, second) = pair;
      DistanceHolder<Node<R>> &edge{this->neighboringMatrix(first, second)};
      if (edge.plan.empty()) {
        runRRT(&edge, iter);
      } 
      newDist += edge.distance;
    }

    solved = (newDist >= prevDist - SFF_TOLERANCE && newDist <= prevDist + SFF_TOLERANCE);
  }
  watches.Stop();

  if (SaveRoadmap <= this->problem.SaveOpt) {
    this->savePaths(this->problem.fileNames[SaveRoadmap], selectedEdges);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    this->saveParams(this->problem.fileNames[SaveParams], iter, solved, watches.GetElapsed(), selectedEdges);
  }

  if (SaveTSP <= this->problem.saveOptions) {
    this->saveTsp(this->problem.fileNames[SaveTSP]);
  }
}

template <class R>
void LazyTSP<R>::getPaths() {

}

#endif 
