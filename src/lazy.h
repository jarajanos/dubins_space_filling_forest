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
    
    void runRRT(DistanceHolder<Node<R>> *edge, int &iterations);
    void processResults(std::string &line, std::deque<std::tuple<int,int>> &edgePairs, double &pathLength);

    void getPaths() override;
    void saveTsp(const FileStruct file) override;
};

template<class R>
LazyTSP<R>::LazyTSP(Problem<R> &problem) : Solver<R>(problem) {
  // create adjacency matrix
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    Node<R> &node{rootNodes.emplace_back(this->problem.Roots[i], nullptr, nullptr, 0, 0, 0)};
    this->allNodes.push_back(&node);
  }
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    for (int j{0}; j < this->problem.GetNumRoots(); ++j) {
      if (this->neighboringMatrix(i, j).Exists()) {
        continue;
      }
      this->neighboringMatrix(i, j) = DistanceHolder<Node<R>>(&(rootNodes[i]), &(rootNodes[j]), rootNodes[i].Position.Distance(rootNodes[j].Position));
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
    std::string id{"id_" + std::to_string(this->problem.Repetition) + "_"};
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
      if (edge.Plan.empty()) {
        runRRT(&edge, iter);
      } 
      newDist += edge.Distance;
    }

    solved = (newDist >= prevDist - SFF_TOLERANCE && newDist <= prevDist + SFF_TOLERANCE);
  }
  watches.Stop();

  if (SaveRoadmap <= this->problem.SaveOpt) {
    //this->savePaths(this->problem.FileNames[SaveRoadmap], selectedEdges);
  }

  if (SaveParams <= this->problem.SaveOpt) {
    //this->saveParams(this->problem.FileNames[SaveParams], iter, solved, watches.GetElapsed(), selectedEdges);
  }

  if (SaveTSP <= this->problem.SaveOpt) {
    this->saveTsp(this->problem.FileNames[SaveTSP]);
  }
}

template <class R>
void LazyTSP<R>::getPaths() {

}

template <class R>
void LazyTSP<R>::runRRT(DistanceHolder<Node<R>> *edge, int &iterations) {
  Node<R> *goal{edge->Node2};
  Tree<Node<R>> *rrtTree{new Tree<Node<R>>};
  treesToDel.push_back(rrtTree);
  
  Node<R> &start{rrtTree->Leaves.emplace_back(edge->Node1->Position, rrtTree, nullptr, 0, 0, 0)};
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
      double bestDist{newPoint.Distance(nearest->Position) + nearest->DistanceToRoot};
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
        double neighDist{newPoint.Distance(neighbor.Position) + neighbor.DistanceToRoot};
        if (neighDist < bestDist - SFF_TOLERANCE && this->isPathFree(neighbor.Position, newPoint)) {
          bestDist = neighDist;
          nearest = &neighbor;
        }
      }

      newNode = &(rrtTree->Leaves.emplace_back(newPoint, rrtTree, nearest, nearest->Position.Distance(newPoint), bestDist, iter));
      nearest->Children.push_back(newNode);

      for (int &ind : indRow) {
        Node<R> &neighbor{rrtTree->Leaves[ind]}; // offset goal node
        double newPointDist{neighbor.Position.Distance(newPoint)};
        double proposedDist{bestDist + newPointDist};
        if (proposedDist < neighbor.DistanceToRoot - SFF_TOLERANCE && this->isPathFree(newPoint, neighbor.Position)) {
          // rewire
          std::deque<Node<R> *> &children{neighbor.Closest->Children};
          auto iter{find(children.begin(), children.end(), &neighbor)};
          if (iter == children.end()) {
            ERROR("Fatal error: Node not in children");
            exit(1);
          }
          neighbor.Closest->Children.erase(iter);
          neighbor.Closest = newNode;
          neighbor.Root = newNode->Root;
          neighbor.DistanceToClosest = newPointDist;
          neighbor.DistanceToRoot = proposedDist;
          newNode->Children.push_back(&neighbor);
        }
      }
    } else {
      newNode = &(rrtTree->Leaves.emplace_back(newPoint, rrtTree, nearest, this->problem.SamplingDist, nearest->DistanceToRoot + this->problem.SamplingDist, iter));
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
    double goalDistance{goal->Position.Distance(newNode->Position)};
    if (goalDistance < this->problem.DistTree && this->isPathFree(newNode->Position, goal->Position)) {
      solved = true;
      edge->Distance = goalDistance + newNode->DistanceToRoot;
      
      // fill-in plan
      edge->Plan.push_front(goal);
      Node<R> *nodeToPush{newNode};
      edge->Plan.push_front(nodeToPush);
      while (!nodeToPush->IsRoot()) {
        nodeToPush = nodeToPush->Closest;
        edge->Plan.push_front(nodeToPush);
      }
    }
  }

  if (!solved) {
    edge->Distance = std::numeric_limits<double>::max();
  }

  iterations += iter;
}

template <class R>
void LazyTSP<R>::processResults(std::string &line, std::deque<std::tuple<int,int>> &edgePairs, double &pathLength) {
  std::string parsedPart;
  ParseString(line, parsedPart, line, resultDelimiter);
  pathLength = std::stod(parsedPart);

  ParseString(line, parsedPart, line, resultDelimiter);
  int prevPoint{std::stoi(parsedPart)};
  for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
    ParseString(line, parsedPart, line, resultDelimiter);
    int actPoint{std::stoi(parsedPart)};
    edgePairs.push_back(std::tuple<int, int>(prevPoint, actPoint));
    prevPoint = actPoint;
  }
}

template <class R>
void LazyTSP<R>::saveTsp(const FileStruct file) {
  INFO("Saving TSP file");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;
    ERROR(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << "NAME: " << this->problem.ID << "\n";
    fileStream << "COMMENT:\n";
    fileStream << "TYPE: TSP\n";
    fileStream << "DIMENSION: " << this->problem.GetNumRoots() << "\n";
    fileStream << "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    fileStream << "EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW\n";

    fileStream << "EDGE_WEIGHT_SECTION\n";
    for (int i{0}; i < this->problem.GetNumRoots(); ++i) {
      for (int j{0}; j < i; ++j) {
        fileStream << this->neighboringMatrix(i, j).Distance / this->problem.Env.ScaleFactor << TSP_DELIMITER;
      }
      fileStream << "0\n";
    }
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    ERROR(message.str());
  }
}

#endif 
