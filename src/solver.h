/**
 * @file solver.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 06. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <deque>
#include <chrono>
#include "problem.h"
#include "random-generator.h"
#include "common.h"
#include "tsp-handler.h"

#include "opendubins/dubins.h"

template<class R>
class Solver {
  public:
    Solver(Problem<R> &problem);

    virtual ~Solver() {
    }

    virtual void Solve() = 0;
  protected:
    Problem<R> &problem;
    RandomGenerator<R> rnd;

    std::deque<Node<R> *> allNodes;
    std::deque<Tree<R>> trees;
    DistanceMatrix<DistanceHolder<R>> neighboringMatrix;
    std::deque<Tree<R> *> connectedTrees;
    TSPOrder tspSolution;
    GATSPOrder gatspSolution;

    bool isPathFree(const R start, const R finish);

    void initNeighboringMatrix();

    virtual void getPaths() = 0;
    virtual void getAllPaths();
    void checkPlan();
    void checkDistances(std::deque<Node<R> *> &plan, double distanceToCheck);
    double computeDistance(std::deque<R> &plan);

    void computeTsp();

    virtual void checkIterationSaves(const int iter);
    virtual void saveTrees(const FileStruct file);
    virtual void saveCities(const FileStruct file);
    virtual void savePaths(const FileStruct file);
    virtual void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);
    virtual void saveTsp(const FileStruct file);
    virtual void saveTspPaths(const FileStruct file);
};

template <> void Solver<Point2DDubins>::initNeighboringMatrix();
template <> void Solver<Point2DDubins>::getAllPaths();

template <> void Solver<Point2DDubins>::computeTsp();
template <> void Solver<Point2DDubins>::saveTspPaths(const FileStruct file);

template <> void Solver<Point2DDubins>::saveTrees(const FileStruct file);
template <> void Solver<Point2DDubins>::saveTsp(const FileStruct file);
template <> void Solver<Point2DDubins>::savePaths(const FileStruct file);
template <> void Solver<Point2DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);

template<class R>
Solver<R>::Solver(Problem<R> &problem) : problem{problem}, rnd{problem.Env.Limits} {
  initNeighboringMatrix();
}

template<class R>
void Solver<R>::initNeighboringMatrix() {
  this->neighboringMatrix = DistanceMatrix<DistanceHolder<R>>(problem.GetNumRoots());
}

template <class R>
double Solver<R>::computeDistance(std::deque<R> &plan) {
  R previous;
  bool first{true};
  double distance{0};
  for (R pos : plan) {
    if (!first) {
      distance += previous.Distance(pos);
    } else {
      first = false;
    }
    
    previous = pos;
  }

  return distance;
}

template<class R>
void Solver<R>::getAllPaths() {
  int numRoots{this->problem.GetNumRoots()};
  for (int k{0}; k < this->connectedTrees.size(); ++k) {
    int id3{connectedTrees[k]->Root->ID};
    for (int i{0}; i < this->connectedTrees.size(); ++i) {
      int id1{connectedTrees[i]->Root->ID};
      if (i == k || !this->neighboringMatrix.Exists(id1, id3)) {
        continue;
      }
      DistanceHolder<R> &holder1{this->neighboringMatrix(id1, id3)};
      for (int j{0}; j < this->connectedTrees.size(); ++j) {
        int id2{connectedTrees[j]->Root->ID};
        if (i == j || !this->neighboringMatrix.Exists(id2, id3)) {
          continue;
        }

        DistanceHolder<R> &link{this->neighboringMatrix(id1, id2)};
        DistanceHolder<R> &holder2{this->neighboringMatrix(id3, id2)};

        Node<R> *node1;
        Node<R> *node2;
        std::deque<R> plan1;
        std::deque<R> plan2;
        std::deque<R> finalPlan;

        bool reversed1{false};
        bool reversed2{false};
        plan1 = holder1.Plan;
        if (holder1.Node1->SourceTree->Root->ID == id1) {
          node1 = holder1.Node1;
        } else {
          if (problem.Dimension == D2Dubins) {
            ERROR("Error in neighboring matrix: DistanceHolder's plan can't be reversed in Dubins scenario - problem is not symmetric!");
          }
          node1 = holder1.Node2;
          std::reverse(plan1.begin(), plan1.end());
          reversed1 = true;
        }

        plan2 = holder2.Plan;
        if (holder2.Node1->SourceTree->Root->ID == id2) {
          node2 = holder2.Node1;
        } else {
          if (problem.Dimension == D2Dubins) {
            ERROR("Error in neighboring matrix: DistanceHolder's plan can't be reversed in Dubins scenario - problem is not symmetric!");
          }
          node2 = holder2.Node2;
          std::reverse(plan2.begin(), plan2.end());
          reversed2 = true;
        }

        R last;
        while (!plan1.empty() && !plan2.empty() && plan1.back() == plan2.back()) {
          last = plan1.back();
          plan1.pop_back();
          plan2.pop_back();

        }

        while (plan1.size() > 0) {
          finalPlan.push_back(plan1.front());
          plan1.pop_front();
        }
        finalPlan.push_back(last);
        while (plan2.size() > 0) {
          finalPlan.push_back(plan2.back());
          plan2.pop_back();
        }

        double distance{computeDistance(finalPlan)};
        if (distance < link.Distance - SFF_TOLERANCE) {     // for nonexisting connections always true (infinite distance at init)
          this->neighboringMatrix(id1, id2) = DistanceHolder<R>(node1, node2, distance, finalPlan);
        }
      }
    }
  }
}

template<class R>
void Solver<R>::computeTsp() {
  INFO("Computing TSP");
  TSPMatrix<R> tsp{this->problem, this->neighboringMatrix};

  // already TSP, run desired solver
  if (this->problem.TspType == Concorde) {
    tspSolution = tsp.SolveByConcorde();
  } else if (this->problem.TspType == LKH) {
    tspSolution = tsp.SolveByLKH();
  } else {
    ERROR("TSP solver not implemented");
    exit(1);
  }

  if (tspSolution.size() == 0) {
    WARN("TSP not solved");
    if (SaveTSPPaths <= this->problem.SaveOpt) {
      this->problem.SaveOpt = this->problem.SaveOpt - SaveTSPPaths;
    }
  }
}

template<class R>
void Solver<R>::checkIterationSaves(const int iter) {
  if (this->problem.SaveFreq[SaveGoals] != 0 && iter % this->problem.SaveFreq[SaveGoals] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->saveCities(PrefixFileName(this->problem.FileNames[SaveGoals], prefix));
  }
  
  if (this->problem.SaveFreq[SaveTree] != 0 && iter % this->problem.SaveFreq[SaveTree] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->saveTrees(PrefixFileName(this->problem.FileNames[SaveTree], prefix));
  }

  if (this->problem.SaveFreq[SaveRoadmap] != 0 && iter % this->problem.SaveFreq[SaveRoadmap] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->savePaths(PrefixFileName(this->problem.FileNames[SaveRoadmap], prefix));
  }

  if (this->problem.SaveFreq[SaveParams] != 0 && iter % this->problem.SaveFreq[SaveParams] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->saveParams(PrefixFileName(this->problem.FileNames[SaveParams], prefix), iter, false, std::chrono::duration<double>());
  }

  if (this->problem.SaveFreq[SaveTSPFile] != 0 && iter % this->problem.SaveFreq[SaveTSPFile] == 0) {
    std::string prefix{"iter_" + std::to_string(iter) + "_"};
    this->saveTsp(PrefixFileName(this->problem.FileNames[SaveParams], prefix));
  }
}

template<class R>
void Solver<R>::saveCities(const FileStruct file) {
  INFO("Saving points");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    if (file.type == Obj) {
      fileStream << "o Points\n";;
      for (Node<R>* &node : this->allNodes){
        fileStream << "v" << DELIMITER_OUT << node->Position / problem.Env.ScaleFactor << "\n";
      }
    } else if (file.type == Map) {
      fileStream << "#Cities" << DELIMITER_OUT << problem.Dimension << "\n";
      for(Node<R>* &node : this->allNodes) {
        fileStream << node->Position / problem.Env.ScaleFactor << "\n";
      }
    } else {
      throw std::string("Unimplemented file type!");
    }
    
    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}

template<class R>
void Solver<R>::saveTrees(const FileStruct file) {
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
      for (int i{0}; i < this->allNodes.size(); ++i) {
        R temp{this->allNodes[i]->Position / problem.Env.ScaleFactor};
        fileStream << "v" << DELIMITER_OUT;
        temp.PrintPosition(fileStream);
        fileStream << "\n";
      }

      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<R> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            fileStream << "l" << DELIMITER_OUT << node.ID + 1 << DELIMITER_OUT << node.Closest->ID + 1 << "\n";
          }
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Trees" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<R> &node : this->trees[i].Leaves) {
          if (!node.IsRoot()) {
            fileStream << node.Position / problem.Env.ScaleFactor << DELIMITER_OUT << node.Closest->Position / problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
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

template<class R>
void Solver<R>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
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
  
    // print out compressed lower-diagonal matrix for ALL nodes
    int numRoots{this->problem.GetNumRoots()};
    for (int i{0}; i < numRoots; ++i) {
      for (int j{0}; j < i; ++j) {
        double dist{neighboringMatrix(i,j).Distance};
        if (dist == std::numeric_limits<double>::max()) {
          fileStream << CSV_NO_PATH;
        } else {
          fileStream << dist / problem.Env.ScaleFactor;
        }
        if (i + 1 != numRoots || j + 1 != i) {
          fileStream << CSV_DELIMITER_2;
        }
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

template <class R>
void Solver<R>::saveTsp(const FileStruct file) {
  INFO("Saving TSP file");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    fileStream << "NAME: " << this->problem.ID << "\n";
    fileStream << "COMMENT: ";
    for (int i{0}; i < this->connectedTrees.size(); ++i) {
      fileStream << this->connectedTrees[i]->Root->ID;
      if (i + 1 != connectedTrees.size()) {
        fileStream << TSP_DELIMITER;
      }
    }
    fileStream << "\n";
    fileStream << "TYPE: TSP\n";
    fileStream << "DIMENSION: " << this->connectedTrees.size() << "\n";
    fileStream << "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    fileStream << "EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW\n";

    fileStream << "EDGE_WEIGHT_SECTION\n";
    int numRoots{(int)this->connectedTrees.size()};
    for (int i{0}; i < numRoots; ++i) {
      for (int j{0}; j < i; ++j) {
        int id1{this->connectedTrees[i]->Root->ID};
        int id2{this->connectedTrees[j]->Root->ID};
        fileStream << neighboringMatrix(id1, id2).Distance / problem.Env.ScaleFactor << TSP_DELIMITER;
      }
      fileStream << "0\n";
    }

    fileStream.flush();
    fileStream.close();
  } else {
    std::stringstream message;
    message << "Cannot open file at: " << file.fileName;
    WARN(message.str());
  }
}

template<class R>
void Solver<R>::savePaths(const FileStruct file) {
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
      fileStream << "o Paths\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<R> &plan{holder.Plan};
          for (int k{0}; k < plan.size(); ++k) {
            R temp{plan[k] / problem.Env.ScaleFactor};
            fileStream << "v" << DELIMITER_OUT;
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }
        }
      }
      
      uint32_t pointCounter{1};
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<R> &plan{holder.Plan};
          for (int k{0}; k < plan.size() - 1; ++k) {
            fileStream << "l" << DELIMITER_OUT << pointCounter << DELIMITER_OUT << ++pointCounter << "\n";
          }
          ++pointCounter;
        }
      }
    } else if (file.type == Map) {
      fileStream << "#Paths" << DELIMITER_OUT << problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<R> &plan{holder.Plan};
          for (int k{0}; k < plan.size() - 1; ++k) {
            fileStream << plan[k] / problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / problem.Env.ScaleFactor << "\n";
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

template<class R>
void Solver<R>::saveTspPaths(const FileStruct file) {
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
    if (file.type == Obj) {
      fileStream << "o TspPaths\n";
      for (int i{0}; i < numRoots; ++i) {
        int actNode{tspSolution[i]};
        int nextNode{tspSolution[(i + 1) % numRoots]};
        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode)};
        
        if (holder.Node1 == NULL) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size(); ++k) {
          R temp{plan[k] / problem.Env.ScaleFactor};
          fileStream << "v" << DELIMITER_OUT;
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        }
      }
      
      uint32_t pointCounter{1};
      for (int i{0}; i < numRoots; ++i) {
        int actNode{tspSolution[i]};
        int nextNode{tspSolution[(i + 1) % numRoots]};
        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode)};

        if (holder.Node1 == NULL) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          fileStream << "l" << DELIMITER_OUT << pointCounter << DELIMITER_OUT << ++pointCounter << "\n";
        }
        ++pointCounter;
      }
    } else if (file.type == Map) {
      fileStream << "#TspPaths" << DELIMITER_OUT << problem.Dimension << "\n";;
      for (int i{0}; i < numRoots; ++i) {
        int actNode{tspSolution[i]};
        int nextNode{tspSolution[(i + 1) % numRoots]};
        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode)};

        if (holder.Node1 == NULL) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          fileStream << plan[k] / problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / problem.Env.ScaleFactor << "\n";
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

template <class R>
void Solver<R>::checkDistances(std::deque<Node<R> *> &plan, double distanceToCheck) {
  Node<R> *previous{nullptr};
  double distance{0};
  for (Node<R> *node : plan) {
    if (previous != nullptr) {
      if (!IsPathFree(node->Position, previous->Position)) {
        ERROR("Path not feasible!");
        exit(1);
      }
      distance += previous->Position.Distance(node->Position);
    }
    previous = node;
  }

  if (distance < distanceToCheck - SFF_TOLERANCE || distance > distanceToCheck + SFF_TOLERANCE) {
    ERROR("Distances mismatch");
    exit(1);
  }
}

#endif /* __SOLVER_H__ */
