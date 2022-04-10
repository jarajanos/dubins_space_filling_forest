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
#include "opendubins/dubins3D.h"

#include "RapidTrajectoryGenerator.h"
using namespace RapidQuadrocopterTrajectoryGenerator;

template<class R>
class SolverBase {
  public:
    SolverBase(Problem<R> &problem);

    virtual ~SolverBase() {
    }

    virtual void Solve() = 0;

    DistanceMatrix<DistanceHolder<R>> &GetNeighboringMatrix();
    TSPOrder &GetTSPSolution();

  protected:
    Problem<R> &problem;
    RandomGenerator<R> rnd;

    std::deque<Node<R> *> allNodes;
    std::deque<Tree<R>> trees;
    DistanceMatrix<DistanceHolder<R>> neighboringMatrix;
    std::deque<Tree<R> *> connectedTrees;
    TSPOrder tspSolution;
    GATSPOrder gatspSolution;

    std::vector<bool> connected;

    bool isPathFree(const R start, const R finish);

    virtual void getPaths() = 0;
    virtual void getConnected();
    virtual void getAllPaths() = 0;
    void checkPlan();
    void checkDistances(std::deque<Node<R> *> &plan, double distanceToCheck);
    virtual double computeDistance(std::deque<R> &plan) = 0;

    virtual void computeTsp() = 0;

    virtual void checkIterationSaves(const int iter);
    virtual void saveCities(const FileStruct file);
    virtual void saveTrees(const FileStruct file) = 0;
    virtual void savePaths(const FileStruct file) = 0;
    virtual void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) = 0;
    virtual void saveTsp(const FileStruct file) = 0;
    virtual void saveTspPaths(const FileStruct file) = 0;
};

template<class R, bool = isDubins<R>::value>
class Solver : public SolverBase<R> {
};

template<class R>
class Solver<R, false> : public SolverBase<R> {
  public:
    Solver(Problem<R> &problem);
  
  protected:
    void getAllPaths();

    double computeDistance(std::deque<R> &plan);
    void computeTsp();

    void saveTrees(const FileStruct file);
    void savePaths(const FileStruct file);
    void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);
    void saveTsp(const FileStruct file);
    void saveTspPaths(const FileStruct file);
};

template<> void Solver<Point3DPolynom, false>::saveTrees(const FileStruct file);
template<> void Solver<Point3DPolynom, false>::savePaths(const FileStruct file);

template<class R>
class Solver<R, true> : public SolverBase<R> {
  public:
    Solver(Problem<R> &problem);
      
  protected:
    void getAllPaths();

    double computeDistance(std::deque<R> &plan);
    void computeTsp();

    void saveTrees(const FileStruct file);
    void savePaths(const FileStruct file);
    void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);
    void saveTsp(const FileStruct file);
    void saveTspPaths(const FileStruct file);
};

template<class R>
SolverBase<R>::SolverBase(Problem<R> &problem) : problem{problem}, rnd{problem} {
  this->connected = std::vector<bool>((size_t)problem.GetNumRoots(), false);
}

template<class R>
Solver<R, false>::Solver(Problem<R> &problem) : SolverBase<R>(problem) {
  this->neighboringMatrix = DistanceMatrix<DistanceHolder<R>>(problem.GetNumRoots());
}

template<class R>
Solver<R, true>::Solver(Problem<R> &problem) : SolverBase<R>(problem) {
  this->neighboringMatrix = DistanceMatrix<DistanceHolder<R>>(problem.GetNumRoots(), problem.DubinsResolution);
}


template <class R>
double Solver<R, false>::computeDistance(std::deque<R> &plan) {
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
void SolverBase<R>::getConnected() {
  for (int i{0}; i < connectedTrees.size(); ++i) {
    this->connected[this->connectedTrees[i]->Root->ID] = true;
  }
}

template<class R>
void Solver<R, false>::getAllPaths() {
  int numRoots{this->problem.GetNumRoots()};
  for (int k{0}; k < this->connectedTrees.size(); ++k) {
    int id3{this->connectedTrees[k]->Root->ID};
    for (int i{0}; i < this->connectedTrees.size(); ++i) {
      int id1{this->connectedTrees[i]->Root->ID};
      if (i == k || !this->neighboringMatrix.Exists(id1, id3)) {
        continue;
      }
      DistanceHolder<R> &holder1{this->neighboringMatrix(id1, id3)};
      for (int j{0}; j < this->connectedTrees.size(); ++j) {
        int id2{this->connectedTrees[j]->Root->ID};
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
          if (this->problem.Dimension == D2Dubins) {
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
          if (this->problem.Dimension == D2Dubins) {
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
void Solver<R, true>::getAllPaths() {
  int numRoots{this->problem.GetNumRoots()};
  int numAngles{this->problem.DubinsResolution};
  for (int k{0}; k < this->connectedTrees.size(); ++k) {
    int id3{this->connectedTrees[k]->Root->ID};
    for (int angle3{0}; angle3 < numAngles; ++angle3) {
      for (int i{0}; i < this->connectedTrees.size(); ++i) {
        int id1{this->connectedTrees[i]->Root->ID};
        if (i == k) {
          continue;
        }

        for (int angle1{0}; angle1 < numAngles; ++angle1) {
          if (!this->neighboringMatrix.Exists(id1, id3, angle1, angle3)) {
            continue;
          }

          DistanceHolder<R> &holder1{this->neighboringMatrix(id1, id3, angle1, angle3)};
          for (int j{0}; j < this->connectedTrees.size(); ++j) {
            int id2{this->connectedTrees[j]->Root->ID};
            if (i == j) {
              continue;
            }
            for (int angle2{0}; angle2 < numAngles; ++angle2) {
              if (!this->neighboringMatrix.Exists(id3, id2, angle3, angle2)) {
                continue;
              }

              DistanceHolder<R> &holder2{this->neighboringMatrix(id3, id2, angle3, angle2)};

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

              R last{plan1.back()};   // root of id3
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
                DistanceHolder<R> newLink{node1, node2, distance, finalPlan};
                this->neighboringMatrix.AddLink(newLink, id1, id2, angle1, angle2, true);
              } else {
                DistanceHolder<R> &link{this->neighboringMatrix(id1, id2, angle1, angle2)};
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

template<class R>
void Solver<R, false>::computeTsp() {
  INFO("Computing TSP");
  TSPMatrix<R> tsp{this->problem, this->neighboringMatrix};

  // already TSP, run desired solver
  if (this->problem.TspType == Concorde) {
    this->tspSolution = tsp.SolveByConcorde();
  } else if (this->problem.TspType == LKH) {
    this->tspSolution = tsp.SolveByLKH();
  } else {
    ERROR("TSP solver not implemented");
    exit(1);
  }

  if (this->tspSolution.size() == 0) {
    WARN("TSP not solved");
    if (SaveTSPPaths <= this->problem.SaveOpt) {
      this->problem.SaveOpt = this->problem.SaveOpt - SaveTSPPaths;
    }
  }
}

template <class R> 
void Solver<R, true>::computeTsp() {
  INFO("Computing TSP");
  TSPMatrix<R> gatsp{this->problem, this->neighboringMatrix};
  TSPMatrix<R> atsp{gatsp.TransformGATSPtoATSP()};

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

  this->gatspSolution = atsp.TransformATSPSolToGATSP(tempSol);

  // convert also to tsp solution
  for (auto &pair : this->gatspSolution) {
    auto [ nodeId, angleId ] = pair;
    this->tspSolution.push_back(nodeId);
  }
}

template<class R>
void SolverBase<R>::checkIterationSaves(const int iter) {
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
void SolverBase<R>::saveCities(const FileStruct file) {
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
        fileStream << "v" << DELIMITER_OUT << node->Position / this->problem.Env.ScaleFactor << "\n";
      }
    } else if (file.type == Map) {
      fileStream << "#Cities" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for(Node<R>* &node : this->allNodes) {
        fileStream << node->Position / this->problem.Env.ScaleFactor << "\n";
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
void Solver<R, false>::saveTrees(const FileStruct file) {
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
        R temp{this->allNodes[i]->Position / this->problem.Env.ScaleFactor};
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
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < this->trees.size(); ++i) {
        for (Node<R> &node : this->trees[i].Leaves) {
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

template<class R>
void Solver<R, true>::saveTrees(const FileStruct file) {
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
        Node<R> &node{*this->allNodes[i]};
        if (node.IsRoot()) {
          continue;
        }
        
        if (!node.Closest->IsRoot()) {
          unsigned startingInd{vertexInd};
          R actPoint{node.Position};
          R closestPoint{node.Closest->Position};
          auto path{closestPoint.SampleDubinsPathTo(actPoint, this->problem.CollisionDist)};

          for (auto &point : path) {
            fileStream << "v" << DELIMITER_OUT;
            R temp{point / this->problem.Env.ScaleFactor}; 
            temp.PrintPosition(fileStream);
            fileStream << "\n";
          }

          vertexInd += path.size();

          vertexRanges.push_back(std::tuple<unsigned, unsigned>(startingInd, vertexInd - 1));
        } else {
          for (auto &angle : node.GetExpandedAngles()) {
            unsigned startingInd{vertexInd};
            R actPoint{node.Position};
            R closestPoint{node.Closest->Position};
            closestPoint.SetHeading(angle, this->problem.DubinsResolution);

            auto path{closestPoint.SampleDubinsPathTo(actPoint, this->problem.CollisionDist)};
            for (auto &point : path) {
              fileStream << "v" << DELIMITER_OUT;
              R temp{point / this->problem.Env.ScaleFactor}; 
              temp.PrintPosition(fileStream);
              fileStream << "\n";
            }

            vertexInd += path.size();

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
      fileStream << "#Trees" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < this->allNodes.size(); ++i) {
        Node<R> &node{*this->allNodes[i]};
        if (node.IsRoot()) {
          continue;
        }
            
        if (!node.Closest->IsRoot()) {
          R actPoint{node.Position};
          R closestPoint{node.Closest->Position};
          auto path{closestPoint.SampleDubinsPathTo(actPoint, this->problem.CollisionDist)};

          for (int j{0}; j < path.size() - 1; ++j) {
            fileStream << path[j] / this->problem.Env.ScaleFactor << DELIMITER_OUT << path[j + 1] / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
          }
        } else {
          for (auto &angle : node.GetExpandedAngles()) {
            R actPoint{node.Position};
            R closestPoint{node.Closest->Position};
            closestPoint.SetHeading(angle, this->problem.DubinsResolution);

            auto path{closestPoint.SampleDubinsPathTo(actPoint, this->problem.CollisionDist)};
            for (int j{0}; j < path.size() - 1; ++j) {
              fileStream << path[j] / this->problem.Env.ScaleFactor << DELIMITER_OUT << path[j + 1] / this->problem.Env.ScaleFactor << DELIMITER_OUT << node.SourceTree->Root->ID << DELIMITER_OUT << node.GetAge() << "\n";
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

template<class R>
void Solver<R, false>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
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
  
    // print out compressed lower-diagonal matrix for ALL nodes
    int numRoots{this->problem.GetNumRoots()};
    for (int i{0}; i < numRoots; ++i) {
      for (int j{0}; j < i; ++j) {
        double dist{this->neighboringMatrix(i,j).Distance};
        if (dist == std::numeric_limits<double>::max()) {
          fileStream << CSV_NO_PATH;
        } else {
          fileStream << dist / this->problem.Env.ScaleFactor;
        }
        if (i + 1 != numRoots || j + 1 != i) {
          fileStream << CSV_DELIMITER_2;
        }
      }
    }
    fileStream << "]" << CSV_DELIMITER;

    double tspLength{0};
    if (this->tspSolution.size() == 0) {
      tspLength = -1;
    }
    for (int tspOrder{0}; tspOrder < numRoots && tspLength != -1; ++tspOrder) {
      int node1{this->tspSolution[tspOrder]};
      int node2{this->tspSolution[(tspOrder + 1) % numRoots]};
      if (!this->neighboringMatrix.Exists(node1, node2)) {
        tspLength = -1;
      }

      double dist{this->neighboringMatrix(node1, node2).Distance / this->problem.Env.ScaleFactor};
      if (dist != std::numeric_limits<double>::max()) {
        tspLength += dist;
      } else {
        tspLength = -1;
      }
    }
    
    if (tspLength != -1) {
      fileStream << tspLength << CSV_DELIMITER;
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

template<class R>
void Solver<R, true>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {
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
        if (this->connected[i]) {
            if (!first) {
              fileStream << CSV_DELIMITER_2;
            }
         
            fileStream << i;
            first = false;
        }
    }
    fileStream << "]" << CSV_DELIMITER << "[";
  
    // do not print distance matrix as for non-dubins case -- would be too big (except single-goal)
    if (this->problem.HasGoal) {
      double dist;
      bool exists{false};
      if (this->neighboringMatrix.Exists(0, 1, 0, 0)) {
        dist = this->neighboringMatrix(0, 1, 0, 0).Distance;
        exists = true;
      } else if (this->neighboringMatrix.Exists(1, 0, 0, 0)) {
        dist = this->neighboringMatrix(1, 0, 0, 0).Distance;
        exists = true;
      }

      if (exists && dist != std::numeric_limits<double>::max()) {
        fileStream << dist / this->problem.Env.ScaleFactor;
      } else {
        fileStream << CSV_NO_PATH;
      }
    } else {
      fileStream << CSV_NO_PATH;
    }

    fileStream << "]" << CSV_DELIMITER;

    double tspLength{0};
    int numRoots{this->problem.GetNumRoots()};
    if (this->gatspSolution.size() == 0) {
      tspLength = -1;
    }
    for (int tspOrder{0}; tspOrder < numRoots && tspLength != -1; ++tspOrder) {
      auto [ node1, angle1 ] = this->gatspSolution[tspOrder];
      auto [ node2, angle2 ] = this->gatspSolution[(tspOrder + 1) % numRoots];
      if (!this->neighboringMatrix.Exists(node1, node2, angle1, angle2)) {
        tspLength = -1;
        break;
      }

      double dist{this->neighboringMatrix(node1, node2, angle1, angle2).Distance / this->problem.Env.ScaleFactor};
      if (dist != std::numeric_limits<double>::max()) {
        tspLength += dist;
      } else {
        tspLength = -1;
      }
    }
    
    if (tspLength != -1) {
      fileStream << tspLength << CSV_DELIMITER;
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

template <class R>
void Solver<R, false>::saveTsp(const FileStruct file) {
  INFO("Saving TSP file");
  std::ofstream fileStream{file.fileName.c_str()};
  if (!fileStream.good()) {
    std::stringstream message;
    message << "Cannot create file at: " << file.fileName;

    WARN(message.str());
    return;
  }

  if (fileStream.is_open()) {
    double resolution{TSP_MAX * this->problem.Env.ScaleFactor / this->neighboringMatrix.GetMaximum() / 2};

    fileStream << "NAME: " << this->problem.ID << "\n";
    fileStream << "COMMENT: ";
    bool first{true};
    for (int i{0}; i < this->connected.size(); ++i) {
      if (!first) {
        fileStream << TSP_DELIMITER;
      }

      if (this->connected[i]) {
        fileStream << i;
        first = false;
      }
    }
    fileStream << ", resolution: " << resolution << "\n";
    fileStream << "TYPE: TSP\n";
    fileStream << "DIMENSION: " << this->connected.size() << "\n";
    fileStream << "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    fileStream << "EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW\n";

    fileStream << "EDGE_WEIGHT_SECTION\n";
    int numRoots{(int)this->connected.size()};
    for (int i{0}; i < numRoots; ++i) {
      if (!this->connected[i]) {
        continue;
      }

      for (int j{0}; j < i; ++j) {
        if (!this->connected[j]) {
          continue;
        }
        
        double dist{this->neighboringMatrix(i, j).Distance}; 
        if (dist != std::numeric_limits<double>::max()) {
          fileStream << (int)(dist / this->problem.Env.ScaleFactor * resolution) << TSP_DELIMITER;        
        } else {
          fileStream << TSP_MAX << TSP_DELIMITER;
        }
      }
      fileStream << TSP_MAX << "\n";
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
void Solver<R, true>::saveTsp(const FileStruct file) {
  ERROR("Saving DTSP file is not supported");
}

template<class R>
void Solver<R, false>::savePaths(const FileStruct file) {
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
            R temp{plan[k] / this->problem.Env.ScaleFactor};
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
      fileStream << "#Paths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{i + 1}; j < numRoots; ++j) {
          DistanceHolder<R> &holder{this->neighboringMatrix(i, j)};
          if (holder.Node1 == NULL) {
            continue;
          }

          std::deque<R> &plan{holder.Plan};
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

template<class R>
void Solver<R, true>::savePaths(const FileStruct file) {
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
              DistanceHolder<R> &holder{this->neighboringMatrix(i, j, k, l)};
              std::deque<R> &plan{holder.Plan};
              unsigned startingInd{vertexInd};
              for (int m{0}; m < plan.size(); ++m) {
                R actPoint{plan[m]};

                fileStream << "v" << DELIMITER_OUT;
                R temp{actPoint / this->problem.Env.ScaleFactor}; 
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
      fileStream << "#Paths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        for (int j{0}; j < numRoots; ++j) {
          for (int k{0}; k < numAngles; ++k) {
            for (int l{0}; l < numAngles; ++l) {
              if (!this->neighboringMatrix.Exists(i, j, k, l)) {
                continue;  
              }
              DistanceHolder<R> &holder{this->neighboringMatrix(i, j, k, l)};
              std::deque<R> &plan{holder.Plan};
              for (int m{0}; m < plan.size() - 1; ++m) {
                R actPoint{plan[m]};
                R lastPoint{plan[m + 1]};
                if (actPoint == lastPoint) {
                  continue;
                }

                fileStream << actPoint / this->problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / this->problem.Env.ScaleFactor << "\n";   
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

template<class R>
void Solver<R, false>::saveTspPaths(const FileStruct file) {
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
        int actNode{this->tspSolution[i]};
        int nextNode{this->tspSolution[(i + 1) % numRoots]};
        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode)};
        
        if (holder.Node1 == NULL) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size(); ++k) {
          R temp{plan[k] / this->problem.Env.ScaleFactor};
          fileStream << "v" << DELIMITER_OUT;
          temp.PrintPosition(fileStream);
          fileStream << "\n";
        }
      }
      
      uint32_t pointCounter{1};
      for (int i{0}; i < numRoots; ++i) {
        int actNode{this->tspSolution[i]};
        int nextNode{this->tspSolution[(i + 1) % numRoots]};
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
      fileStream << "#TspPaths" << DELIMITER_OUT << this->problem.Dimension << "\n";;
      for (int i{0}; i < numRoots; ++i) {
        int actNode{this->tspSolution[i]};
        int nextNode{this->tspSolution[(i + 1) % numRoots]};
        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode)};

        if (holder.Node1 == NULL) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        std::deque<R> &plan{holder.Plan};
        for (int k{0}; k < plan.size() - 1; ++k) {
          fileStream << plan[k] / this->problem.Env.ScaleFactor << DELIMITER_OUT << plan[k+1] / this->problem.Env.ScaleFactor << "\n";
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
void Solver<R, true>::saveTspPaths(const FileStruct file) {
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
      fileStream << "o TspPaths\n";
      for (int i{0}; i < numRoots; ++i) {
        auto [ actNode, actAngle ] = this->gatspSolution[i];
        auto [ nextNode, nextAngle ] = this->gatspSolution[(i + 1) % numRoots];
        if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
        std::deque<R> &plan{holder.Plan};
        unsigned startingInd{vertexInd};
        for (int m{0}; m < plan.size(); ++m) {
          R actPoint{plan[m]};

          fileStream << "v" << DELIMITER_OUT;
          R temp{actPoint / this->problem.Env.ScaleFactor}; 
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
      fileStream << "#TspPaths" << DELIMITER_OUT << this->problem.Dimension << "\n";
      for (int i{0}; i < numRoots; ++i) {
        auto [ actNode, actAngle ] = this->gatspSolution[i];
        auto [ nextNode, nextAngle ] = this->gatspSolution[(i + 1) % numRoots];
        if (!this->neighboringMatrix.Exists(actNode, nextNode, actAngle, nextAngle)) {
          ERROR("Invalid TSP solution");
          exit(1);
        }

        DistanceHolder<R> &holder{this->neighboringMatrix(actNode, nextNode, actAngle, nextAngle)};
        std::deque<R> &plan{holder.Plan};
        for (int m{0}; m < plan.size() - 1; ++m) {
          R actPoint{plan[m]};
          R lastPoint{plan[m + 1]};
          if (actPoint == lastPoint) {
            continue;
          }

          fileStream << actPoint / this->problem.Env.ScaleFactor << DELIMITER_OUT << lastPoint / this->problem.Env.ScaleFactor << "\n";
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

template <class R>
void SolverBase<R>::checkDistances(std::deque<Node<R> *> &plan, double distanceToCheck) {
  Node<R> *previous{nullptr};
  double distance{0};
  for (Node<R> *node : plan) {
    if (previous != nullptr) {
      if (!IsPathFree(node->Position, previous->Position)) {
        ERROR("Path not feasible!");
        exit(1);
      }
      distance += previous->Distance(*node);
    }
    previous = node;
  }

  if (distance < distanceToCheck - SFF_TOLERANCE || distance > distanceToCheck + SFF_TOLERANCE) {
    ERROR("Distances mismatch");
    exit(1);
  }
}

template<class R>
DistanceMatrix<DistanceHolder<R>> &SolverBase<R>::GetNeighboringMatrix() {
  return this->neighboringMatrix;
}

template<class R>
TSPOrder &SolverBase<R>::GetTSPSolution() {
  return this->tspSolution;
}

#endif /* __SOLVER_H__ */
