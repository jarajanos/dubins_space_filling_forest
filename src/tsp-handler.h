/**
 * @file tsp-handler.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz) -- original code by Robert
 * Penicka
 * @brief
 * @version 1.0
 * @date 07. 11. 2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#ifndef __TSP_HANDLER_H__
#define __TSP_HANDLER_H__

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "common.h"
#include "constants.h"
#include "problem.h"

typedef std::vector<std::vector<double>> LengthMatrix;
typedef std::vector<std::vector<int>> IntegerMatrix;
typedef std::deque<int> TSPOrder;
typedef std::deque<std::tuple<int, int>> GATSPOrder;

template <class R> class TSPMatrix {
public:
  LengthMatrix data;

  TSPMatrix() {}
  TSPMatrix(int size);
  TSPMatrix(int size, double value);
  TSPMatrix(Problem<R> &problem, DistanceMatrix<DistanceHolder<R>> &neighboringMatrix);

  TSPMatrix TransformGATSPtoATSP();
  GATSPOrder TransformATSPSolToGATSP(TSPOrder &atspSeq);

  TSPMatrix TransformATSPtoTSP();
  TSPOrder TransformTSPSolToATSP(TSPOrder &tspSolution);

  TSPOrder SolveByConcorde();
  TSPOrder SolveByLKH();

protected:
  Problem<R> *problem;

  int dubinsResolution{-1};
  int size;
  double maxLength{TSP_TOLERANCE};
  double bigM;

  bool isATSP{false};
  bool isTSP{false};

  IntegerMatrix getIntegerMatrix();
};

template<> TSPMatrix<Point2DDubins>::TSPMatrix(Problem<Point2DDubins> &problem, DistanceMatrix<DistanceHolder<Point2DDubins>> &neighboringMatrix);

class Wrapper {
public: 
  Wrapper(std::string solver_dir, std::string run_dir);
  std::string WriteTSPLIBFile(std::string &fname_basis,
                                 IntegerMatrix &valuematrix,
                                 std::string &user_comment,
                                 bool isATSP);
  virtual TSPOrder ReadResultCmd(std::string fname_basis) = 0;
  virtual bool RunSolverCmd(std::string fname_basis, bool silent) = 0;
  bool RmSolutionFileCmd(std::string fname_basis);
protected:
  std::string solver_dir;
  std::string tsplib_dir;
};

class LKHWrapper : public Wrapper {
public:
  using Wrapper::Wrapper;

  TSPOrder ReadResultCmd(std::string fname_basis) override;
  bool RunSolverCmd(std::string fname_basis, bool silent) override;
};

class ConcordeWrapper : public Wrapper {
public:
  using Wrapper::Wrapper;

  TSPOrder ReadResultCmd(std::string fname_basis) override;
  bool RunSolverCmd(std::string fname_basis, bool silent) override;
};

template <class R> TSPMatrix<R>::TSPMatrix(int size) : size{size} {
  data.resize(size);
  for (int i{0}; i < size; ++i) {
    data[i].resize(size);
  }
}

template <class R> TSPMatrix<R>::TSPMatrix(int size, double value) : size{size} {
  data.resize(size);
  for (int i{0}; i < size; ++i) {
    data[i].resize(size, value);
  }
}

template <class R>
TSPMatrix<R>::TSPMatrix(Problem<R> &problem, DistanceMatrix<DistanceHolder<R>> &neighboringMatrix)
    : TSPMatrix(neighboringMatrix.GetSize(), -1) {
      this->problem = &problem;
      
      isTSP = true;
      for (int i{0}; i < size; ++i) {
        for (int j{0}; j < i; ++j) {
          if (neighboringMatrix.Exists(i, j)) {
            double dist{neighboringMatrix(i, j).Distance};
            data[i][j] = dist;
            data[j][i] = dist;

            maxLength = MAX(maxLength, dist);
          }
        }
      }

      bigM = maxLength * size;
      // replate all nonexisting paths with M
      for (int i{0}; i < size; ++i) {
        for (int j{0}; j < size; ++j) {
          if (data[i][j] == -1) {
            data[i][j] = bigM;
          }
        }
      }
    }


// Noon-Bean transformation
template<class R>
TSPMatrix<R> TSPMatrix<R>::TransformGATSPtoATSP() {
  if (isATSP || isTSP) {
    ERROR("Trying to convert to ATSP a matrix which is already ATSP/TSP");
    exit(1);
  } else if (dubinsResolution == -1) {
    ERROR("Missing resolution of Dubins");
    exit(1);
  }

  // for nonexisting paths, calculate even bigger M
  double biggerM{bigM * size};

  // init new matrix
  TSPMatrix<R> atsp(size, bigM);
  atsp.problem = this->problem;
  atsp.dubinsResolution = dubinsResolution;
  atsp.isATSP = true;
  atsp.maxLength = bigM + maxLength;
  atsp.bigM = atsp.maxLength * size;

  int dimension{size / dubinsResolution};
  for (int i{0}; i < dimension; ++i) {
    for (int k{0}; k < dubinsResolution; ++k) {
      int fromId{i * dubinsResolution + k};
      for (int j{0}; j < dimension; ++j) {
        for (int l{0}; l < dubinsResolution; ++l) {
          int toId{j * dubinsResolution + l};
          int toTfId{j * dubinsResolution + ((l + 1) % dubinsResolution)};

          if (i == j && k == l) {
            // zero cycle
            atsp.data[fromId][toTfId] = 0;
          } else if (i != j) {
            atsp.data[fromId][toTfId] = this->data[fromId][toId] + bigM;  // + bigM to assure zero cycle is covered first
          }
        }
      }
    }
  }

  return atsp;
}

template<class R>
GATSPOrder TSPMatrix<R>::TransformATSPSolToGATSP(TSPOrder &atspSeq) {
  GATSPOrder gatspSolution;

  int seqSize{(int)atspSeq.size()};
  int last{atspSeq[0] / dubinsResolution};
  for (int i{seqSize - 1}; i >= 0; --i) {
    int val{atspSeq[i]};
    int clusterInd{val / dubinsResolution};
    if (clusterInd != last) {
      gatspSolution.push_front(std::make_pair(val / dubinsResolution, val % dubinsResolution));
      last = clusterInd;
    }
  }

  return gatspSolution;
}

template<class R>
TSPMatrix<R> TSPMatrix<R>::TransformATSPtoTSP() {
  if (!isATSP || isTSP) {
    ERROR("Trying to convert to TSP a matrix which is GATSP/TSP");
    exit(1);
  }

  double biggerM{3 * this->bigM};
  TSPMatrix tsp(3 * size, biggerM);
  tsp.problem = this->problem;
  tsp.dubinsResolution = dubinsResolution;
  tsp.isTSP = true;
  tsp.maxLength = biggerM;
  tsp.bigM = biggerM;

  for (int x{0}; x < size; ++x) {
    for (int y{0}; y < size; ++y) {
      tsp.data[2 * size + x][y] = this->data[x][y];
      tsp.data[y][2 * size + x] = this->data[x][y];
    }
  }

  for (int x{0}; x < size; ++x) {
    tsp.data[2 * size + x][x] = biggerM;
    tsp.data[x][2 * size + x] = biggerM;
  
    tsp.data[size + x][x] = 0;
    tsp.data[2 * size + x][size + x] = 0;
    tsp.data[x][size + x] = 0;
    tsp.data[size + x][2 * size + x] = 0;
  }

  return tsp;
}

template<class R>
TSPOrder TSPMatrix<R>::TransformTSPSolToATSP(TSPOrder &tspSolution) {
  TSPOrder atspSolution;
  for (auto &a : tspSolution) {
    if (a < tspSolution.size() / 3) {
      atspSolution.push_back(a);
    }
  }

  if (std::abs(tspSolution[0] - tspSolution[1]) == tspSolution.size() / 3) {
    if (tspSolution[0] > tspSolution[1]) {
      std::reverse(atspSolution.begin(), atspSolution.end());
    }
  } else {
    if (tspSolution[1] > tspSolution[2]) {
      std::reverse(atspSolution.begin(), atspSolution.end());
    }
  }

  return atspSolution;
}

template<class R>
TSPOrder TSPMatrix<R>::SolveByConcorde() {
  TSPMatrix<R> *tsp;
  TSPMatrix<R> dummy;
  if (isTSP) {
    tsp = this;
  } else if (isATSP) {
    dummy = TransformATSPtoTSP();
    tsp = &dummy;
  } else {
    ERROR("Trying to solve TSP from GATSP -- convert matrix to ATSP/TSP first");
    exit(1);
  }

  ConcordeWrapper wrapper(this->problem->TspSolver, TEMP_DIR);

  FileStruct tempTsp;
  tempTsp.fileName = TEMP_TSP;
  std::string id{"id_" + std::to_string(this->problem->Repetition) + "_"};
  FileStruct runFile{PrefixFileName(tempTsp, id)};

  IntegerMatrix valueMatrix{tsp->getIntegerMatrix()};

  wrapper.WriteTSPLIBFile(runFile.fileName, valueMatrix, id, false);
  bool run_ret{wrapper.RunSolverCmd(runFile.fileName, false)};
  TSPOrder solution{wrapper.ReadResultCmd(runFile.fileName)};
  bool delete_ret{wrapper.RmSolutionFileCmd(runFile.fileName)};

  if (isTSP) {
    return solution;
  } else {
    return TransformTSPSolToATSP(solution);
  }
}

template<class R>
TSPOrder TSPMatrix<R>::SolveByLKH() {
  TSPMatrix *tsp;
  TSPMatrix dummy;
  if (!isTSP && !isATSP) {
    ERROR("Trying to solve TSP from GATSP -- convert matrix to ATSP/TSP first");
    exit(1);
  } else {
    tsp = this;
  }

  LKHWrapper wrapper(this->problem->TspSolver, TEMP_DIR);

  FileStruct tempTsp;
  tempTsp.fileName = TEMP_TSP;
  std::string id{"id_" + std::to_string(this->problem->Repetition) + "_"};
  FileStruct runFile{PrefixFileName(tempTsp, id)};

  IntegerMatrix valueMatrix{tsp->getIntegerMatrix()};

  wrapper.WriteTSPLIBFile(runFile.fileName, valueMatrix, id, tsp->isATSP);
  bool run_ret{wrapper.RunSolverCmd(runFile.fileName, false)};
  TSPOrder solution{wrapper.ReadResultCmd(runFile.fileName)};
  bool delete_ret{wrapper.RmSolutionFileCmd(runFile.fileName)};

  return solution;
}

template<class R>
IntegerMatrix TSPMatrix<R>::getIntegerMatrix() {
  IntegerMatrix valueMatrix;

  int maxValue;
  if (this->size <= 10) {
    maxValue = 1 << 10;
  } else if (this->size <= 100) {
    maxValue = 1 << 12;
  } else {
    maxValue = TSP_MAX;
  }

  valueMatrix.resize(this->size);
  double resolution{maxValue / maxLength};
  for (int rowid = 0; rowid < this->size; ++rowid) {
    valueMatrix[rowid].resize(this->size);
    for (int numid = 0; numid < this->size; ++numid) {
      valueMatrix[rowid][numid] = (int)(this->data[rowid][numid] * resolution);
    }
  }

  return valueMatrix;
}

#endif
