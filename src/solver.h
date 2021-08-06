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
    bool optimize;
    bool usePriority;
    double treeDistance;
    double collisionSampleSize{0.1};

    Environment<R> &env;

    //std::deque<Node<T, R> *> allNodes;
    //std::deque<Tree<T, Node<T, R>>> trees;
    //SymmetricMatrix<DistanceHolder<T, Node<T, R>>> neighboringMatrix;
    //std::deque<Tree<T, Node<T, R>> *> connectedTrees;

    bool isPathFree(R &start, R &finish);

    virtual void getPaths() = 0;
    virtual void getAllPaths();
    virtual void smoothPaths() = 0;
    void checkPlan();
    // void checkPlan(std::deque<Node<T, R> *> &plan);
    // void checkDistances(std::deque<Node<T, R> *> &plan, T distanceToCheck);
    // double computeDistance(std::deque<Node<T, R> *> &plan);

    virtual void saveIterCheck(const int iter);
    virtual void saveTrees(const FileStruct file);
    virtual void saveCities(const FileStruct file);
    virtual void savePaths(const FileStruct file);
    virtual void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);
    virtual void saveTsp(const FileStruct file);
};

#endif /* __SOLVER_H__ */
