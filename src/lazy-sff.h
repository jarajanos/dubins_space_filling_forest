/**
 * @file lazy-sff.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 1.0
 * @date 15. 02. 2022
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __LAZY_SFF__
#define __LAZY_SFF__

#include "solver.h"
#include "forest.h"
#include "opendubins/dubins.h"
#include "opendubins/dubins3D.h"

template<class R>
class LazySpaceForest : public Solver<R> {
  using Solver<R>::Solver;
  public:
    void Solve() override;

  protected:
    std::deque<Point3DDubins> pathPoints;
    double finalDistance{std::numeric_limits<double>::max()};

    void getPaths() override;    

    void saveTspPaths(const FileStruct file) override;
    void saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);
};

template<> void LazySpaceForest<Point3DDubins>::Solve();
template<> void LazySpaceForest<Point3DDubins>::saveTspPaths(const FileStruct file);
template<> void LazySpaceForest<Point3DDubins>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime);

template<class R>
void LazySpaceForest<R>::Solve() {
  ERROR("Lazy-SFF defined for 3D Dubins problems only!");
}

template<class R>
void LazySpaceForest<R>::getPaths() {

}

template<class R>
void LazySpaceForest<R>::saveTspPaths(const FileStruct file) {

}

template<class R>
void LazySpaceForest<R>::saveParams(const FileStruct file, const int iterations, const bool solved, const std::chrono::duration<double> elapsedTime) {

}

#endif
