#pragma once
/**
 * @file random-generator.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __RANDOM_GENERATOR_H__
#define __RANDOM_GENERATOR_H__

#include <random>
#include <math.h>
#include <vector>
#include <chrono>

#include "common.h"
#include "problem.h"

#include "RapidTrajectoryGenerator.h"
#include "RapidTrajectoryGenerator2D.h"
using namespace RapidQuadrocopterTrajectoryGenerator;

#define TO_SIGN(A)  ((A) * 2 - 1)

typedef std::mt19937_64 randomEngine;   // set used random engine

template <class R>
class RandomGenerator {
  public:
    RandomGenerator(Problem<R> &problem);

    bool RandomPointInDistance(const R& center, R& point, const double distance);
    void RandomPointInSpace(R& point);
    int RandomIntMinMax(const int min, const int max);
    double RandomProbability();

  private:
    Problem<R> &problem;
    Range limits;
    int maxIter;
    randomEngine rndEng;

    std::uniform_real_distribution<double> uniDistAngle;
    std::uniform_real_distribution<double> uniDistPitch;
    std::uniform_real_distribution<double> uniSpaceX;
    std::uniform_real_distribution<double> uniSpaceY;
    std::uniform_real_distribution<double> uniSpaceZ;
    std::uniform_real_distribution<double> uniProb;
    std::normal_distribution<double> normProb;
};

template<class R>
RandomGenerator<R>::RandomGenerator(Problem<R> &problem) 
  : limits{problem.Env.Limits}, maxIter{problem.MaxMisses}, problem{problem} {
  // seed with actual time
  std::chrono::time_point<std::chrono::high_resolution_clock, std::chrono::nanoseconds> tSeed{std::chrono::high_resolution_clock::now()};
  std::uint_fast64_t uSeed{static_cast<uint_fast64_t>(tSeed.time_since_epoch().count())};

  rndEng = randomEngine(uSeed);
  //INFO(uSeed);
  //rndEng = randomEngine(1644768989840546946); // deterministic behaviour

  uniDistAngle = std::uniform_real_distribution<double>(-M_PI, M_PI);
  uniSpaceX = std::uniform_real_distribution<double>(limits.mins[0], limits.maxs[0]);
  uniSpaceY = std::uniform_real_distribution<double>(limits.mins[1], limits.maxs[1]);
  uniSpaceZ = std::uniform_real_distribution<double>(limits.mins[2], limits.maxs[2]);
  uniDistPitch = std::uniform_real_distribution<double>(limits.mins[3], limits.maxs[3]);
  uniProb = std::uniform_real_distribution<double>(0, 1);
  normProb = std::normal_distribution<double>(0, 1);
}

template<class R>
int RandomGenerator<R>::RandomIntMinMax(const int min, const int max) {
  std::uniform_int_distribution<int> uniInt(min, max);
  return uniInt(rndEng);
}

template <class R>
double RandomGenerator<R>::RandomProbability() {
  return uniProb(rndEng);
}

#endif /*__RANDOM_GENERATOR_H__*/
