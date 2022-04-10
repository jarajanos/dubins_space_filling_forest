#pragma once
/**
 * @file constants.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 04. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <string>

#ifndef DELIMITER_OUT
#define DELIMITER_OUT         (" ")
#endif

#ifndef CSV_DELIMITER
#define CSV_DELIMITER         (",")
#endif

#ifndef CSV_DELIMITER_2
#define CSV_DELIMITER_2       (";")
#endif

#ifndef TSP_DELIMITER
#define TSP_DELIMITER         (" ")
#endif

#define CSV_NO_PATH           "NaN"

// used for LazyTSP
#ifndef TEMP_TSP
#define TEMP_TSP              "tempTsp"
#endif

#ifndef TEMP_DIR
#define TEMP_DIR              "temp/"
#endif

#define TSP_TOLERANCE         1e-5
#define TSP_MAX_10            (1 << 10)
#define TSP_MAX_100           (1 << 12)
#define TSP_MAX               (1 << 20)
#define DEFAULT_CONCORDE      "solver/concorde/concorde/TSP/concorde"
#define DEFUALT_LKH           "solver/lkh/LKH-3.0.6/LKH"

#define DEFAULT_THRES_MISS    3
#define SFF_TOLERANCE         1e-9

#define EQ_TOLERANCE          1e-6

#define DEFAULT_SAMP_DIST     1
#define DEFAULT_MIN_THRUST    0
#define DEFAULT_MAX_THRUST    1
#define DEFAULT_MAX_ROTSPEED  1
#define DEFAULT_CTRL_INT      0.01
#define AVG_VEL_STD_MULT      0.5

#define DEF_GRAVITY           9.81

#define FLANN_NUM_KD_TREES    4
#define FLANN_NUM_SEARCHES    128
#define FLANN_PREC_MULTIPLIER 1.1

#define CHECK_CONNECTED_ITER  50

const std::string WHITESPACE = " \n\r\t\f\v";

#endif 
