#pragma once
/**
 * @file main.h
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __MAIN_H__
#define __MAIN_H__

#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include "yaml-cpp/yaml.h"

#include "common.h"
#include "environment.h"
#include "problem.h"

int main(int argc, char *argv[]);
void ParseFile(const std::string &fileName, Problem<double> &problem);
bool GetFile(YAML::Node node, FileStruct &file, int iteration=0, bool includeIter=true);

#endif /*__MAIN_H__*/
