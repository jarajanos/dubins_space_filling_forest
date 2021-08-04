/**
 * @file main.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz)
 * @brief 
 * @version 2.0
 * @date 02. 08. 2021
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "main.h"

int main(int argc, char *argv[]) {
  YAML::Node config;
  return 0;
}

void ParseFile(const std::string &fileName, Problem<double> &problem) {
  
}

bool GetFile(YAML::Node *node, FileStruct &file, int iteration, bool includeIter) {
//   rapidxml::xml_attribute<> *attr;
//   if (node == nullptr) {
//     return true;
//   }
  
//   attr = node->first_attribute("file");
//   if (attr == nullptr) {
//     return true;
//   }  
//   file.fileName = attr->value();
//   if (iteration != 0 && includeIter) {
//     size_t iter{file.fileName.find_last_of('.')};
//     file.fileName.insert(iter, '_' + std::to_string(iteration));
//   }

//   attr = node->first_attribute("is_obj");
//   if (attr == nullptr || !strcmp(attr->value(), "false")) {
//     file.type = Map;
//   } else if (!strcmp(attr->value(), "true")) {
//     file.type = Obj;
//   } else {
//     throw std::invalid_argument("invalid attribute isObj in file node!");
//   }

  return false;
}
