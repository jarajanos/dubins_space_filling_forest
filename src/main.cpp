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
  try {
    config = YAML::LoadFile(".//test_building.yaml");
  } catch (...) {
    ERROR("Error loading configuration file!");
    exit(1);
  }

  std::string dim{config["problem"]["dimension"].as<std::string>()};
  if (dim == "2D") {
    Problem<Point2D<double>> problem;
    problem.dimension = D2;


  } else if (dim == "2DDubins") {
    Problem<Point2DDubins<double>> problem;
    problem.dimension = D2Dubins;


  } else if (dim == "3D") {
    Problem<Point3D<double>> problem;
    problem.dimension = D3;

  } else if (dim == "3DDubins") {

  } else {
    ERROR("Invalid value of \"dimension\" node in \"problem\" root node.");
    exit(1);
  }

  return 0;
}

template<class R>
void ParseFile(YAML::Node &config, Problem<R> &problem) {
  double scale;
  R pos;
  try {
    // problem
    YAML::Node node{config["problem"]};
    if (node.IsNull()) {
      throw std::invalid_argument("invalid problem node");
    }

    YAML::Node subNode{node["solver"]};
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid solver node in \"problem\" root node");
    }
    std::string value{subNode.as<std::string>()};
    if (value == "sff") {
      problem.solver = SFF;
    } else if (value == "rrt") {
      problem.solver = RRT;
    } else if (value == "lazy") {
      problem.solver = Lazy;
    } else {
      throw std::invalid_argument("unknown solver type in \"problem\" root node, use either sff, rrt or lazy");
    }

    subNode = node["optimize"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid optimize node in \"problem\" root node");
    }
    problem.optimize = subNode.as<bool>();

    subNode = node["iterations"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid iterations node in \"problem\" root node");
    }
    problem.maxIterations = subNode.as<int>();

    // delimiters
    node = config["delimiters"];
    if (!node.IsNull()) {
      subNode = node["standard"];
      if (!subNode.IsNull()) {
        Obstacle<R>::Delimiter = subNode.as<std::string>();
      }

      subNode = node["name"];
      if (!subNode.IsNull()) {
        Obstacle<R>::NameDelimiter = subNode.as<std::string>();
      }
    }

    // TSP solver for Lazy-RRT
    node = config["TSP-solver"];
    if (node.IsNull() && problem.solver == Lazy) {
      throw std::invalid_argument("missing TSP solver parameters for Lazy solver!");
    } else if (!node.IsNull()) {
      if (problem.solver != Lazy) {
        WARN("TSP solver is called only in Lazy solver algorithm -- defined TSP parameters are redundant and will not be used.");
      }

      subNode = node["path"];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid path node in \"TSP-solver\" root node");
      }
      problem.tspSolver = subNode.as<std::string>();

      subNode = node["type"];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid type node in \"TSP-solver\" root node");
      }
      problem.tspType = subNode.as<std::string>();
    }

    // algorithm node
    node = config["algorithm"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"algorithm\" root node");
    }
    subNode = node["m2r-ratio"];
    if (subNode.IsNull()) {
      throw std::invalid_argument("invalid m2r-ratio node in \"algorithm\" root node");
    }
    scale = subNode.as<double>();
    subNode = node["misses"];
    if (!subNode.IsNull()) {
      // TODO
      //Node<double>::ThresholdMisses = subNode.as<int>();
    }
    subNode = node["bias"];
    if (!subNode.IsNull()) {
      problem.priorityBias = subNode.as<double>();
    }

    // robot node
    node = config["robot"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalid \"robot\" root node");
    }
    FileStruct tempFile;
    if (GetFile(node, tempFile)) {
      throw std::invalid_argument("invalid file node in \"robot\" root node");
    }
    problem.environment.Robot = new Obstacle<R>(tempFile.fileName, tempFile.type == Obj, scale);

    // parse range
    node = config["range"];
    if (node.IsNull()) {
      throw std::invalid_argument("invalide \"range\" root node");
    }
    subNode = node["autoDetect"];
    problem.autoRange = (!subNode.IsNull() && subNode.as<bool>());

    if (!problem.autoRange) {
      subNode = node['x'];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid x node in \"range\" root node");
      }
      problem.environment.limits.Parse(subNode.as<std::string>(), scale, 0);

      subNode = node['y'];
      if (subNode.IsNull()) {
        throw std::invalid_argument("invalid y node in \"range\" root node");
      }
      problem.environment.limits.Parse(subNode.as<std::string>(), scale, 1);

      subNode = node['z'];
      if (subNode.IsNull() && (problem.dimension == D3 || problem.dimension == D3Dubins)) {
        throw std::invalid_argument("invalid z node in \"range\" root node");
      } else if (!subNode.IsNull()) {
        problem.environment.limits.Parse(subNode.as<std::string>(), scale, 2);
      }
    }

    // parse obstacles
    node = config["obstacles"];
    if (node.IsNull()) {
      problem.environment.HasMap = false;
    } else {
      problem.environment.ScaleFactor = scale;

      for (YAML::iterator it=node.begin(); it != node.end(); ++it) {
        FileStruct file;
        if (GetFile((*it), file)) {
          throw std::invalid_argument("invalid file attribute in \"obstacles\" root sequence");
        }

        subNode = (*it)["position"];
        if (subNode.IsNull()) {
          pos = R();
        } else {
          pos = R(subNode.as<std::string>());
        }

        Obstacle<R> &obst{problem.environment.Obstacles.emplace_back(file.fileName, file.type, pos, scale)};

        if (problem.autoRange) {
          problem.environment.processLimits(obst.GetRange());
        }
      }
    }

    // parse cities

    // parse goal

    // parse distances

    // parse save

    if (!problem.hasGoal && problem.priorityBias != 0 && problem.solver == RRT) {
      throw std::invalid_argument("Multi-T-RRT with bias is undefined!");
    } else if (problem.solver == Lazy && problem.priorityBias != 0) {
      throw std::invalid_argument("priority bias for Lazy solver is not implemented!");
    }
  } catch (const std::invalid_argument e) {

  }
}

bool GetFile(YAML::Node &node, FileStruct &file, int iteration, bool includeIter) {
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
