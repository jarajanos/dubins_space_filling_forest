/**
 * @file tsp-handler.cpp
 * @author Jaroslav Janos (janosjar@fel.cvut.cz) -- original code by Robert
 * Penicka
 * @brief
 * @version 1.0
 * @date 07. 11. 2021
 *
 * @copyright Copyright (c) 2021
 *
 */

#include "tsp-handler.h"

template <>
TSPMatrix<Point2DDubins>::TSPMatrix(Problem<Point2DDubins> &problem,
    DistanceMatrix<DistanceHolder<Point2DDubins>> &neighboringMatrix)
    : TSPMatrix(neighboringMatrix.GetSize() * neighboringMatrix.GetResolution(), -1) {
        this->dubinsResolution = neighboringMatrix.GetResolution();
        this->problem = &problem;

        int dimension{neighboringMatrix.GetSize()};
        int resolution{neighboringMatrix.GetResolution()};
        
        for (int i{0}; i < dimension; ++i) {
          for (int j{0}; j < dimension; ++j) {
            for (int k{0}; k < resolution; ++k) {
              for (int l{0}; l < resolution; ++l) {
                if (neighboringMatrix.Exists(i, j, k, l)) {
                  double dist{neighboringMatrix(i, j, k, l).Distance};
                  if (dist != std::numeric_limits<double>::max()) {
                    data[i * resolution + k][j * resolution + l] = dist;

                    maxLength = MAX(maxLength, dist);
                  }
                }
              }
            }
          }
        }

        bigM = maxLength * size;
        // replace nonexisting paths with M
        for (int i{0}; i < size; ++i) {
          for (int j{0}; j < size; ++j) {
            if (data[i][j] == -1) {
              data[i][j] = bigM;
            }
          }
        }
      }

template <>
TSPMatrix<Point3DDubins>::TSPMatrix(Problem<Point3DDubins> &problem,
    DistanceMatrix<DistanceHolder<Point3DDubins>> &neighboringMatrix)
    : TSPMatrix(neighboringMatrix.GetSize() * neighboringMatrix.GetResolution(), -1) {
        this->dubinsResolution = neighboringMatrix.GetResolution();
        this->problem = &problem;

        int dimension{neighboringMatrix.GetSize()};
        int resolution{neighboringMatrix.GetResolution()};
        
        for (int i{0}; i < dimension; ++i) {
          for (int j{0}; j < dimension; ++j) {
            for (int k{0}; k < resolution; ++k) {
              for (int l{0}; l < resolution; ++l) {
                if (neighboringMatrix.Exists(i, j, k, l)) {
                  double dist{neighboringMatrix(i, j, k, l).Distance};
                  if (dist != std::numeric_limits<double>::max()) {
                    data[i * resolution + k][j * resolution + l] = dist;

                    maxLength = MAX(maxLength, dist);
                  }
                }
              }
            }
          }
        }

        bigM = maxLength * size;
        // replace nonexisting paths with M
        for (int i{0}; i < size; ++i) {
          for (int j{0}; j < size; ++j) {
            if (data[i][j] == -1) {
              data[i][j] = bigM;
            }
          }
        }
      }

Wrapper::Wrapper(std::string solver_dir, std::string run_dir) {
  std::filesystem::path current_path = std::filesystem::current_path();

  std::filesystem::path full_lkh_dir = current_path;
  full_lkh_dir += "/" + solver_dir;
  this->solver_dir = full_lkh_dir.string(); // relative path
  this->tsplib_dir = run_dir;        // absolute path
}

std::string Wrapper::WriteTSPLIBFile(std::string &fname_basis,
                               std::vector<std::vector<int>> &valuematrix,
                               std::string &user_comment,
                               bool isATSP) {
  if (!std::filesystem::exists(tsplib_dir)) {
    std::filesystem::create_directory(tsplib_dir);
  }
  int dims_tsp = valuematrix.size();
  std::string name_line = "NAME : " + fname_basis + "\n";
  std::string comment_line = "COMMENT : " + user_comment + "\n";
  std::string tsp_line = isATSP ? "TYPE: ATSP\n" : "TYPE: TSP\n";
  std::string dimension_line = "DIMENSION : " + std::to_string(dims_tsp) + "\n";
  std::string edge_weight_type_line =
      "EDGE_WEIGHT_TYPE : EXPLICIT\n"; // explicit only
  std::string edge_weight_format_line = "EDGE_WEIGHT_FORMAT: FULL_MATRIX\n";
  std::string display_data_type_line =
      "DISPLAY_DATA_TYPE: NO_DISPLAY\n"; // 'NO_DISPLAY'
  std::string edge_weight_section_line = "EDGE_WEIGHT_SECTION\n";
  std::string eof_line = "EOF\n";
  std::stringstream Cost_Matrix_STRline;
  for (int i = 0; i < dims_tsp; ++i) {
    std::stringstream cost_matrix_strline;
    int j = 0;
    for (; j < dims_tsp - 1; ++j) {
      cost_matrix_strline << valuematrix[i][j] << " ";
    }
    j = dims_tsp - 1;
    cost_matrix_strline << valuematrix[i][j] << '\n';
    Cost_Matrix_STRline << cost_matrix_strline.str();
  }

  // write the tsp file
  std::stringstream tsp_out_fname;
  tsp_out_fname << tsplib_dir << fname_basis << ".tsp";
  std::ofstream outfile(tsp_out_fname.str());
  outfile << name_line;
  outfile << comment_line;
  outfile << tsp_line;
  outfile << dimension_line;
  outfile << edge_weight_type_line;
  outfile << edge_weight_format_line;
  outfile << edge_weight_section_line;
  outfile << Cost_Matrix_STRline.str();
  outfile << eof_line;
  outfile.close();

  // write the param file
  std::stringstream tsp_outpar_fname;
  tsp_outpar_fname << tsplib_dir << fname_basis << ".par";
  std::ofstream param_outfile(tsp_outpar_fname.str());

  std::string problem_file_line = "PROBLEM_FILE = " + tsplib_dir +
                                  fname_basis +
                                  ".tsp\n"; // remove pwd + tsplib_dir
  std::string move_type_line = "MOVE_TYPE = 5\n";
  std::string patching_c_line = "PATCHING_C = 3\n";
  std::string patching_a_line = "PATCHING_A = 2\n";
  std::string runs_line = "RUNS = 10\n";
  std::string tour_file_line =
      "TOUR_FILE = " + tsplib_dir + fname_basis + ".txt" + "\n";

  param_outfile << problem_file_line;

  param_outfile << move_type_line;
  param_outfile << patching_c_line;
  param_outfile << patching_a_line;
  param_outfile << runs_line;
  param_outfile << tour_file_line;
  param_outfile.close();
  return tour_file_line;
}


bool Wrapper::RmSolutionFileCmd(std::string fname_basis) {
  std::stringstream rm_sol_cmd;
  rm_sol_cmd << "rm " << tsplib_dir << fname_basis << ".txt " << tsplib_dir << fname_basis << ".par " << tsplib_dir << fname_basis << ".tsp";
  int retval = std::system(rm_sol_cmd.str().c_str());
  return retval == 0;
}

///////////////////////////////// LKH section ////////////////////////////

bool LKHWrapper::RunSolverCmd(std::string fname_basis,
                                   bool silent = false) {
  std::stringstream run_lkh_cmd;
  run_lkh_cmd << solver_dir << " " << tsplib_dir << fname_basis
              << ".par";

  if (silent) {
    run_lkh_cmd << " > /dev/null 2>&1";
  }

  int retval = std::system(run_lkh_cmd.str().c_str());
  return retval == 0;
}

TSPOrder LKHWrapper::ReadResultCmd(std::string fname_basis) {
  TSPOrder sequence;
  std::stringstream outfilename;
  outfilename << tsplib_dir << fname_basis << ".txt";
  std::ifstream fileout(outfilename.str());
  std::string line;
  int linei = 0;
  if (fileout.is_open()) {
    while (std::getline(fileout, line)) {
      if (linei >= 6) { // first lines are useless
        int seq;
        std::istringstream ss(line);
        ss >> seq;
        if (seq == -1) {
          break;
        }
        sequence.push_back(seq - 1);
      }
      linei += 1;
    }
    fileout.close();
  }

  return sequence;
}

//////////////////// Concorde section //////////////////////////////////

bool ConcordeWrapper::RunSolverCmd(std::string fname_basis,
                                   bool silent = false) {
  std::stringstream run_conc_cmd;
  run_conc_cmd << solver_dir << " -o " << tsplib_dir << fname_basis << ".txt -x " << tsplib_dir << fname_basis
              << ".tsp";

  if (silent)
  {
    run_conc_cmd << " > /dev/null 2>&1";
  }

  int retval = std::system(run_conc_cmd.str().c_str());
  return retval == 0;
}

TSPOrder ConcordeWrapper::ReadResultCmd(std::string fname_basis) {
  TSPOrder sequence;
  std::stringstream outfilename;
  outfilename << tsplib_dir << fname_basis << ".txt";
  std::ifstream fileout(outfilename.str());

  if (fileout.is_open()) {
    int size;
    fileout >> size;
    
    int a;
    for (int i{0}; i < size; ++i) {
      fileout >> a;
      sequence.push_back(a);
    }
    fileout.close();
  }

  return sequence;
}

