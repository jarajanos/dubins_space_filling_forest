# Space Filling Forest Star for planning with Dubins vehicles or on polynomial trajectories
The library of planning algorithms allows planning with SFF, NR-SFF\* and SFF\* with all extensions, sPRM and PRM\*, RRT and RRT\* (for single-goal planning only), Multi-T-RRT (for multi-goal planning only) and Lazy TSP with RRT or RRT\* as low-level planners (for multi-goal planning only). All the single-goal planners are able to plan in 2D and 3D Euclidean spaces (these are covered by OMPL and would not be further tested), to plan for Dubins car (covered by OMPL too) and Dubins airplane model, as well as to plan polynomial trajectories in 2D and 3D. All multi-goal planners are able to plan in 2D and 3D Euclidean spaces; SFF, NR-SFF\*, SFF\*, sPRM, PRM\* and Lazy TSP can additionally plan for Dubins car and Dubins airplane models.

# Dependencies
Some of the dependencies are already included in the `lib` folder (e. g., GDIP library, which had to be modified for proper functionality), however there are other dependencies, which you must install on your own:

 * **GCC compiler** -- CLANG can be used as well, you only need to change the `makefile`. Version of the used code is c++17, please download the proper versions of compilers on your own.
 * **GMP** -- More info [here](https://gmplib.org/).
 * **RAPID library** -- This library could not be included for its licensing terms, even though it is free to use for non-commercial projects. You can request your own copy [here](http://gamma.cs.unc.edu/OBB/). After obtaining one, insert the source code to the `lib/rapid-2.0.1` folder, it will be automatically compiled in the next step.
 * **FLANN library** -- This library was not included due to licensing terms, it can be downloaded from its official [GitHub repository](https://github.com/flann-lib/flann). Please insert the repository into the `lib/flann` folder. Compilation is then done automatically in the next step.
 * **YAML-CPP library** -- This library is used for parsing of the input configuraion file. It was not included due to licensing terms as well, but it can be downloaded from its official [GitHub repository](https://github.com/jbeder/yaml-cpp) as well. Please insert the repository into the `lib/yaml-cpp` folder. Compilation is then done automatically in the next step.
 * **TSP solvers** -- This software needs at least one TSP solver for proper functionality. This can be either LKH, or Concorde. Both can be downloaded and installed automatically by calling `make lkh`, or `make concorde`, or `make solver` for installation of both. If installing the Concorde solver, please make sure you have CGAL installed, and specify its location in the `makefile` file. You can also use external installations of the solvers, in that case specify their locations in the configuration file for each solved problem.
  
# How to compile
After satisfying all requirements given in the previous step, run compilation of all dependencies by

```
make install
```
After successfull compilation of dependencies, run

```
make
```

to create the executable `main` file in the folder `release`. If you need to debug the library, run

```
make debug
```

to create the executable `main` adjusted to debugging, in the folder `debug`. 

# How to run the compiled program
The compiled executable might be executed by calling

```
./release/main configuration-file [iteration-id]
```
The `iteration-id` argument is optional and it is used only for proper file indexing. Structure of the required configuration file is described in the dedicated manual file. You can also use prepared example files in the `tests` folder.

# Configuration file
The description of the structure of the configuration file can be found in the attached `config-file_manual.pdf` file, along with a description of the structure used for definition of points of various dimensions.
