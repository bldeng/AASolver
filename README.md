## AASolver: Anderon Acceleration for Geometry Optimization and Physics Simulation

This is the source code for the examples in the following paper:

* Yue Peng, Bailin Deng, Juyong Zhang, Fanyu Geng, Wenjie Qin, and Ligang Liu. 2018. [Anderson acceleration for geometry optimization and physics simulation](https://arxiv.org/abs/1805.05715). ACM Trans. Graph. 37, 4, Article 42 (July 2018).


### Compilation

The code has been tested with the following systems and compilers:

* Ubuntu 18.04 using GCC 7.3;
* macOS High Sierra 10.13.6, using Apple Clang 9.1 (without OpenMP support) and GCC 8.2 (with/without OpenMP support).


Compile the code using the following steps:

1. Before compiling the code, make sure [Eigen](http://http://eigen.tuxfamily.org) (version 3.3+) is installed, using either of the following methods:

	* Download the package from the Eigen homepage, unzip it into a folder `eigen` within the folder `AndersonOpt/external`. Make sure the file `external/eigen/Eigen/Dense` can be found.
	* Alternatively, use a package manager of your operating system to install Eigen, e.g.:
		* on Ubuntu, use the command `apt-get install libeigen3-dev`;
		* on macOS, install [homebrew](https://brew.sh/) and run `brew install eigen`.

2. GLUT is required to compile the simulation code.

3. Create a folder `build` within the folder `AndersonOpt`.

4. From within the `build` folder, run the following commands:  
	
	```
	$ cmake -DCMAKE_BUILD_TYPE=Release ..
	$ make
	```
	
	* On mac, the default Apple clang compiler does not support OpenMP. To enable OpenMP, first install the homebrew GCC using command `brew install gcc`. This should install the latest gcc with associated commands gcc-X and g++-X where X is the version number, e.g. gcc-8 and g++-8. Then use the following commands to compile the code with homebrew GCC:
	
	```
	$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 ..
	$ make
	```
	
5. Afterwards, there should be five executables in the `build` folder: `PlanarityOpt`, `WireMeshOpt`, `ARAP2D`, `ARAP3D` and `TetMeshSimulation`.

6. Test models are provided in folder `Test_files`.


### Usage

* Use the following commands for planar quad mesh and wire mesh optimization respectively:
	```
	$ PlanarityOpt INPUT_POLY_MESH REF_TRI_MESH OPTION_FILE OUPUT_MESH
	```
	```
	$ WireMeshOpt INPUT_POLY_MESH REF_TRI_MESH OPTION_FILE OUTPUT_MESH
	```
	*  `INPUT_POLY_MESH` is the input mesh to be optimized.
	*  `REF_TRI_MESH` is the triangle mesh for the reference shape.
	*  `OPTION_FILE` is a text file storing the algortihm parameters. An example file is provided as `Options.txt`.
	*  `OUTPUT_MESH` is the output mesh file name.

* To perform 2D as-rigid-as-possible optimization:
	```
	$ ARAP2D PROBLEM_DEF_MESH INIT_MESH HANDLE_IDX_FILE HANDLE_COORDS_FILE OPTION_FILE OUTPUT_MESH
	```
	*  `PROBLEM_DEF_MESH` is a mesh that defines the rest shape.
	*  `INIT_MESH` is the initial shape of the mesh to be optimized.
	*  `HANDLE_IDX_FILE` is a text file storing the handle vertex indices.
	*  `HANDLE_COORDS_FILE` is a file that stores the 2D coordinates of handle vertices.

* To perform 3D as-rigid-as-possible optimization:
	```
	$ ARAP3D PROBLEM_DEF_MODEL_NAME INIT_COORD_MODEL_NAME HANDLE_IDX_FILE HANDLE_COORDS_FILE OPTION_FILE OUTPUT_MESH
	```
	* `PROBLEM_DEF_MODEL_NAME` is the name (without suffix) for the ele and node files that describe the tet mesh of the rest shape.
	* `INIT_COORD_MODEL_NAME` is the name (without suffix) for the ele and node files that describe the initial shape of the tet mesh to be optimized.
	* `HANDLE_IDX_FILE` is a text file storing the handle vertex indices.
	* `HANDLE_COORDS_FILE` is a file that stores the 3D coordinates of handle vertices.
	*  `OUTPUT_MESH` is a mesh that stores the boundary surface of the optimized tet mesh.

* To perform physics simulation for elastic tet mesh:
	```
	$ TetMeshSimulation MODEL_NAME HANDLE_IDX_FILE OPTION_FILE
	``` 
	* `INIT_COORD_MODEL_NAME` is the name (without suffix) for the ele and node files that describe the initial shape of the tet mesh.
	* `HANDLE_IDX_FILE` is a text file storing the handle vertex indices. The handle vertices are fixed during simulation.
	
	The simulation is controlled using the following keys:
	
	* Press `s` to simulate one frame.
	* Press `1` to start the demo.
	* Press `e` to suspend demo. 
	* Press `Esc` to exit the program.

* Below are example commands using the data provided in `Test_files`.
	* `$ ARAP2D bar.obj bar_Init.obj bar_Handle.txt xy_bar.txt Options.txt output.obj`
	* `$ ARAP3D boy boy_Init boy_handle.txt xyz_boy.txt Options.txt output.obj`
	* `$ PlanarityOpt bunny_poly.obj bunny_ref.obj Options.txt output.obj`
	* `$ WireMeshOpt MaleTorso.obj MaleTorso_target.obj Options.txt output.obj`
	* `$ TetMeshSimulation armadillo armadilloH.txt Options.txt`

### License
The code is released under BSD 3-Clause License.
	
### Contact
Please contact Bailin Deng <<bldeng@gmail.com>> or Yue Peng <<echoyue@mail.ustc.edu.cn>> if you have any comments or questions.	
 

	