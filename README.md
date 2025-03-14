# Picasso: Memory Efficient Parllel Graph Coloring for Quantum Chemistry


## Introduction

This repository provides the sequential and GPU source code of the graph coloring algorithm proposed in the following publication.

    @inproceedings{ferdous2024picasso,
      title={Picasso: Memory-Efficient Graph Coloring Using Palettes With Applications in Quantum Computing},
      author={Ferdous, SM and Neff, Reece and Peng, Bo and Shuvo, Salman and Minutoli, Marco and Mukherjee, Sayak and Kowalski, Karol and Becchi, Michela and Halappanavar, Mahantesh},
      booktitle={2024 IEEE International Parallel and Distributed Processing Symposium (IPDPS)},
      pages={241--252},
      year={2024},
      organization={IEEE}

}
The graph coloring problem is motivated from a quantum chemistry application, where unitary partioning of a set of Pauli strings are used for efficient measurement. Please cite the above publication if you use Picasso in your research. 


## Installation
To use this program, you need to have the following dependencies installed:
- C++ compiler (e.g., g++)
- CMake >=3.18
- CUDA >=11 (tested, earlier versions may be compatible)

Once you have the dependencies installed, you can follow these steps to install the program:

1. Clone the repository:
   ```bash
   git clone git@github.com:smferdous1/Picasso.git --recurse-submodules
   cd Picasso
   ```
  
  2. Create the build directory, run cmake and  make:
		    
		```bash
		# Create and enter build directory
		mkdir build
		cd build

		# Run CMake configuration
		cmake ..

		# Compile the project
		make
		```

