# Picasso: Memory Efficient Parallel Graph Coloring for Quantum Chemistry

This repository provides sequential and GPU source code for the graph coloring algorithm introduced in the following publication:

```bibtex
@inproceedings{ferdous2024picasso,
  title={Picasso: Memory-Efficient Graph Coloring Using Palettes With Applications in Quantum Computing},
  author={Ferdous, SM and Neff, Reece and Peng, Bo and Shuvo, Salman and Minutoli, Marco and Mukherjee, Sayak and Kowalski, Karol and Becchi, Michela and Halappanavar, Mahantesh},
  booktitle={2024 IEEE International Parallel and Distributed Processing Symposium (IPDPS)},
  pages={241--252},
  year={2024},
  organization={IEEE}
}
```

The graph coloring problem addressed here is motivated by quantum chemistry applications, specifically for efficient measurement through unitary partitioning of Pauli strings. Please cite the above publication if you use Picasso in your research.

## Installation

### Dependencies
Ensure you have the following dependencies installed:
- C++ compiler (e.g., g++)
- CMake (version >= 3.18)
- CUDA (version >= 11; tested, earlier versions may also work)

### Building Picasso
Clone and build the repository:

```bash
git clone git@github.com:smferdous1/Picasso.git --recurse-submodules
cd Picasso

# Create and enter build directory
mkdir build
cd build

# Run CMake configuration
cmake ..

# Compile the project
make
```

## Dataset

A sample dataset is included in JSON format under `data/pauli_ket_ccsd_data`. This dataset represents the hydrogen molecule (H) with varying numbers of bases and basis functions. The complete dataset used in our paper will soon be made available, with a link provided here.

## Running the Code

Two primary executables are provided:
- CPU-only version: `build/apps/palcolEr` corresponds to `apps/appPalColEncRecDir.cc` 
- GPU-accelerated version: `build/apps/palcolGr` corresponds to `apps/appPalColGpuMemRec.cc`

Use the `-h` option with either executable to see available parameters and usage instructions.

Example 1:
```bash
./build/apps/palcolGr -h
```
Example 2:
This example shows  how to execute the CPU-only coloring code to generate grouping for the `Pauli_ket_ccsd_H2_631g.json` file. We will use the recursive version (`-r`) and will also output the grouping into a new json file (`--out`). We execute the following command into the base directory. Here, our initial target color value is 15 (The actual coloring will be more). 
```bash
./build/apps/palcolEr -t 15 --in data/pauli_ket_ccsd_data/Pauli_ket_ccsd_H2_631g.json --out data/group_results/H2_631g_out1.json
```
The expected output on the screen could be

```
The grouping output: data/group_results/H2_631g_out1.json
using list greedy ordering for conflict coloring

***********Level 0*******
Num Nodes: 89
Num Edges: 1988
Avg. Deg.: 44.6742
Palette Size: 15
List Size: 4
Num Conflict Edges: 1539
Conflict to Edge (%): 77.4145
Num Colors: 15
Assign Time: 2.80626e-05
Conf. Build Time: 0.000160687
Conf. Color Time: 0.000172731

Final Num invalid Vert: 23
Naive Color TIme: 8.6166e-06
# of Final colors: 26
```
The corresponding grouping JSON file should be created in the `data/group_results` directory (make sure that the directory exists!)

## Contact

For questions or further information, please contact:

- **S M Ferdous:** [sm.ferdous@pnnl.gov](mailto:sm.ferdous@pnnl.gov), [ferdous.csebuet@gmail.com](mailto:ferdous.csebuet@gmail.com)
- **Reece Neff:** [rwneff@ncsu.edu](mailto:rwneff@ncsu.edu), [reece.neff@pnnl.gov](mailto:reece.neff@pnnl.gov)
