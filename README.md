# Picasso: Memory Efficient Parallel Graph Coloring for Quantum Chemistry

This repository provides sequential and GPU source code for the graph coloring algorithm introduced in the following publication:

```bibtex
@inproceedings{ferdous2024picasso,
  title={Picasso: Memory-Efficient Graph Coloring Using Palettes With Applications in Quantum Computing},
  author={Ferdous, S M and Neff, Reece and Peng, Bo and Shuvo, Salman and Minutoli, Marco and Mukherjee, Sayak and Kowalski, Karol and Becchi, Michela and Halappanavar, Mahantesh},
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
- C++ compiler (e.g., g++) with C++17 support
  - g++ 9.1.0 is the earliest tested version known to work
- CMake (version >= 3.1)
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
- **CPU-only version**: `build/apps/palcolEr` (source: `apps/appPalColEncRecDir.cc`)
- **GPU-accelerated version**: `build/apps/palcolGr` (source: `apps/appPalColGpuMemRec.cc`)

Use the `-h` option with either executable to view available parameters and usage instructions.

### Examples

**Example 1:** GPU version help
```bash
./build/apps/palcolGr -h
```

**Example 2:** CPU-only execution

Generate grouping for the `Pauli_ket_ccsd_H2_631g.json` dataset using the recursive option (`-r`) with an initial target of 15 colors:

```bash
./build/apps/palcolEr -t 15 -r --in data/pauli_ket_ccsd_data/Pauli_ket_ccsd_H2_631g.json --out data/group_results/H2_631g_out1.json
```

Make sure the `data/group_results` directory exists before running this command.

Expected screen output:
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

**Example 3:** GPU execution

Generate grouping for the `Pauli_ket_ccsd_H2_631g.json` dataset using GPU with
the recursive option (`-r`) with an initial target of 15 colors:

```bash
./build/apps/palcolGr -t 15 -r --in data/pauli_ket_ccsd_data/Pauli_ket_ccsd_H2_631g.json
```

Expected screen output:
```
Using 32-bit offsets
Temp Storage Bytes: 767
Current Bytes: 360
Fits: 12312 < 31328
***********Level 0*******
Num Nodes: 89
Num Edges: 0
Avg. Deg.: 0
Palette Size: 15
List Size: 4
Num Conflict Edges: 1539
Conflict to Edge (%): inf
Num Colors: 15
Assign Time: 1.27645
Conf. Build Time: 0.0173272
Conf. Color Time: 0.000178425

Final Num invalid Vert: 23
Naive Color TIme: 8.87737e-06
# of Final colors: 26
```

**Example 4:** Large GPU execution

Generate grouping for the `Pauli_ket_ccsd_H4_1D_631g.json` dataset using GPU with
the recursive option (`-r`), alpha=2 (`-a 2`), and an initial target of 1% of the total number
of nodes (`-t 0.01`).

```bash
./build/apps/palcolGr -t 0.01 -a 2 -r --in data/pauli_ket_ccsd_data/Pauli_ket_ccsd_H4_1D_631g.json
```

Expected screen output of the last level:
```
***********Level 67*******
Num Nodes: 106
Num Edges: 0
Avg. Deg.: 0
Palette Size: 1
List Size: 1
Num Conflict Edges: 2852
Conflict to Edge (%): inf
Num Colors: 3943
Assign Time: 1.46385e-05
Conf. Build Time: 0.0126362
Conf. Color Time: 3.08678e-05

Final Num invalid Vert: 99
Naive Color TIme: 0.00393682
# of Final colors: 4004
```

## Contact

For questions or further information, please contact:

- **S M Ferdous:** [sm.ferdous@pnnl.gov](mailto:sm.ferdous@pnnl.gov), [ferdous.csebuet@gmail.com](mailto:ferdous.csebuet@gmail.com)
- **Reece Neff:** [rwneff@ncsu.edu](mailto:rwneff@ncsu.edu), [reece.neff@pnnl.gov](mailto:reece.neff@pnnl.gov)
