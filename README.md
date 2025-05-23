# BiqBin: A Solver for Max-Cut and Unconstrained Binary QP

**Copyright © 2021 — The BiqBin Project**  
Funded by **FWF (I 3199-N31)** and **ARRS (P2-0162)**

This program implements the algorithm presented in the publication:

> _Nicolò Gusmeroli, Timotej Hrga, Borut Lužar, Janez Povh, Melanie Siebenhofer, and Angelika Wiegele._  
> **"BiqBin: A Parallel Branch-and-Bound Solver for Binary Quadratic Problems with Linear Constraints"**  
> _ACM Trans. Math. Softw. 48, 2, Article 15 (June 2022), 31 pages._  
> [https://doi.org/10.1145/3514039](https://doi.org/10.1145/3514039)  
> [Also available on arXiv](https://arxiv.org/abs/2009.06240)

---

## License

This program is **Free Software**, released under the terms of the **GNU General Public License (GPL) v3**  
or (at your option) **any later version**.

You are free to:

- Use
- Modify
- Distribute

**without any warranty**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**.

For more details, refer to the [GNU General Public License](https://www.gnu.org/licenses/).

---

## 🔧 BiqBin Setup Requirements (Ubuntu 22.04)

### 📦 System Dependencies (Needs updating):

- `build-essential`
- `libopenblas-dev`
- `mpich`
- `python3`
- `python3-pip`

### 📦 Python Packages:

- `scipy`

## ⚙️ Installation

1. Open the `Makefile`.
2. Set the **compiler** and **BLAS/LAPACK** package paths according to your system.
3. Run:

```bash
make
```

<!-- ### 🔧 Docker

Dockerfile available for running the solver in a docker container, image can be created using Makefile command

```bash
make docker
```

--- -->

## 🚀 Usage

### Original Biqbin Maxcut Parallel solver - C only

Run with

```bash
mpirun -n num_processess ./biqbin instance_file params
```

- `num_processes`: number of processes to run the program using MPI, program needs at least 2 (1 master, and 1 worker process) to be used.
- `instance_file`: A file containing the graph (in edge list format).
- `params`: The parameter file used to configure the solver.

To run the g05_60.0 instance use this make command

---

### Python Wrapper for Biqbin Maxcut Parallel solver

Can be run with

```bash
mpiexec -n num_processes python3 run_example.py instance_file params
```

- `num_processes`: number of processes to run the program using MPI, program needs at least 2 (1 master, and 1 worker process) to be used.
- `instance_file`: A file containing the graph (in edge list format).
- `params`: The parameter file used to configure the solver.

To run the g05_60.0 instance use this make command


---

### Examples

Please check the following Python files to find how to setup biqbin solver through Python

- `run_example.py`: Example on how to run the default version of biqbin.
- `run_heuristic_example`: Example on how to inject your own heuristic function into biqbin
- `run_qubo_example`: Example on how to run QUBO instances on biqbin and transform the solution back to a QUBO solution


---

> ⚠️ **NOTE:** The default **maximum problem size** is **1024 nodes**.  
> If you need more, **edit** the value of `NMAX` in `biqbin.h` and `biqbin_data_objects.py`.

---

## 🛠️ Explanation of Parameters

_(Default values can be found in `biqbin.h`)_

| Parameter           | Description                                                                    |
| ------------------- | ------------------------------------------------------------------------------ |
| `init_bundle_iter`  | Initial number of iterations for the **bundle method**                         |
| `max_bundle_iter`   | Maximum number of iterations for the bundle method                             |
| `triag_iter`        | Number of iterations for **triangle inequality separation**                    |
| `pent_iter`         | Number of iterations for **pentagonal inequality separation**                  |
| `hept_iter`         | Number of iterations for **heptagonal inequality separation**                  |
| `max_outer_iter`    | Maximum number of **cutting plane algorithm** iterations                       |
| `extra_iter`        | Additional iterations for refinement or fallback in cutting plane algorithm    |
| `violated_TriIneq`  | Threshold for **triangle inequality violation**: `B(X) - 1 > violated_TriIneq` |
| `TriIneq`           | Maximum number of **triangle inequalities** added during separation            |
| `adjust_TriIneq`    | Whether to **adjust triangle inequalities dynamically** (0 or 1)               |
| `PentIneq`          | Number of **pentagonal inequalities** to add (usually `3 * Pent_Trials`)       |
| `HeptaIneq`         | Number of **heptagonal inequalities** to add (usually `4 * Hepta_Trials`)      |
| `Pent_Trials`       | Number of **simulated annealing trials** for pentagonal inequalities           |
| `Hepta_Trials`      | Number of **simulated annealing trials** for heptagonal inequalities           |
| `include_Pent`      | Include **pentagonal inequalities** in SDP bound (0 or 1)                      |
| `include_Hepta`     | Include **heptagonal inequalities** in SDP bound (0 or 1)                      |
| `root`              | If `1`, compute only **SDP bound at the root node**                            |
| `use_diff`          | If `1`, **only add cutting planes** when necessary to speed up B&B             |
| `time_limit`        | Maximum runtime in **seconds**. If `0`, runs until optimal solution is found   |
| `branchingStrategy` | Branching strategy:<br>`0 = LEAST_FRACTIONAL`<br>`1 = MOST_FRACTIONAL`         |
---