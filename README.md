# Simulation Study: Asymptotic Fluctuations of Submatrices of Wigner Matrix Powers

## Overview
This repository contains the simulation codebase (R and C++) developed for my MSTAT dissertation at the Indian Statistical Institute (ISI) under the supervision of Dr. Soumendu Sundar Mukherjee. The theoretical core of the dissertation investigates the universal Gaussian fluctuations of complex, multivariate statistics derived from random matrices—specifically, the non-spectral observables of powers of Gaussian Wigner matrices. 

This simulation was built to empirically validate the Central Limit Theorems (CLT) derived in the dissertation, specifically visualizing the asymptotic normality of linear combinations (sums) of entries within principal submatrices of $X_N^p$, where $X_N$ is a Gaussian Wigner matrix.

## Objectives: What We Extract and Visualize
The primary goal of this codebase is to translate the theoretical bounds (derived via Stein's Method and the quantitative Cramér-Wold device) into empirical visualizations. Specifically, the simulation is designed to:

1. **Extract the Submatrix Sum Statistic:** Generate an $N \times N$ Gaussian Orthogonal Ensemble (GOE) matrix, raise it to a power $p$, extract the principal $k \times k$ submatrix, and compute the total sum of its entries.
2. **Test Submatrix Dimension Constraints ($k$ vs $N$):** Observe the statistical behavior of the sum as a function of the submatrix growth rate $k$.
3. **Visualize Asymptotic Normality:** Construct empirical distributions (histograms) over thousands of Monte Carlo replications and overlay the theoretical standard normal density curve to visually assess convergence.
4. **Statistical Verification:** Compute the Shapiro-Wilk test $p$-value to quantify the departure from normality.
## Repository Contents

### 1. `wigner_submatrix.cpp` (C++ Backend)
A high-performance C++ script utilizing the `RcppArmadillo` library to handle heavy matrix operations. 
* **Matrix Generation:** Simulates an $N \times N$ matrix of i.i.d. standard normals and symmetrizes it to create a properly scaled GOE matrix ($X_{ij} \sim \mathcal{N}(0, 1/N)$ and $X_{ii} \sim \mathcal{N}(0, 2/N)$ ).
* **Eigen-decomposition:** Efficiently computes the $p$-th power of the matrix using symmetric eigen-decomposition ($W^p = V D^p V^T$).
* **Extraction:** Dynamically defines the submatrix dimension $k$ (e.g., $k = \lfloor\sqrt{N}\rfloor$) and extracts the sum of the $k \times k$ principal block, returning it to R.

### 2. `Wigner power submatrix sums.R` (R Frontend)
The R script responsible for driving the simulation, managing memory, and plotting.
* **Parameters:** Sets the matrix dimension ($N = 1000$), matrix power ($p = 4$), and Monte Carlo iterations ($iterations = 1000$).
* **Replication:** Uses `replicate()` to repeatedly call the C++ backend and build a large sample vector of submatrix sums.
* **Standardization:** Centers and scales the resulting vector of sums to a mean of $0$ and standard deviation of $1$.
* **Outputs:** Generates the final overlaid histogram and outputs the Shapiro-Wilk $p$-value to verify Gaussianity.

## Dependencies
To run this code, you will need the following R packages installed:
* `Rcpp`
* `RcppArmadillo`

*Note: You must have a working C++ compiler installed to compile the Armadillo matrix algebra code.*

## 🚀 How to Run
1. Clone this repository.
2. Update the `setwd()` path in the R script to match your local directory.
3. Modify the submatrix size $k$ inside `wigner_submatrix.cpp` if you wish to test different growth limits (e.g., change `int k = (int) std::sqrt((double)N);` to `int k = N / 10;`).
4. Run `Wigner power submatrix sums.R`.
