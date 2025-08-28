# P0-Radical-Pair-Recombination

Fortran simulations of **spin dynamics** in the **central spin problem** and **radical pair** systems.

This repo aims to **reproduce Fig. 2 and Fig. 3** from:

> D. E. Manolopoulos & P. J. Hore, *An improved semiclassical theory of radical pair recombination reactions*, **J. Chem. Phys.** 139, 124106 (2013).

---

## Overview

- Language/standard: **Fortran 2023**
- Parallelism: **OpenMP** (configured for **5 performance cores**)
- Purpose: verification/reproduction of published results; comparison of exact QM vs semiclassical dynamics

## Files

- **`QM_nN_fast.f90`** – Fastest and most memory-efficient exact QM dynamics.  
  Needed for large systems (e.g. **16 nuclei**); avoids explicit dense matrices.  
  The Hilbert space for 16 nuclei + 1 electron is **2^17 = 131072** dimensions, which exceeds typical memory without efficient algorithms.

- **`QM_n4.f90`** – Standard exact QM approach for small systems (e.g. 4 nuclei).  
  Builds explicit operator matrices and propagates with matrix algebra.

- **`SC_n16.f90`** – **Manolopoulos-Hore** improved semiclassical dynamics.  
  Cost grows **linearly** with number of nuclei (vs **exponential** for exact QM).  
  Typically faster than QM for **≥12 nuclei**; coherent oscillations become less prominent as N increases.
  
- **`SW_n16.f90`** – **Schulten-Wolynes** semiclassical dynamics.  
  Closed form expression, but fails to capture the correct long-time behaviour of the electron spin dynamics. 

---

## Build

You’ll need a recent Fortran compiler with OpenMP (e.g. **gfortran 14+**).

```bash
# Exact QM (large N, fast path)
gfortran -O2 -std=f2023 -fopenmp -o QM_nN_fast QM_nN_fast.f90

# Exact QM (small N, educational)
gfortran -O2 -std=f2023 -o QM_n4 QM_n4.f90

# Semiclassical
gfortran -O2 -std=f2023 -o SC_n16 SC_n16.f90
