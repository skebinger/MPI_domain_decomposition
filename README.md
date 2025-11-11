# 📦 MPIDCL — MPI Decomposition Library

**MPIDCL** (MPI Decomposition Library) is a lightweight **Fortran** library for structured domain decomposition in **1D, 2D, and 3D** using **MPI**. It’s designed for **CFD simulations** and other stencil-based solvers on structured grids, with built-in ghost cell exchange and consistent subdomain indexing across ranks.

Whether you're new to MPI or developing a scalable parallel solver, `mpidcl` offers a clear, modular, and efficient way to handle distributed grids.

---

## ✨ Features

- 🧩 **Block Decomposition**  
  Static rectangular domain splitting in 1D, 2D, or 3D.

- 🔁 **Ghost Cell Exchange**  
  Built-in halo exchange routines for structured grids.

- 📦 **Consistent Indexing**  
  Subdomains are assigned globally consistent index bounds.

- 📐 **Modular Design**  
  Clean separation between halo packing, buffer management, and communication routines.

- 🚀 **Scalable**  
  Supports arbitrary Cartesian process layouts using MPI.

---

## 🛠 Requirements

- **Fortran 2008** (or newer)
- **MPI library** (e.g., OpenMPI, MPICH, Intel MPI)
- **CMake**
- **[pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit)** (optional, for unit testing)
- **[FORD](https://github.com/Fortran-FOSS-Programmers/ford)** (optional for generating documentation)

> **Note:**  
> A setup script `install_pfUnit_MPI.sh` is provided under `./test` to install **pFUnit** using your chosen MPI implementation (currently Intel MPI or OpenMPI).
> The project uses **FORD** for generating documentation from the provided code. I can be easily installed via pip. For more detailed installation instructions please visit the respective repository page.

---

## ⚙️ Building the Library

To build the library, a helper script `makeall.sh` is provided at the repository root. The library can be compiled as either a **static** or **shared** library.

Compiler flags and MPI settings are managed via dedicated CMake extension files under `./config`.

> The final library artifact (`libmpidcl`) will be placed under `./build/lib`.

### 🔧 Example

To compile a **static** library **without tests**, using **Release** flags:

```bash
./makeall.sh Release OFF ON OFF
