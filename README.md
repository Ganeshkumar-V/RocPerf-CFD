# Propellant Regression Framework

![Simulation GIF](media.gif)

An OpenFOAM-based collection of applications, solvers, utilities, and functionObjects for modelling propellant regression, particle-laden flows, and rocket motor-related post-processing. This repository bundles example applications (e.g., the `rocketMotor` application), a set of custom `functionObjects`, `fvPatchFields`, interfacial models, and utilities commonly used for propellant combustion/regression studies.

This codebase is provided as source and is intended to be built and run on Linux platforms with a compatible C++ toolchain and OpenFOAM build tools.

## ✨ Key Features

* **Custom Application:** `applications/rocketMotor/rocketMotor` - a
  specialised application for rocket motor simulations and propellant
  regression handling.

* **Post-Processing:** A collection of `functionObjects` useful for post-processing
  (thrust, vorticity/stream function, etc.).

* **Utilities:** Tools under `utilities/` for mapping fields, post-processing
  rocket-related quantities, and initializing rocket-specific fields.

* **Build System:** Platform-specific makefiles (`platforms/`) and per-target
  Makefiles (`Make/`) for `wmake`.

## 📸 Screenshots / Demo

*(Optional: Add simulation results, graphs, or GIFs here.)*

`![Simulation Result](https://placehold.co/600x400/222/fff?text=Simulation+Result+Screenshot)`

## 🛠️ Requirements

* **Platform:** Linux 

* **Core:** OpenFOAM (developed and tested on **v2112**. This version is recommended, though it may work on other versions.)

* **Compiler:** C++ (GCC or Intel)

* **Build System:** `make` and OpenFOAM toolchain (wmake, Allwmake)

Note: This repository contains a `bashrc` at the project root to set up the environment. The `platforms/` directory contains the library and the applications build on this project

## 🚀 Quick Start — Compilation

1. Ensure OpenFOAM v2112 is installed and source its environment before building this project.

2. Clone the repository in your preferred parent directory (e.g., your projects folder), then enter the project directory:

    ```sh
    git clone https://github.com/Ganeshkumar-V/Propellant-Regression-Framework.git
    cd Propellant-Regression-Framework
    ```

3.  Source the project `bashrc` to set up the environment:

    ```sh
    source bashrc
    ```

4.  Compile the libraries, applications, and utilities:

    ```sh
    ./Allwmake
    ```

    (Optional: run `./Allwclean` before building to clean the project.)

## 🎓 Tutorial Case Files

This repository provides four sample cases under `tutorials/`:

* CD nozzle — gas phase

* CD nozzle — gas + particle two-phase

* Rocket motor — with propellant regression

* Rocket motor — without propellant regression

Each tutorial folder contains a `README` with case-specific setup, runtime, and post-processing instructions.

## 🗺️ Repository Layout

* `applications/` — Custom applications (e.g., `rocketMotor`).

* `src/` — Source modules (`functionObjects`, `fvPatchFields`, `interfacialModels`, etc.).

* `utilities/` — Additional tools and post-processing utilities.

* `platforms/` — Platform-specific compilation configurations.

* `Make/` — Per-module build configurations used by `wmake`/`Allwmake`.

* `bashrc` — Environment setup script.

* `LICENSE` — GNU General Public License v3 (GPL-3.0).
## 🧪 Tests and Validation

A `validation/` directory will be added, containing full validation cases and reference data.

* `validation/JPL-gas-particle-nozzle` — Gas–particle JPL nozzle flow validation case.

* `validation/NASA-CDV` — NASA CDV validation case.

Each case includes a `README` describing setup, run commands, and post-processing.

## 🤝 Contributing

Contributions are welcome. Please open an issue to describe the change, then fork the repository and create a pull request. Follow the existing code style and respect the project's license (GPL-3.0).

## 📄 License

Distributed under the GNU General Public License v3 (GPL-3.0). See the `LICENSE` file for the full text.

This software is an independent project and is not affiliated with, endorsed by, or associated with OpenCFD Ltd or the OpenFOAM® project. OpenFOAM® is a registered trademark of OpenCFD Ltd.

## 🧑‍🔬 Authors and Acknowledgements

* **Original author:** Ganeshkumar V

* **Developed under the supervision of:** Prof. Dilip Srinivas Sundaram

## 📜 Citation

If you use this software in published work, please cite the associated paper:

Venukumar, Ganeshkumar and Sundaram, Dilip Srinivas, "Computational Study of Propulsive Performance of Frozen Nano-Aluminum and Water (ALICE) Mixtures," *Journal of Propulsion and Power*, vol. 41, no. 3, pp. 330–346, 2025. https://doi.org/10.2514/1.B39541

**Suggested BibTeX:**

```bibtex
@article{doi:10.2514/1.B39541,
author = {Venukumar, Ganeshkumar and Sundaram, Dilip Srinivas},
title = {Computational Study of Propulsive Performance of Frozen Nano-Aluminum and Water (ALICE) Mixtures},
journal = {Journal of Propulsion and Power},
volume = {41},
number = {3},
pages = {330-346},
year = {2025},
doi = {10.2514/1.B39541}
}
```

## 📬 Contact

If you encounter problems, please open an issue on the repository.