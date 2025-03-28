# STYX Hydrological Model

STYX is a distributed hydrological model written in modern C++ (C++17). It simulates surface and subsurface water flow, heat, and solute transport and includes a network component (e.g., stormwater or pipe systems). The model uses OpenMP for parallelization, VTK ASCII files for grid input (surface, subsurface, and network), and CSV files for material libraries and initial/boundary conditions. A JSON configuration file is used to control simulation parameters.

---

## Table of Contents

- [Features](#features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Profiling & Optimization](#profiling--optimization)
- [Contributing](#contributing)
- [Roadmap](#roadmap)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Features

- **Modular Architecture:**  
  The code is organized into separate modules for simulation, mesh generation, network connectivity, input/output, and physics.  
- **Multi-Domain Modeling:**  
  Supports surface grids (2D), subsurface grids (3D), and network models (e.g., wells and pipes).
- **Modern C++:**  
  Written in C++17 using RAII, smart pointers, and robust exception handling.
- **Input/Output:**  
  Uses VTK ASCII format for grid geometries and CSV for material libraries and initial/boundary conditions.
- **JSON Configuration:**  
  Simulation settings are provided via a JSON file.
- **Parallel Processing:**  
  OpenMP is used for computational tasks.
- **Profiling Support:**  
  Built-in timing (using custom RAII Timer objects) to profile expensive routines (e.g., mesh connectivity).

---

## Project Structure

```
STYX/
├── include/          # Header files (Simulation, Mesh classes, Network, InputManager, Timer, etc.)
├── src/              # Source files for the simulation and modules
├── config/           # JSON configuration file (e.g., settings.json)
├── data/             # Input files: VTK grids, CSV material libraries, boundary/initial conditions, etc.
└── Makefile          # Linux Makefile for building the project
```

---

## Installation

### Requirements

- **Compiler:** A C++ compiler supporting C++17 (e.g., g++ 7 or later).
- **OpenMP:** Ensure OpenMP is available on your system.
- **Dependencies:**  
  - [nlohmann/json](https://github.com/nlohmann/json) (header-only JSON library)
  
### Build Instructions

1. **Clone the Repository:**

   ```bash
   git clone git clone https://extgit.vtt.fi/lassi.warsta/styx.git
   cd styx
   ```

2. **Build the Project:**

   Use the provided Makefile:
   ```bash
   make
   ```
   This will compile all source files and produce the executable `styx_model`.

---

## Usage

1. **Configure the Simulation:**  
   Edit the JSON file in the `config/` folder (e.g., `config/settings.json`) to set simulation length, time step, file paths, and other parameters.

2. **Run the Simulation:**

   ```bash
   ./styx_model
   ```
   
3. **Visualization:**  
   The simulation outputs grid data in VTK format (and results in CSV), which can be visualized using ParaView or other compatible tools.

---

## Profiling & Optimization

The model includes timer instrumentation (using a custom RAII `Timer` class) in key modules such as mesh connectivity routines and the simulation loop. When you run the simulation, you will see console output with elapsed times for functions like:

- `SurfaceMesh::buildConnectivity`
- `SubsurfaceMesh::buildConnectivity`
- Total simulation run time

These timers help identify bottlenecks so you can further optimize performance. For example, the network connectivity and 3D mesh connectivity routines have been optimized using hash maps rather than \(O(n^2)\) comparisons.

---

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a feature branch (e.g., `feature/my-new-feature`).
3. Commit your changes with clear messages.
4. Submit a merge request.

Please ensure your code adheres to the project's coding standards and include tests if possible.

---

## Roadmap

- **Improve Parallelization:**  
  Explore further OpenMP or alternative parallelism for connectivity and solver routines.
- **Extend Physics Modules:**  
  Enhance heat and solute transport modules with additional state variables and coupling.
- **Robust I/O:**  
  Integrate more robust CSV/VTK parsers if needed.
- **User Interface:**  
  Develop a simple GUI or command-line interface for easier configuration and monitoring.
- **Documentation:**  
  Expand documentation with detailed developer guides and API references.

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Special thanks to the developers of [nlohmann/json](https://github.com/nlohmann/json) for their excellent JSON library.
- Inspired by previous versions of the STYX model and contributions from the hydrological modeling community.

---

*This README is a starting point. Contributions to expand and improve this documentation are welcome!*

