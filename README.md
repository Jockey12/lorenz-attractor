# Lorenz Attractor 3D Visualization

A Python application that visualizes the Lorenz Attractor in 3D with an interactive GUI for modifying initial conditions and system parameters.

## Description

The Lorenz Attractor is a set of chaotic solutions to the Lorenz system of differential equations. This application provides an interactive 3D visualization where you can modify:

- **Initial Conditions**: x0, y0, z0 (starting position in 3D space)
- **Lorenz Parameters**:
  - σ (sigma): Prandtl number (default: 10.0)
  - ρ (rho): Rayleigh number (default: 28.0)
  - β (beta): Geometric factor (default: 8/3 ≈ 2.667)

## Features

- Interactive 3D plot using matplotlib
- GUI with input fields for all parameters
- Real-time updates when parameters change
- Clear visualization of the chaotic butterfly-shaped attractor
- Red marker showing the initial point

## Requirements

- Python 3.6 or higher
- NumPy
- SciPy
- Matplotlib

## Installation

1. Clone this repository:
```bash
git clone https://github.com/Jockey12/lorenz-attractor.git
cd lorenz-attractor
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the application:
```bash
python lorenz_attractor.py
```

The GUI will open with:
- Left panel: Input fields for initial conditions and parameters
- Right panel: 3D visualization of the Lorenz Attractor

To modify the visualization:
1. Enter your desired values in the input fields
2. Click the "Update Plot" button
3. The 3D plot will refresh with the new parameters

## The Lorenz System

The Lorenz system is described by three differential equations:

```
dx/dt = σ(y - x)
dy/dt = x(ρ - z) - y
dz/dt = xy - βz
```

Where:
- σ (sigma) relates to the Prandtl number
- ρ (rho) relates to the Rayleigh number
- β (beta) is a geometric factor

## Examples

### Classic Lorenz Attractor
- Initial conditions: (1.0, 1.0, 1.0)
- Parameters: σ=10, ρ=28, β=2.667

### Different Initial Conditions
Try experimenting with different initial conditions to see how the trajectory changes while still following the same attractor pattern.

## License

This project is open source and available under the MIT License.