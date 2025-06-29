# Ellipsoid Geodesic-Patch Volume Calculator

A Python tool to **compute** and **visualize** the exact volume between a true geodesic patch on an ellipsoid and its straight-edge chord triangle.

![Volume Between Geodesic Patch and Flat Triangle](images/volume_plot.png)

---

## Features

- **Ellipsoid volume**  
  Computes the volume of the curved geodesic patch on an axis-aligned ellipsoid  
  \[
    \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1
  \]
  bounded by three surface points \(A,B,C\).

- **Flat-triangle volume**  
  Calculates the tetrahedral volume of the straight-edge triangle and reports  
  \(\Delta V = V_{\rm patch} - V_{\rm base}\).

- **High-precision integration**  
  Uses `mpmath` for tunable-precision surface integrals, including built-in error estimates.

- **Mesh verification**  
  Builds a closed triangular mesh (patch + base + side-walls) and verifies the result via the divergence theorem.

- **3D visualization**  
  Renders with Matplotlib:
  - **Red mesh:** true geodesic-patch interior  
  - **Blue face:** flat triangle base  
  - **Gray wireframe:** full ellipsoid  

---

## Installation

```bash
git clone https://github.com/yourusername/ellipsoid-volume.git
cd ellipsoid-volume
pip install numpy mpmath matplotlib


## Usage

Run the main script:
```bash
python plot.py

    Adjust ellipsoid axes (a,b,c) and points A, B, C directly in plot.py.

    Configure precision and mesh resolution at the top of the file:

    mp.mp.dps   = 50    # digits of precision
    n_phi, n_theta = 80, 80  # mesh resolution

    The script outputs the computed volumes in the console and opens an interactive 3D plot.

