# Superconducting Fermi Liquid with Spin-Orbit Coupling

This repository contains code and configuration files for a superconducting Fermi liquid with spin-orbit coupling.

### 📌 Calculation of the real-space Green function

The config.txt file is the configuration file where the rows correspond to the Fermi momentum, k<sub>F</sub>, the order parameter, Δ, the Rashba parameter, α, the Zeeman field, B, and the energy, E, in this order. The file poles.jl (written in Julia) calculates the poles of the momentum-space Green function, which have to be avoided during the integration in order to acquire the real-space Green function and outputs them into the file kpoles.dat. This file is then read by the Integration.f90 code, which performs the integrations required to bring the system's Green function into its configuration-space form.

### 📌 Slider

The slider.py is a Python file which builds a slider displaying the energy eigenvalues for several values of the Fermi momentum, the order parameter, the Rashba parameter and the Zeeman field (see the image below).

<p align="center">
  <img width="500" src="http://users.uoa.gr/~srigas/GitHub/Slider.png">
</p>
