<h1 align="center">
Far-field heat transfer and monochromatic thermal currents in a cylindrical nonreciprocal cavity
</h1>

<p align="center">
  <strong>Authors:</strong> Guillem Masdemont, J. Legendre and G. T. Papadakis
  <br>
</p>

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
</p>


> **Note:** This is the source code for the paper "Far-field heat transfer and monochromatic thermal currents in a cylindrical
nonreciprocal cavity", which is currently under review.
[**Read the full paper here**](Link)

## Abstract

We investigate the consequences of breaking Kirchhoff’s law of thermal radiation in a hollow cylindrical cavity operating in the far-field regime, employing a custom specular ray-tracing algorithm. At thermal equilibrium, we show that the violation of reciprocity leads to nonzero heat rectification coefficients between different parts of the cylinder, which can be tuned for perfect rectification and circulation, while internal monochromatic currents vanish due to the intrinsic coupling between emission and absorption at specular surfaces. This constraint is lifted under non-equilibrium conditions, where rotational heat fluxes within the cavity can be precisely controlled by an appropriate composition of reciprocal and nonreciprocal materials. These findings open new avenues for thermal management and provide design principles for nonreciprocal photonic devices.


## Running codes 

To reproduce the results and figures presented in the article without installing local dependencies, we recommend running the code interactively via Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GuillemMasdemont/PRB-Codes/blob/paralleloptimization/Codes/Simulation_transmission_coefficient/notebook_TransmissionSimulation.ipynb)

Alternatively, for research and experimental purposes, you can run the code locally. First clone the repository

```bash
git clone https://github.com/GuillemMasdemont/Far-field-heat-transfer-and-monochromatic-heat-currents-in-a-cylindrical-nonreciprocal-cavity.git PRB-CylindricalNonreciprocalCodes
cd PRB-CylindricalNonreciprocalCodes
```

and create an environment with conda (or venv):

`conda create -n myenv python=3.12`
`conda activate myenv`

Once your environment is active, install the required libraries:

`pip install -r requirements.txt`

## Workflow & Reproducibility
To reproduce the figures and results from the article, follow the steps in the provided Jupyter Notebooks:

1.  **Define Parameters:** Set the emission, absorption, and reflection profiles for all the elements of the cylinder. 
2.  **Run the Simulation:** The notebook runs a specular ray tracing algorithm to compute the transmission coefficient between the elements of the cylinder and the heat currents arising in the far-field regime. 
3.  **Analyze Results:** The final cells output the **Heat Rectification Coefficient** and the **Monochromatic Heat Current**, generating the plots shown in the manuscript.

All intermediate data is handled automatically within the notebook workflow.





