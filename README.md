<h1 align="center">
  Far-field heat transfer and monochromatic heat currents in a cylindrical nonreciprocal cavity
</h1>

<p align="center">
  <strong>Authors:</strong> Guillem Masdemont, J. Legendre and G. T. Papadakis
  <br>
  <em>(To be submitted to Physical Review B)</em>
</p>

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
</p>


## Abstract

> **Note:** This is the source code for the paper (published in [Jornal and Arxiv]).

We investigate the consequences of breaking Kirchhoff’s law of thermal radiation in a hollow cylindrical cavity operating in the far-field regime, employing a custom specular ray-tracing algorithm. At thermal equilibrium, we show that the violation of reciprocity leads to nonzero \hl{heat rectification coefficients} between different parts of the cylinder, \hl{which can be tuned for perfect rectification and circulation}, while internal monochromatic currents vanish due to the intrinsic coupling between emission and absorption at specular surfaces. This constraint is lifted under non-equilibrium conditions, where rotational heat fluxes within the cavity can be precisely controlled by an appropriate composition of reciprocal and nonreciprocal materials. These findings open new avenues for thermal management and provide design principles for nonreciprocal photonic devices.

[**Read the full paper here**](INSERT_YOUR_DOI_OR_ARXIV_LINK_HERE)

We recommend running the codes by clicking on the following colab notebook: 


[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GuillemMasdemont/PRB-Codes/blob/main/monochromatic_currents/monochromatic_currents.ipynb)


Alternatively, for research and experimental purposes, clone the github repository and create an environment with conda (or venv): 

`conda create -n myenv python=3.12`

and install all the dependencies in this environemnt: 

`pip install requirements.txt`





