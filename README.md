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


## Introduction 

> **Note:** This is the source code for the paper (published in [Jornal and Arxiv]).
[**Read the full paper here**](INSERT_YOUR_DOI_OR_ARXIV_LINK_HERE)


## Abstract


We investigate the consequences of breaking Kirchhoff’s law of thermal radiation in a hollow cylindrical cavity operating in the far-field regime, employing a custom specular ray-tracing algorithm. At thermal equilibrium, we show that the violation of reciprocity leads to nonzero heat rectification coefficients between different parts of the cylinder, which can be tuned for perfect rectification and circulation, while internal monochromatic currents vanish due to the intrinsic coupling between emission and absorption at specular surfaces. This constraint is lifted under non-equilibrium conditions, where rotational heat fluxes within the cavity can be precisely controlled by an appropriate composition of reciprocal and nonreciprocal materials. These findings open new avenues for thermal management and provide design principles for nonreciprocal photonic devices.


## Running codes 

To reproduce the results and figures presented in the article without installing local dependencies, we recommend running the code interactively via Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GuillemMasdemont/PRB-Codes/blob/main/monochromatic_currents/monochromatic_currents.ipynb)

Alternatively, for research and experimental purposes, you can run the code locally. First clone the repository

```bash
git clone [https://github.com/GuillemMasdemont/PRB-Codes.git](https://github.com/GuillemMasdemont/PRB-Codes.git)
cd PRB-Codes
```

and create an environment with conda (or venv):

`conda create -n myenv python=3.12`
`conda activate myenv`

Once your environment is active, install the required libraries:

`pip install -r requirements.txt`

## Workflow & Reproducibility
To reproduce the figures and results from the article, follow the steps in the provided Jupyter Notebooks:

1.  **Define Parameters:** Set the emission profile (temperature distributions) and material properties in the initial cells.
2.  **Run the Simulation:** The code calculates the electromagnetic Green's functions for the cylindrical cavity.
3.  **Analyze Results:** The final cells output the **Transmission Coefficient** $\mathcal{T}(\omega)$ and the **Monochromatic Heat Current** $J(\omega)$, generating the plots shown in the manuscript.

All intermediate data is handled automatically within the notebook workflow.

## Citation
If you use this code or data in your own research, please cite our paper:

> **Guillem Masdemont**, Julien Legendre, Georgia Papadakis.  
> *"Far-field heat transfer and monochromatic heat currents in a cylindrical nonreciprocal cavity"* > [Journal Name/arXiv], [Year]. [DOI Link]

**BibTeX:**
```bibtex
@article{YourLastName2026,
  title = {Far-field heat transfer and monochromatic heat currents in a cylindrical nonreciprocal cavity},
  author = {YourLastName, FirstName and CoAuthor, Name},
  journal = {Physical Review B},
  year = {2026},
  doi = {10.1103/PhysRevB.XX.XXXXXX}
}





