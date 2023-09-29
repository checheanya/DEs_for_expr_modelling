# Simulating the gene expression using systems of differential equations and chemical kinetics
_by Anna Chechenina, BII, July-September 2023_

The goal of the project is to build modelling setup for robust simulations of bulk gene expression patterns given initial parameter values. In order to build such model we designed a set of chemical reactions that describes events leading to gene expression activation or repression. By adjusting rate constants and initial concentration of the species in the system, one can get simulation of number of transcripts over time. <br>

## One-equation modifications
Firstly, we attempted to model expression by modifying transcriptional model based on Gillespie stochastic algorithm presented in [1]. The initial system contains three ordinary differential equations with five kinetic parameters and considers the following events: promotor activation, promotor inactivation, transcription with amplification, and RNA decay:

$$\begin{cases}
dA_i(t)/dt = k_{on,i}I_i(t) - k_{off,i}A_i(t) \\
dI_i(t)/dt = k_{off,i}A_i(t) - k_{on,i}I_i(t) \\
dx_i(t)/dt = \phi_i s_i A_i(t) - \delta_i x_i(t)
\end{cases}$$

In our model we want to consider next level of regulation by introducing enhancer and silencer terms as positive and negative regulators of the transcription. We created two modifications of the transcription equation of theinitial model:

First version:
$$\frac{dx_i(t)}{dt} = p_i \cdot s_i \cdot A_i(t) \cdot \frac{E_i(t)^n}{E_i(t)^n + E_{0i}^n} \cdot \frac{S_{0i}^m}{S_i(t)^m + S_{0i}^m} - q_i \cdot x_i(t)$$

Second version:
$$\frac{dx_i(t)}{dt} = p_i \cdot s_i \cdot A_i(t) \cdot \frac{(E_i(t)/K_e1 + 1)(A_i(t)+K_e2)}{(A_i(t)+K_s2)(1+S_i(t)/K_s1)} - q_ix_i(t)$$


build a model using [GillesPy2](https://github.com/StochSS/GillesPy2/tree/main) library in python

This repository contains code for building the model and reproducing plots for results analysis. There are three folders containing separate parts of the project:

**1. ATAC-seq data distribution and noise exploration.** 
This folder contains code to perform analysis on selected datasets (see details in the notebooks) and to explore corresponding distributions.  

### Installation

To get the scripts clone the git repository:

```bash
git clone https://github.com/checheanya/DEs_for_expr_modelling.git && cd DEs_for_expr_modelling
```

Create a `conda/mamba` environment with necessary packages and activate it:

```bash
conda env create -f environment.yml
conda activate DE_dash
```

### Usage

To run the script, just call it from the directory where the tool is located:

```bash
python lineplots_2equations.py 
```
for plotting the line plots and

```bash
python scatterplots_3equations.py 
```
for plotting the scatterplots.


Then proceed to the link that appears on the screen.

Sourses:
1. Piras, V., Tomita, M. & Selvarajoo, K. Transcriptome-wide Variability in Single Embryonic Development Cells. Sci Rep 4, 7137 (2014). https://doi.org/10.1038/srep07137
