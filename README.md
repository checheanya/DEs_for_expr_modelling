# Simulating the gene expression using systems of differential equations and chemical kinetics
_by Anna Chechenina, BII, July-September 2023_

The goal of the project is to build modelling setup for robust simulations of bulk gene expression patterns given initial parameter values. In order to build such model we designed a set of chemical reactions that describes events leading to gene expression activation or repression. By adjusting rate constants and initial concentration of the species in the system, one can get simulation of number of transcripts over time. <br>

## One-equation models
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

For both models we performed simulations and parameter space exploration. However, neith 
You can find full code with runs and comments in the notebook `ODE_modellng_part1.ipynb`. 

In order to better visualise the parameter values we made the Ploty application where the user can select initial parameters for the model, run the simulation and check the trajectory. Follow the instruction below to download and run the script:

#### Interactive parameter exploration. Installation

To get the scripts clone the git repository:

```bash
git clone https://github.com/checheanya/DEs_for_expr_modelling.git && cd expression_modelling/part1
```

Create a `conda/mamba` environment with necessary packages and activate it:

```bash
conda env create -f environment.yml
conda activate expr_modelling
```

### Interactive parameter exploration. Usage

To run the script, just call it from the directory where the tool is located:

```bash
python lineplots_2equations.py 
```
for plotting the line plots and

```bash
python scatterplots_3equations.py 
```
for plotting the scatterplots. Then proceed to the link that appears on the screen.

## Mass action based modelling 

Since the model described above was not giving us results close to the real-world data, we decided to build another more complex model taking into account both chromatin state and regulators. In this model we assume four set of reactions leading to the gene expression activation or silencing. First set of reactions corresponds to the chromatin state. For both regulatory regions and promoters we are assuming two dynamic states: open and closed. Closed state correspond to the hetero-chromatin with nucleosomes occupying the region and so no transcripts can be produced by the polymerase. Open state, on the contrary, corresponds to the euchromatin, when promotors or regulatory elements are located in the inter nucleosomic regions accesible fo the polymerase.

Full set of discussed reactions and equations can be found in the file `part2/Reactions_and_equations.pdf`.

Simulation experiments were performed using Python library [GillesPy2](https://github.com/StochSS/GillesPy2) that allows to construct systems of chemical kinetics reactions and converts it the system of differential equations for each species. 

For full code, simulations and comments, please consult the notebook in the `part2` section. This notebook is quite extensive and covers all analysis from building the model to parameter exploration and actual simulations. You can find the contents in provided notebook helpful for navigating over the document. Since the notebook size is large, you might prefer to view the notebook in [Google Colab](https://drive.google.com/file/d/14BWNZXek6hNR2mW3YLtCV6a9weFNKZT_/view?usp=sharing).
Additionally, you can find a parameter exploration plots [here](https://docs.google.com/document/d/e/2PACX-1vSMz3o24yA6p8Afb96XpxYsBYigMyanBkF8m0hYsVXj2-iLRETh1Beg08Hjb31t0uFrdeMLiaRruXnI/pub).

Sourses:
1. Piras, V., Tomita, M. & Selvarajoo, K. Transcriptome-wide Variability in Single Embryonic Development Cells. Sci Rep 4, 7137 (2014). https://doi.org/10.1038/srep07137
