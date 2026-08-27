<!-- PROJECT LOGO -->
<br />
<h3 align="center">Spherically Embedded Time Series with Unknown Trend and Periodic Components</h3>

</div>

<!-- ABOUT THE PROJECT -->
## Abstract
Spherically embedded time series with values naturally residing on or can be equivalently mapped to the sphere often exhibit complex non-stationarity driven by latent trend and periodic components. Traditional Euclidean methods fail to account for the sphere's intrinsic non-Euclidean geometry, leaving a critical gap in rigorous modeling and forecasting techniques. To address this methodological gap, we propose a unified geometric framework to analyze nonstationary spherically embedded time series. Central to our approach is a novel nonparametric spherical trend-periodicity decomposition model that uses an optimal-transport-based removal operation to sequentially extract the smooth trend and periodic components while preserving spherical topology. The resulting de-trended and de-seasonalized stationary residuals can be further modeled using a spherical autoregressive model, formalizing a novel trend-periodic spherical autoregressive model. Theoretical foundations for the modeling procedure are established on the consistency under temporal dependence. Extensive simulations corroborate these theoretical guarantees and demonstrate the superior finite-sample predictive performance of the trend-periodic spherical autoregressive model compared to existing spherical autoregressive models and neural network models. Finally, we validate the practical utility of our methodology through applications to electricity generation compositions and bike trip volume profiles, yielding significantly enhanced forecasting accuracy while providing interpretable insights into the underlying structural dynamics.

## Setup Instructions for Reproducibility

To ensure full reproducibility of the results presented in the paper, please configure your environment as follows:
### Software Requirements
The analysis was performed using **R version 4.5.1** and **Python 3.13.13**. The following R and Python packages and their specific versions are required to replicate the results:

* **R Packages:** ggplot2 (version 3.5.2), plot3D (version 1.4.1), ggforce (version 0.5.0), patchwork (version 1.3.0), psych (version 2.5.3), osqp (version 0.6.3.3), transport (version 0.15-4), easyCODA (version 0.40.2), tidyr (version 1.3.1), movMF (version 0.2-9), trust (version 0.1-8), parallel (version 4.5.0), foreach (version 1.5.2), doSNOW (version 1.0.20), xtable (version 1.8-4).
* **Python Packages:** Standard library modules for Python 3.13.13 --- copy, os, random, sys, time. Third-party packages --- numpy (version 2.4.4), pandas (version 3.0.2), torch (version 2.12.0), matplotlib (version 3.10.8), scipy (version 1.17.1).



### Code
1. main function/basic_function.R --- basic R functions for numerical studies. It contains codes for the global and local Fréchet regression from <https://github.com/functionaldata/tFrechet>, modified to contain the rolling-window-based bandwidth selection
2. TimePFN --- contain Python code for TimePFN downloaded from <https://github.com/egetaga/TimePFN> 
3. figure/visual_plot.R --- R functions to generate Figures 1, 2, 3, and S1.
4. simulation/prediction.R --- R functions to generate required statistical model results needed for Figures 4, S3, and S4.
5. simulation/python_data_generation.R --- R functions to generate required data for neural network results needed for Figures 4, S3, and S4.
6. simulation/LSTM_sim.ipynb --- Python functions to generate required LSTM model results needed for Figures 4, S3, and S4.
7. simulation/TimePFN_sim.ipynb --- Python functions to generate required TimePFN model results needed for Figures 4, S3, and S4.
8. simulation/prediction_figure.R --- R functions to generate Figures 4, S3, and S4.
9. simulation/estimation.R --- R functions to generate Tables S1, S2, and S3.
10. real data/real_data_engergy.R --- R functions to generate Figures 5, 6, and S2, and to generate required statistical model results needed for Figure S5.
11. real data/real_data_trip.R --- R functions to generate Figures 7, 8, S6, S7, S8, S9, and S10, and to generate required statistical model results needed for Figure S11.
12. real data/LSTM_real.ipynb --- Python functions to generate required LSTM model results needed for Figures S5 and S11.
13. real data/TimePFN_real.ipynb --- Python functions to generate required LSTM model results needed for Figures S5 and S11.
14. real data/prediction_figure.R --- R functions to generate Figures S5 and S11.

### Data
1. real data/dataset/US_energy.csv --- US energy generation composition dataset downloaded from <https://www.eia.gov/electricity/data/browser/>. For details, see Section 6.1 of the paper.
2. real data/dataset/trip.csv--- Preprocessed New York City Citi bike sharing system data with original dataset downloaded from <https://citibikenyc.com/system-data>. For details, see Section 6.2 of the paper.

### Procedure to reproduce Figures 1, 2, 3, and S1
Step 1. Run figure/visual_plot.R to reproduce Figures 1, 2, 3, and S1

### Procedure to reproduce simulation results on prediction (recommend using high performance computer)
Step 1. Run simulation/prediction.R to get required statistical model results needed for Figures 4, S3, and S4.

Step 2. Run simulation/python_data_generation.R to generate required simulated data for neural network results needed for Figures 4, S3, and S4.

Step 3. Run simulation/LSTM_sim.ipynb to generate required LSTM model results needed for Figures 4, S3, and S4.

Step 4. Run simulation/TimePFN_sim.ipynb to generate required TimePFN model results needed for Figures 4, S3, and S4.

Step 5. Run simulation/prediction_figure.R to reproduce Figures 4, S3, and S4.

### Procedure to reproduce simulation results on estimation (recommend using high performance computer)
Step 1. Run simulation/estimation.R to reproduce Tables S1, S2, and S3. 

### Procedure to reproduce real data results
Step 1. Run real data/real_data_engergy.R to reproduce Figures 5, 6, and S2, and to generate required statistical model results needed for Figure S5.

Step 2. Run real data/real_data_trip.R to reproduce Figures 7, 8, S6, S7, S8, S9, and S10, and to generate required statistical model results needed for Figure S11.

Step 3. Run real data/LSTM_real.ipynb to generate required LSTM model results needed for Figures S5 and S11.

Step 4. Run real data/TimePFN_real.ipynb to generate required TimePFN model results needed for Figures S5 and S11.

Step 5. Run real data/prediction_figure.R to reproduce Figures S5 and S11.




