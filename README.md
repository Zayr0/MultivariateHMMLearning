# MultivariateHMMLearning
A novel Multivariate Hidden Markov Model learning method using Coupled Canonical Polyadic Decomposition


# Main Contribution
The main contribution of this work, Coupled CPD, can be found as a function in \TensorDecompositions\CPD_Multivariate_Coupled\multivariate_coupled_cpd_2.m .


# Results
This repo contains a novel method for learning Hidden Markov Models (HMMs) from data using Coupled Canonical Polyadic Decomposition (CPD). Some examples of learning HMMs from synthetic data and real data are given in the examples folder. Furthermore how this method performs with respect to Baum-Welch and Uncoupled CPD can be observed in the mlx files in \PerformanceTests\ (Note that these take long to run).

## Performance Tests
MethodComparison.mlx show the general performance over 100 HMM runs for Baum-Welch, Uncoupled CPD and Coupled CPD, over HMM error and time cost.
Time_Performance_Test.mlx shows how all methods respond to different numbers of samples.
State_Emission_Sequences_Performance.mlx shows how all methods respond to different numbers of states, emissions, and channels.
HyperParameter_Performance.mlx shows how coupled CPD responds to changes in the learning parameters $\alpha$ and $\rho$.

## Examples
Real_Data.mlx shows the application of Coupled CPD on the sleep physionet dataset for a single patient.
NoiseFlatteningNormal.mlx shows the effect of noise on a distribution of data, it will flatten out.
Continuous_Data.mlx shows the effect of noise on coupld CPDs ability to learn HMMs from synthetic data.

# Data
The data is gathered from the Sleep Physionet dataset (https://www.physionet.org/content/sleep-edfx/1.0.0/), and requires the following directory and file format for patient 1, night 1:

"Data\sleep-cassette\SC4001E0-PSG.edf"
"Data\sleep-cassette\SC4001EC-Hypnogram.edf"

Note that here only the sleep-casette part of the dataset is used.

# Other functions

## Error Metrics
Contains some error metrics to compare HMMs.

## Probability Tensors
Contains functions for making and working with Joint Probability Tensors, and data.

### TensorOperations
Contains some elementary functions for working with tensors.

## Utils
### HMM Utils
Contains functions for the creation and permutation of both continuous HMMs and discrete HMMs. It also matches permutations and is capable of generating data sequences from HMMs.

### Utils
Has some utility functions with respect to HMMs, data and distributions

### Plotting
Contains functions for plotting HMMs.

### Other
Contains external packages used.
