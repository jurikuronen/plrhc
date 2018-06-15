# plrhc

Published MATLAB code for the PLRHC-BIC algorithm for learning pairwise Markov network structures with logistic regression.

The logistic regression computations rely on the software package minFunc by Mark Schmidt (2005) (https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html). The complete set of .m and .mex files for the current version (2012) of minFunc are included in this repository. In addition, two files from the L1General package by the same author and two custom loss functions for weighted data are included in the custom folder within.

## How to use

### Data set format

Your data set should be a Nxd matrix in 'double' format, where N is the number of data samples (rows) and d is the number of variables (columns). The outcome space must be {1, 2, ...}. The code supports polytomous variables, but the method has been rigorously tested on only binary variables. Use `test_data.m` to see if your data is in the right format.

### Compile MEX files

Depending on your system, you may not have to do this.

In MATLAB, compile the MEX function from the C++ file `mex_compute_weighted_data.cpp` with:
```
mex CXXFLAGS="\$CXXFLAGS -std=c++0x" mex_compute_weighted_data.cpp -largeArrayDims
```

In MATLAB, compile the MEX functions from the minFunc package with:
```
cd minFunc_2012
mexAll
```

### Learn the Markov network structure from data.

To run the full PLRHC-BIC with gamma = 0.5:
1) Set the following in MATLAB:
```
g = 0.5; % Gamma value
penalty = [0 0]; % Don't use Lq penalty
parallel = 1; % Use MATLAB's parallelization
depth = 3; % Search spaces will contain nodes within this graph-theoretic distance with respect to E^star.
```
2) Run the following:
```
Gstar = plr(data, g, penalty, parallel); % PLR
Ghator = hcor(data, Gstar, g, penalty, depth, parallel); % Solution given by PLRHC_OR
Ghat = hc(data, Ghator, g, penalty, 0); % Solution given by PLRHC-BIC
```

