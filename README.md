# RQA using sampling approach RQA_Samp

[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/pucicu/RQA_Samp/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/pucicu/RQA_Samp)
![file size](https://img.shields.io/github/repo-size/pucicu/RQA_Samp)
![GitHub Release](https://img.shields.io/github/v/release/pucicu/RQA_Samp)

## General

Julia function to approximate some of the measures of the recurrence quantification
analysis (RQA). The calculation is performed in an efficient way by sampling
a limited number of recurrence plot lines directly on the time series
(and without preceeding calculation of the recurrence plot).

## Files

- `RPLineLengths.jl`: Julia functions for estimating line length distibutions and some RQA measures
- `example_calculation.jl`: example script demonstrating the use of the sampling scheme and how to get the RQA measures

## Functions

- `get_hist_diagonal_RP(x, ε)`: Compute the histogram of diagonal line lengths in a recurrence plot.
- `get_hist_diagonal_woRP(x, ε)`: Compute the histogram of diagonal line lengths directly from a time series without 
explicitly building the recurrence plot.
- `get_hist_diagonal_sampled(x, ε, M)`: Compute an approximate histogram of diagonal line lengths by random sampling, reducing computation cost compared to a full scan.
- `rqa(P, N)`: Compute Recurrence Quantification Analysis (RQA) measures from a diagonal line length histogram.

with the arguments
- `x`: timeseries or a phase space vector (needs to be prepared first, e.g. using the function `embed` from the `DynamicalSystems.jl` package)
- `ε`: the recurrence threshold
- `M`: number of samplings of line structures
- `P`: histogram of diagonal lines as returned by the `get_hist_diagonal` functions
- `N`: number of points in the upper (or lower) triangle of the recurrence matrix (number of searches using the sampling approach).

## Usage

```
histL, N_rp = get_hist_diagonal_sampled(x, ε, M);
r = rqa(histL, N_rp)
```

## Copyright

Norbert Marwan\
Potsdam Institute for Climate Impact Research\
3/2026

License: GPLv3
