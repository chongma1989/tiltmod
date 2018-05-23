# R Package tiltmod: Exponential tilting model for reproducible screening of large scale testing problem

The R package tiltmod develops an exponential tilting method to tilt a mixture of two component of Beta distributions for p-values (left tail areas) obtained from a large scale hypotheses, yielding a tilted mixture model by conditioning on the false discovery rate. The package includes the function UBMM (MBM) to fit a mixture of Uniform-Beta (two general Beta distributions) for p-values (left tail areas) of test statistics for a large amount of hypotheses by using a Boosted EM algorithm. The Boosted EM algorithm shares similar ideas from the backtracking line search algorithm, which greatly boosts the EM algorithm. Moreover, the UBMM (MBM) functions are built by using C++, which is quite fast and stable. The package also includes some utility functions to generate the tilted false discovery rates and frequency network. 

## Install tiltmod
Download the R package tiltmod via 
