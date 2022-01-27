# EstimNetDirected

EstimNetDirected implements the Equilibrium Expectation (EE) algorithm for estimating parameters of exponential random graph models (ERGMs) for large directed (or undirected) networks.

The source code is written in C and (optionally) uses MPI to run multiple estimations in parallel. It was written on a Linux cluster system with OpenMPI but should be portable to any system with a standard C compiler (it works on cygwin under Windows for example). It uses the Random123 library for random number generation. 

Also included is a simple Python demonstration implementation. The Python implementation uses the NumPy library for vector and matrix data types and functions. In addition, there are R scripts for estimating standard errors and plotting results from the output.

For more information, and the original implementation for undirected networks, see http://www.estimnet.org/.

As well as some simulated networks, as empirical examples, the political bloggers network from the paper:

Adamic, Lada A, & Glance, Natalie. (2005). The political blogosphere and the 2004 US election: divided they blog. Pages 36-43 of: Proceedings of the 3rd international workshop on link discovery. ACM.

and the network science coauthorship network from the paper:

Newman, M. E. (2006). Finding community structure in networks using the eigenvectors of matrices. Physical review E, 74(3), 036104.

are included. These networks were downloaded from Mark Newman's network data page.

If you use this software (or any alternative implementation of the algorithms described in the references), please cite the papers below (specifically Stivala, Robins, & Lomi (2020) for the EstimNetDirected software, and Byshkin et al. (2018) for the EE algorithm) in any resulting publications.

The "citation ERGM" (cERGM) model variant, which can also be estimated with this software, is described by:

Schmid, C., Chen, T., & Desmarais, B. (2021). Generative Dynamics of Supreme Court Citations: Analysis with a New Statistical Network Model. Political Analysis, 1-20. doi:10.1017/pan.2021.20

The original statnet (http://statnet.org/) R implementation of cERGM is available from https://github.com/schmid86/cERGM/.


## Funding

Development of the EstimNetDirected software was funded by the Swiss National Science Foundation project numbers 167326 (NRP 75) and 200778.

## References

Borisenko, A., Byshkin, M., & Lomi, A. (2019). A Simple Algorithm for Scalable Monte Carlo Inference. arXiv preprint arXiv:1901.00533. https://arxiv.org/abs/1901.00533

Byshkin, M., Stivala, A., Mira, A., Krause, R., Robins, G., & Lomi, A. (2016). Auxiliary parameter MCMC for exponential random graph models. *Journal of Statistical Physics*, 165(4), 740-754. https://doi.org/10.1007/s10955-016-1650-5

Byshkin, M., Stivala, A., Mira, A., Robins, G., & Lomi, A. (2018). [Fast Maximum Likelihood Estimation via Equilibrium Expectation for Large Network Data](https://www.nature.com/articles/s41598-018-29725-8). *Scientific Reports* 8:11509. https://doi.org/10.1038/s41598-018-29725-8

Stivala, A., Robins, G., & Lomi, A. (2020). [Exponential random graph model parameter estimation for very large directed networks](https://doi.org/10.1371/journal.pone.0227804). *PloS ONE*, 15(1), e0227804. https://arxiv.org/abs/1904.08063

