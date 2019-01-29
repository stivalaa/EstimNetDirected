# EstimNetDirected

EstimNetDirected implements the Equilibrium Expectation (EE) algorithm for estimating parameters of exponential random graph models (ERGMs) for large directed networks.

The source code is written in C and (optionally) uses MPI to run multiple estimations in parallel. It was written on a Linux cluster system with OpenMPI but should be portable to any system with a standard C compiler (it works on cygwin under Windows for example). It uses the Random123 library for random number generation. 

Also included is a simple Python demonstration implementation. The Python implementation uses the NumPy library for vector and matrix data types and functions. In addition, there are R scripts for estimating standard errors and plotting results from the output.

For more information, and the original implementation for undirected networks, see http://www.estimnet.org/.

As well as some simulated networks, as an empirical example, the political bloggers network from the paper:

Adamic, Lada A, & Glance, Natalie. (2005). The political blogosphere and the 2004 US election: divided they blog. Pages 36-43 of: Proceedings of the 3rd international workshop on link discovery. ACM.

is included. This network was downloaded from Mark Newman's network data page.

If you use this software (or any alternative implementation of the algorithms described in the references), please cite the papers below in any resulting publications.

## References

Borisenko, A., Byshkin, M., & Lomi, A. (2019). A Simple Algorithm for Scalable Monte Carlo Inference. arXiv preprint arXiv:1901.00533. https://arxiv.org/abs/1901.00533

Byshkin, M., Stivala, A., Mira, A., Robins, G., & Lomi, A. (2018). [Fast Maximum Likelihood Estimation via Equilibrium Expectation for Large Network Data](https://www.nature.com/articles/s41598-018-29725-8). *Scientific Reports* 8:11509. https://doi.org/10.1038/s41598-018-29725-8

Byshkin, M., Stivala, A., Mira, A., Krause, R., Robins, G., & Lomi, A. (2016). Auxiliary parameter MCMC for exponential random graph models. *Journal of Statistical Physics*, 165(4), 740-754. https://doi.org/10.1007/s10955-016-1650-5

