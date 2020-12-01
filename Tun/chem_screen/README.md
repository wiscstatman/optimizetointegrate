This repo contains standard Bayesian and random-weighting analysis for a chemical screening data. 

The publicly-available dataset is hosted on [Zenodo](https://zenodo.org/record/1411506#.X8bP881KhPY).

Rationale of the experiment was explained in Tony Gitter's [JCIM paper](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00363). Briefly, PriA-SSB interaction is imperative to bacterial DNA replication. We want to virtually screen potential chemical compounds that could effectively target PriA-SSB interaction. 

Response = \% inhibition of PriA-SSB.

Predictors consist of chemical substructures of a drug represented by binary entries of length 1024 via Morgan-fingerprinting. So, dimension p = 1024.

Training set sample size = 72,423.
Test set sample size = 22,434.