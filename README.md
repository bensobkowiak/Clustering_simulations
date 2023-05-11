# Clustering Simulations

This repository contains the scripts needed to simulate outbreaks and perform clustering used in the paper Sobkowiak et. al. The utility of genomic sequence data for informative clustering under different outbreak epidemiological scenarios and sampling strategies (2023). In preparation.

The outbreak simulation function is based on an altered version of simulateoutbreaks.R from the package SEEDY (Version 1.3) (Worby et.al. 2022 https://doi.org/10.1371/journal.pone.0129745) by Pouya Haghmaram, documented here: https://github.com/Pouya-Haghmaram/Clustering-outbreaks. The clustering used is based on cov2clusters (Sobkowiak et.al. 2022 https://doi.org/10.1186/s12864-022-08936-4).

The main changes from this code are:

1. The inside and outside outbreaks can have different infection rate (inf.rate.in and inf.rate.out) and different mutation rates (mut.prob.site.in and mut.prob.site.out).

2. Either a lag time in days (time.lag.outside), or a minimum percentage of outside susceptibles infected (min.perc.outside.inf), can be specified before any infection can be introduced to the inside population. N.B. both options can be null for no lag but both options can not be specified together.

3. Instead of randomly chosing the source of a new infection from all currently infected individuals, the source of infection will be chosen from the pool of currently infected individuals that are closest to an estimated infection time of the newly infected host, which is drawn from a gamma distribution with a given mean and rate (mean.infect and rate.infect).



