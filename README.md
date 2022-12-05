# Clustering Simulations

This repository contains the scripts needed to simulate outbreaks and perform clustering based on cov2clusters (Sobkowiak et.al. 2022 https://doi.org/10.1186/s12864-022-08936-4).

The clustering function is based on an altered version of simulateoutbreaks.R from the package SEEDY (Version 1.3) (Worby et.al. 2022 https://doi.org/10.1371/journal.pone.0129745) by Pouya Haghmaram here: https://github.com/Pouya-Haghmaram/Clustering-outbreaks.

The main changes from this code are:

1. The inside and outside outbreaks can have different infection rate (inf.rate.in and inf.rate.out) and different mutation rates (mut.prob.site.in and mut.prob.site.out).

2. Either a lag time in days, or a minimum proportion of outside susceptibles infected, can be specified before any infection can be introduced to the inside population.

3. Inside of randomly chosing the source of a new infection from all currently infected individuals, the source of infection will be chosen from the pool of currently infected individuals that are closest to a putative infection time, which is drawn from a gamma distribution with a specified mean (mean.infect).



