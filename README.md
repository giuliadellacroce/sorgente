# Code for 'Estimating the contagiousness ratio between two viral strains' 
_by Della Croce di Dojola G., Mastrantonio G. Cerutti F., Ghisetti V., Gasparini M., Bibbona E._

#### FILES IN THE REPOSITORY

_code.R_: the main file containing the R function `sorgente`, which performs both the deterministic and the stochastic method described in the paper

_workspace_dati_piemonte.RData_: R workspace with the data from the _Piemonte_ region, whose results are described in the paper

_dataset_det_sim.RData_: R workspace containing the deterministic simulated dataset, whose results are in the Supplementary Material

_dataset_tauleap_k5_eps2.RData_: R workspace containing the first stochastic simulated dataset, whose results are in the Supplementary Material

_dataset_I032.RData_: R workspace containing the second stochastic simulated dataset, whose results are in the Supplementary Material


#### FUNCTION `sorgente`

The R function `sorgente` is mainly divided in two parts:

**Deterministic method** (argument `stoc`=0), in which the deterministic method is implemented, i.e. maximum likelihood estimation of the parameters and computation 
of both profile likelihood and parametric bootstrap 95% C.I.

**Stochastic method** (argument `stoc`=1), in which the stochastic method is implemented, i.e. the MCMC algorithm.

###### Input arguments

`stoc` =0 to apply deterministic method, = 1 to apply stochastic method

`kstart` starting value of parameter $k$ in the maximisation of the likelihood (_det. method_)
 
`epsstart` starting value of parameter $I_0^2$ (often referred to as _epsilon_) in the maximisation of the likelihood (_det. method_)

`stepk` starting gridwith to compute profile C.I. for parameter $k$ (_det. method_)

`stepeps` starting gridwith to compute profile C.I. for parameter $I_0^2$ (_det. method_)

`max_k` maximum value of $k$ in the computation of the profile likelihood C.I. for parameter $I_0^2$ (_det. method_)

`x_axis` vector containing the days in which sequencing has been performed (i.e. $n_j \neq 0$) for the considered dataset (_det. method_)

`incr_k` used to define the starting point for parameter $k$ in the likelihood maximisation when applying parametric bootstrap (_det. method_)

`incr_eps` used to define the starting point for parameter $I_0^2$ in the likelihood maximisation when applying parametric bootstrap (_det. method_)

`num_it` total number of iterations of the MCMC (_stoc. method_)

`Burnin` burn-in (_stoc. method_)
  
`Thin` thin (_stoc. method_)

`Nn` vector of increments used in the uniform proposals (MH) of the $Y_j^2$ terms (_stoc. method_)

`Nrho` vector of increments used in the uniform proposals (MH) of the $I_j^2$ terms (_stoc. method_)

`N_eps` increment used in the uniform proposal (MH) of parameter $I_0^2$ (_stoc. method_)

`kstart_MCMC` starting value for parameter $k$ in the MCMC (_stoc. method_)

`epsstart_MCMC` starting value for parameter $I_0^2$ in the MCMC (_stoc. method_)

###### Output arguments

A matrix called _sol_, structured as follows

If the **deterministic method** is applied: the first row of _sol_ contains the ML estimate, the 95% profile C.I. and the parametric bootstrap 95% C.I. 
for parameter $k$; the second row contains the same for parameter $I_0^2$ (_sol_ is a matrix 2*5).

If the **stochastic method** is applied: _sol_ will have a number of rows equal to $N$ (number of iterations stored, i.e. burn-in and thin already applied);
the first column contains the samples of parameter $k$, columns from $2$ to $n+1$ (with $n$=number of days) contain the samples of the $Y_j^2$ terms 
(we are not interested in column $2$), columns from $n+2$ to $2n+1$ contain the samples of the $I_j^2$ terms, in particular column $n+2$ contains the samples of
parameter $I_0^2$ (we are not interested in column $2n+1$).

#### Reproducibility of the results

We outline, for each available dataset, the input arguments passed to the function `sorgente`, in order to obtain the results shown in the paper and in the
supplementary material.

###### Real data from the Piemonte region (file _workspace_dati_piemonte.RData_)

`stoc` = $0/1$, `kstart` = $3$, `epsstart` = $15$, `stepk` = $0.1$, `stepeps` = $0.1$, `max_k` = $10$, `x_axis` = $c(2,16,30,44,58)$, `incr_k` = $0.2$, `incr_eps` = $-3$, `num_it` = $6 \cdot 10^7$, `Burnin` = $3 \cdot 10^6$, `Thin` = $10^4$, `Nn` = $c(rep(10,16),rep(20,3),rep(40,10),rep(70,5),rep(100,10),rep(150,14))$,
`Nrho` = $c(rep(20,19),rep(50,10),rep(130,11),rep(150,18))$, `N_eps` = $5$, `kstart_MCMC` = $2$, `epsstart_MCMC` = $10$.

###### Deterministic simulated dataset (file _dataset_det_sim.RData_)

`stoc` = $0$, `kstart` = $2.5$, `epsstart` = $10$, `stepk` = $0.1$, `stepeps` = $0.1$, `max_k` = $3$, `x_axis` = $c(2,8,16,24,30)$, `incr_k` = $0.2$, `incr_eps` = $-3$, `num_it` = -, `Burnin` = -, `Thin` = -, `Nn` = -, `Nrho` = -, `N_eps` = -, `kstart_MCMC` = -, `epsstart_MCMC` = -.

###### First stochastic simulated dataset (file _dataset_tauleap_k5_eps2.RData_)

`stoc` = $0/1$, `kstart` = $3$, `epsstart` = $10$, `stepk` = $0.1$, `stepeps` = $0.1$, `max_k` = $10$, `x_axis` = $c(2,9,16,23,29)$, `incr_k` = $-0.2$, `incr_eps` = $3$, `num_it` = $3 \cdot 10^7$, `Burnin` = $3 \cdot 10^6$, `Thin` = $3 \cdot 10^3$, `Nn` = $c(rep(20,5),rep(30,4),rep(50,7),rep(100,3),rep(150,10))$,
`Nrho` = $c(rep(20,5),rep(30,4),rep(50,7),rep(100,3),rep(150,10))$, `N_eps` = $20$, `kstart_MCMC` = $2$, `epsstart_MCMC` = $10$.

###### Second stochastic simulated dataset (file _dataset_I032.RData_)

`stoc` = $0/1$, `kstart` = $2.5$, `epsstart` = $40$, `stepk` = $0.1$, `stepeps` = $0.1$, `max_k` = $10$, `x_axis` = $c(2,13,24,35,46)$, `incr_k` = $-0.2$, `incr_eps` = $4$, `num_it` = $8 \cdot 10^7$, `Burnin` = $3 \cdot 10^6$ (additional burn-in is needed a-posteriori), `Thin` = $1.5 \cdot 10^4$, 
`Nn` = $c(rep(10,16),rep(20,3),rep(40,10),rep(70,5),rep(100,10),rep(150,2))$,
`Nrho` = $c(rep(20,19),rep(50,10),rep(130,11),rep(150,6))$, `N_eps` = $10$, `kstart_MCMC` = $2$, `epsstart_MCMC` = $10$.





