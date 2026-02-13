<img width="100" height="100" alt="logo" src="https://github.com/user-attachments/assets/d6bc24e3-fb68-4d8a-9cab-8cc4149ad2b4" />

# eVCGsampler R package

eVCGsampler provides a principled framework for sampling of VCG (Virtual Control Group), using energy distance-based covariate balancing. 
The package includes visualization tools for assessing covariate balance, as well as a permutation test to evaluate the statistical significance of the deviations.

## Example.

Test for 3 covariates before balancing, comparison of the treated groups (TG) with the data pool (POOL), shows high imbalance.

Distance permutation test TG vs POOL:

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/29ed667e-06db-4f4f-aa1c-ab4d90ceebed" />


By running the function:  VCG_sampler(treated ~ cov1 + cov2 + cov3, data=dat, n=10)

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/0582bb36-a56f-4c01-a136-163413bce578" />

Distance permutation test TG vs VCG:

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/7a652740-9d37-4191-9ccb-db665d6c62ec" />

Plot specifically for the variable cov3: plot_var(dat_out, what='cov3’)

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/b809aebc-d61f-44a6-ad44-64fc5444ea29" />

## Best VCG size (exploratory)

With BestVCGsize(treat ~ cov1 + cov2 + cov3, data=dat), you can explore the best size for VCG with the best balance of covariates.
It may not necessarily be the best size in terms of power or validity of the study.

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/417ff769-3db2-4c8c-8080-6a8ce2832ab1" />

## Multiple VCG samples

If multiple VCG samples are required, use: multiSampler(treat~cov1+cov2+cov3, n=10, Nsamples=10, data=dat)

Overview of sample overlapping:

<img width="500" height="500" alt="image" src="https://github.com/user-attachments/assets/3ceefe1f-5be4-457b-a1f4-27fd3432ee7c" />






