# bayrhm
 
The `bayrhm` package implements 3 different bayesian regularization priors inspired by [heuclin et al.(2023)](https://institut-agro-montpellier.hal.science/hal-04238536/) on a multiple random effects model used in genetic cartography called RHM, named by  [nagamine et al.(2012)](https://f1000research.com/assets/download/1092003)

## RHM model :

The RHM model we use is as follow :

$$
y = \sum_l \lambda_l u_l + \epsilon \quad \quad u_l \sim \mathcal N(0,A_l)\quad \quad \epsilon \sim \mathcal N(0,\sigma^2)
$$

With $y$ the phenotype, $u_l$ the random effect associated with a similarity matrix $A_l$ computed on the genetic region indexed by $l$, and $\epsilon$ the residual error.

## Regularizations :

All the regularization priors are used to shrink non-pertinent $\lambda_l$ coefficients to zero.

+ The horseshoe prior doesn't assume similarity between adjacent genetic region and shrink all coefficients the same way.

$$
\lambda_l \sim \mathcal N(0,\sigma^2\tau^2W_l^2) \quad \quad \tau^2 \sim \mathcal C^+(0,1) \quad \quad w_l \sim \mathcal C^+(0,1)
$$

+ The fusion-horseshoe prior assume strong similarity between adjacent genetic regions and shrink only the regions at the end and beginning of chromosomes.

$$
\lambda_l - \lambda_{l-1} \sim \mathcal N(0,\sigma^2\tau^2\omega_l^2) \quad \quad \omega_l \sim \mathcal C^+(0,1)
$$

$$
\lambda_1 \sim \mathcal N(0,\sigma^2w_1^2) \quad \quad \lambda_L \sim \mathcal N(0,\sigma^2w_2^2) \quad \quad w_1,w_2 \sim \mathcal C^+(0,1)
$$

+ The fused-horseshoe prior is a mix between the first and second prior , assuming similarity between adjacent regions, and shrinking all coefficients the same way.

$$
\lambda_l - \lambda_{l-1} \sim \mathcal N(0,\sigma^2\tau^2\omega_l^2) \quad \quad \omega_l \sim \mathcal C^+(0,1)
$$

$$
\lambda_l|\lambda_{l-1}-\lambda_l \sim \mathcal N(0,\sigma^2w_l^2) \quad \quad w_l \sim \mathcal C^+(0,1)
$$

## Installation

Use `devtools::install_github("nicolasleroy1/bayrhm")` in R console.
