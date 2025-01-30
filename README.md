# bayrhm
 
The `bayrhm` package implements 3 different bayesian regularization priors inspired by [heuclin et al.(2023)](https://institut-agro-montpellier.hal.science/hal-04238536/) on a multiple random effects model used in genetic cartography called RHM, named by  [nagamine et al.(2012)](https://f1000research.com/assets/download/1092003)

## Installation

Use `devtools::install_github("nicolasleroy1/bayrhm")` in R console.

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
\lambda_l\mid\sigma^2,\tau^2,W^2 \sim \mathcal N(0,\sigma^2\tau^2W_l^2) \quad \quad \tau^2 \sim \mathcal C^+(0,1) \quad \quad w_l^2 \sim \mathcal C^+(0,1)
$$

+ The fusion-horseshoe prior assume strong similarity between adjacent genetic regions and shrink only the regions at the end and beginning of chromosomes.

$$
\lambda_l - \lambda_{l-1} \mid \sigma^2,\tau^2,\omega_l^2 \sim \mathcal N(0,\sigma^2\tau^2\omega_l^2) \quad \quad \omega_l^2 \sim \mathcal C^+(0,1)
$$

$$
\lambda_1\mid\sigma^2,w_1^2 \sim \mathcal N(0,\sigma^2w_1^2) \quad \quad \lambda_L\mid\sigma^2,w_2^2 \sim \mathcal N(0,\sigma^2w_2^2) \quad \quad w_1^2,w_2^2 \sim \mathcal C^+(0,1)
$$

+ The fused-horseshoe prior is a mix between the first and second prior , assuming similarity between adjacent regions, and shrinking all coefficients the same way.

$$
\lambda_l - \lambda_{l-1}\mid\sigma^2\tau^2\omega_l^2 \sim \mathcal N(0,\sigma^2\tau^2\omega_l^2) \quad \quad \omega_l^2 \sim \mathcal C^+(0,1)
$$

$$
\lambda_l\mid\lambda_{l-1}-\lambda_l,\sigma^2,w_l^2 \sim \mathcal N(0,\sigma^2w_l^2) \quad \quad w_l^2 \sim \mathcal C^+(0,1)
$$


## Cauchy prior as a mixture of inverse-gammas :

We use following mixture of inverse-gammas to sample the $\tau^2$ hyperparameter of the horseshoe prior :

$$
\tau^2 \sim \mathcal C^+(0,1) \quad \quad \tau\mid\xi = \frac{1}{\xi^2} \quad \quad \xi \sim \mathcal IG(\frac{1}{2},\frac{1}{2})
$$

## Marginal posterior distributions :

### Horseshoe regularization :

The global marginal posterior opposite log density is :

$$
U(y\mid *) + \sum_l U(u_l) + \sum_l U(\lambda_l\mid\sigma^2,w_l^2,\tau^2) + U(\sigma^2) + \sum_l U(w_l^2) +U(\tau^2)
$$

$$
U(\lambda,\tau^2,(u_l),\sigma^2,(w_l)) = \frac{\|y - \sum_l \lambda_l u_l\|^2}{\sigma^2} + \frac{\|u_l\|^2_{A_l^{-1}}} + \frac{\|\lambda\|^2_{W_l^2}}{\tau^2\sigma^2} + \frac{L+n}\log(\sigma^2) + \sum_l \log(w_l^2) + L\log(\tau^2) + \log(1+\tau^2) + \sum_l \log(1+w_l^2)
$$



