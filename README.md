# bayrhm
 
The `bayrhm` package implements 3 different bayesian regularization of a multiple random effects model for genetic cartography.

## RHM model :

The standard RHM model is as follow :

$$
y = \sum_l \lambda_l u_l + \epsilon \quad \quad u_l \sim \mathcal N(0,A_l)\quad \quad \epsilon \sim \mathcal N(0,\sigma^2)
$$

With $y$ the phenotype, $u_l$ the random effect associated with a similarity matrix $A_l$ computed on the genetic region indexed by $l$, and $\epsilon$ the residual error.

## Regularizations :

All the regularizations prior are used to shrink non-pertinent $lambda_l$ coefficients to zero.

+ The horseshoe prior doesn't assume similarity between adjacent genetic region and shrink all coefficients the same way.
+ The fusion-horseshoe prio assume strong similarity between adjacent genetic region and shrink only the regions situated at the end and beginning of chromosomes.
+ The fused-horseshoe prior is a mix between the first and second prior , assuming similarity between adjacent regions, and shrinking all coefficients the same way.


