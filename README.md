<h1 align="center" style="font-weight: bold;">scCOSMiX</h1>

<p align="center">
  <align="center">A <ins>MiX</ins>ed-Effects Framework for Differential <ins>CO</ins>expre<ins>S</ins>sion and Transcriptional Interactions Modeling in <ins>sc</ins>RNA-Seq
</p>

## The Model

Let $Y_1, Y_2$ represent the expression levels of two genes in $n$ cells. Let $X$ and $Z$ be design matrices for fixed effects and random effects, respectively. In the scCOSMiX framework, the joint expression of $Y_{i1}, Y_{i2}$ is modeled using a zero-inflated bivariate Gaussian copula with negative binomial marginals. Every parameter in the model can be modeled as a linear combination of relevant fixed and random effects via an appropriate one-to-one transformation.

The joint distribution of $Y_1, Y_2$ is given by 

$$
\begin{aligned}
&f(y_{1}, y_{2} ; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho, p_{1}, p_{2}) =\\
&(1-p_{1})(1-p_{2})g(y_{1}, y_{2}; \  \mu_1, \sigma_1, \mu_2, \sigma_2, \rho) + \\
&\begin{cases} 
    0 & \text{if } y_{1}> 0 \text{, } y_{2}> 0 \\
    (1-p_{1})p_{2}f_{NB}(y_{1}; \ \mu_1, \sigma_1)  & \text{if } y_{1}> 0 \text{, } y_{2}= 0 \\
     p_{1}(1-p_{2})f_{NB}(y_{2}; \ \mu_2, \sigma_2) & \text{if } y_{1}= 0 \text{, } y_{2}> 0 \\
     p_{1}p_{2} + (1-p_{1})p_{2}f_{NB}(y_{1}; \ \mu_1, \sigma_1) + p_{1}(1-p_{2})f_{NB}(y_{2}; \ \mu_2, \sigma_2) & \text{if } y_{1}= 0 \text{, } y_{2}=  0   
  \end{cases}
  \end{aligned}
$$

where $g$ is the distribution function from a bivariate copula with association parameter $\rho$, and $p_1$, $p_2$ are parameters representing the probability of an observation from their respective marginal being zeroed out by a drop-out event. Parameters can be made dependent on covariates in the following way 

$$
\begin{aligned}
\log(\mu_{i1}) &= X_{i,\mu_1}\beta_1 + Z_{i,\mu_1}u_{\mu_1} + \log\left(S_i\right)\\
\log(\mu_{i2}) &= X_{i,\mu_2}\beta_2 + Z_{i,\mu_2}u_{\mu_2} + \log\left(S_i\right)\\
\log(\sigma_{i1}) &= X_{i,\sigma_1}\alpha_1 + Z_{i,\sigma_1}u_{\sigma_1}\\
\log(\sigma_{i2}) &= X_{i,\sigma_2}\alpha_2 + Z_{i,\sigma_2}u_{\sigma_2}\\
atanh(\rho_i) &= X_{i,\rho}\tau + Z_{i,\rho}u_{\rho}\\
logit(p_{i1}) &= X_{i,p_1}\kappa_1 + Z_{i,p_1}u_{p_1}\\
logit(p_{i2}) &= X_{i,p_2}\kappa_2 + Z_{i,p_2}u_{p_2},
\end{aligned}
$$

where the $\log, \text{atanh},$ and $\text{logit}$ transformations on the left-hand side ensure that each parameter remains within its valid parameter space. $X_{\theta}$ and $Z_{\theta}$ are submatrices of the full design matrices $X$ and $Z$, containing only the columns relevant to parameter $\theta$. $\log(S_i)$ term is an offset equal to the logarithm of the total sequencing depth of cell $i$.


## Fitting Function

The fitting function for the scCOSMiX model with negative binomial marginals is called `sccosmix`. 

The function works as follows:

```{r}
source(".../functions/scCOSMiX_functions.R")

sccosmix(formula, data, offset1, offset2)
```
Parameters:
* `formula`: A list of seven formula objects specifying the covariate-dependence of the different parameters.
  * `[[1]]`: Covariate-dependence of $\mu_1$
  * `[[2]]`: Covariate-dependence of $\mu_2$
  * `[[3]]`: Covariate-dependence of $\sigma_1$
  * `[[4]]`: Covariate-dependence of $\sigma_2$
  * `[[5]]`: Covariate-dependence of $\rho$
  * `[[6]]`: Covariate-dependence of $p_1$
  * `[[7]]`: Covariate-dependence of $p_2$
* `data`: A data.frame whose column names correspond to the variables referenced in the formula list.
* `offset1`, `offset2`: Vectors containing the $\log(S_i)$ offset terms.



The function `summary.custom()` outputs a table containing estimates, standard errors, and p-values for the model's fixed effects.

## Simulation Studies

We conducted simulation studies on the power, coverage, robustness, FDR control, and precision-recall of the scCOSMiX method with negative binomial marginals and a Gaussian copula. These studies can be found in the `simulations` folders of this repository.

## Real Data Analysis

We performed differential co-expression analysis on triple-negative breast cancer data from gene expression omnibus GSE266919 and colorectal cancer data from gene expression omnibus GSE108989. Code for the analysis can be found in the `real data analysis` folders of this repository.





