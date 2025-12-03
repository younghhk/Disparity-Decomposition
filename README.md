[![Back to Hub](https://img.shields.io/badge/⬅️%20Back%20to%20Hub-2962FF?style=for-the-badge)](https://github.com/younghhk/NCI)

<a id="dd-proportion"></a>
### Peters–Belson Decomposition of the Proportion of Death

### Motivation

Among women with ER-positive breast cancer, **Black patients experience higher mortality than White patients**. A key question is: **how much of the mortality gap is explained by measured factors** (e.g., stage at diagnosis, access to care, adherence to endocrine therapy) and **how much remains unexplained** by those factors.



### Peters–Belson decomposition (explained vs unexplained)

**Goal.** Partition an outcome gap between groups into:
- an **explained** portion due to differences in covariate distributions, and
- an **unexplained** portion that remains after accounting for those covariates.

### Peters–Belson decomposition on the **proportion of death** (Black vs White)

**Goal.** Split the **risk difference in death** between Black and White patients into:
- an **explained** part due to measured covariates (e.g., stage, access, adherence), and
- an **unexplained** part that remains after accounting for those covariates.

This follows the same logic as the worked examples in the literature (see the PMC article you cited), but focuses explicitly on the **proportion (risk) of death**.



#### Steps
1) Fit a **risk model** for death in the **reference group** (commonly Whites):  
$m(X;\hat\beta_W)=\Pr(Y=1 \mid X,\, G=\text{White})$  
where $Y=1$ indicates death by a chosen time horizon (e.g., 5-year death), and $X$ includes covariates like **stage at diagnosis**, **access**, **adherence**, etc.

2) Apply this White-model to the **Black** covariate distribution to get a *counterfactual Black risk under White coefficients*:

$$
\hat\mu_{B\mid W}=\mathbb{E}\left[m(X_B;\hat\beta_W)\right]     
$$


(Use survey weights for the expectation if the data are from a complex sample.)

3) Compute the observed risks (proportions of death) in each group:

$$
\hat\mu_W=\Pr(Y=1 \mid G=\text{White})
\quad\text{and}\quad
\hat\mu_B=\Pr(Y=1 \mid G=\text{Black})
$$

4) Decompose the **total difference**:

$$
\Delta_{\text{total}}=\hat\mu_B-\hat\mu_W
=\underbrace{\big(\hat\mu_{B\mid W}-\hat\mu_W\big)}_{\Delta_{\text{expl}}\ \text{(explained)}}
+\underbrace{\big(\hat\mu_B-\hat\mu_{B\mid W}\big)}_{\Delta_{\text{unexpl}}\ \text{(unexplained)}}
$$

- **Explained disparity (risk-difference scale):**  
  $\Delta_{\text{expl}}=\hat\mu_{B\mid W}-\hat\mu_W$

- **Unexplained disparity:**  
  $\Delta_{\text{unexpl}}=\hat\mu_B-\hat\mu_{B\mid W}$

- **Percent explained (optional):**  
  $\%\text{ explained}=100\times\dfrac{\Delta_{\text{expl}}}{\Delta_{\text{total}}}$

> This is the Peters–Belson logic applied to **proportions of death**: compare the observed Black risk to the *counterfactual* risk Black patients would have if they faced the **same covariate–outcome relationship** as Whites.


#### Numeric llustration 

Suppose by the chosen horizon you estimate:
- $\hat\mu_W=0.12$ (12% died among Whites)  
- $\hat\mu_B=0.18$ (18% died among Blacks)  
- Counterfactual under White coefficients: $\hat\mu_{B\mid W}=0.15$

Then:
- Total disparity: $\Delta_{\text{total}}=0.18-0.12=0.06$ (6 percentage points)  
- **Explained**: $\Delta_{\text{expl}}=0.15-0.12=0.03$  
- **Unexplained**: $\Delta_{\text{unexpl}}=0.18-0.15=0.03$  
- **Percent explained**: $0.03/0.06=50\%$

Interpretation: **Half** of the Black–White risk difference in death is explained by measured factors (stage, access, adherence, …) and **half** remains unexplained by those covariates.




####  R sketch 

```r
## Peters–Belson on the proportion of death (no Cox)

# Reference model in Whites (logistic for death by horizon)
fit_ref <- glm(death_by5yr ~ stage + access + adherence + x_d,
               data = subset(dat, race == "White"),
               family = binomial())

# Observed proportions
mu_W <- with(subset(dat, race == "White"), mean(death_by5yr == 1))
mu_B <- with(subset(dat, race == "Black"), mean(death_by5yr == 1))

# Counterfactual Black risk under White coefficients
mu_B_given_Wbeta <- mean(predict(fit_ref,
                                 newdata = subset(dat, race == "Black"),
                                 type = "response"))

# Decomposition on the risk-difference scale
delta_total       <- mu_B - mu_W
delta_explained   <- mu_B_given_Wbeta - mu_W
delta_unexplained <- mu_B - mu_B_given_Wbeta
percent_explained <- 100 * delta_explained / delta_total
```

**Complex survey support.** Use the **Peters–Belson method for complex surveys** (weights/strata/clusters):

> Graubard, B. I., Rao, R. S., & Gastwirth, J. L. (2005). Using the Peters–Belson method to measure health care disparities from complex survey data. *Statistics in Medicine, 24*(17), 2659–2668. https://doi.org/10.1002/sim.2135

This paper details design-based estimation and inference for explained/unexplained components.




