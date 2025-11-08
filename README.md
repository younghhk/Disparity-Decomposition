## Contents
- [1. Peters–Belson Decomposition of the Proportion of Death](#dd-proportion)
- [2. Quantile-Based Decomposition (Lower-Tail BMI / Telomere)](#dd-quantile)
- [3. Time-Dependent Disparity Decomposition](#dd-time)

<a id="dd-proportion"></a>
### 1. Peters–Belson Decomposition of the Proportion of Death

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



---

<a id="dd-quantile"></a>
### 2. Quantile-Based Disparity Decomposition (focus on **lower BMI**)

When the **outcome is continuous** and you care about the **lower tail** (e.g., “lower BMI” among White vs Black patients), use a **quantile-based Peters–Belson / Oaxaca–Blinder** decomposition. This uses **all data** (no subgroup fitting beyond race) and targets a chosen quantile level $\tau$ (e.g., $\tau=0.10$ or $0.25$).


Let $Y$ be BMI, $G\in\{\text{White},\text{Black}\}$, and $X$ covariates (e.g., stage, access, adherence).  
For a quantile level $\tau\in(0,1)$, define the **marginal** BMI quantiles:
- $Q_W(\tau)$ = $\tau$-quantile of BMI among Whites,
- $Q_B(\tau)$ = $\tau$-quantile of BMI among Blacks.

The **total disparity at the lower tail** is
$\Delta_{\text{total}}(\tau)\;=\;Q_B(\tau)-Q_W(\tau).$

### Steps

- **1)  Reference-group quantile model**
Fit a **quantile regression** in the **reference group** (commonly Whites) at the same $\tau$:
$Q_{Y\mid X,G=\text{White}}(\tau)\;=\;q_\tau\\big(X;\hat\beta_W(\tau)\big).$

- **2) Counterfactual lower-tail BMI for Blacks under White coefficients**

  Apply the White model to the **Black covariate distribution** to get the **counterfactual** $\tau$-quantile for Blacks (as if the conditional BMI–covariate relationship were the same as Whites):

$$
Q_{B\mid W}(\tau)=Q_{\tau}\{\,q_{\tau}(X_{Bi};\hat\beta_W(\tau))\,\}_i
$$

- **3) Decomposition at the chosen lower-tail quantile**

Split the lower-tail disparity into **explained** (covariates) and **unexplained** (differences in conditional relationships):


$$\Delta_{\mathrm{total}}(\tau)=\Delta_{\mathrm{expl}}(\tau)+\Delta_{\mathrm{unexpl}}(\tau)$$



with

$$
\Delta_{\mathrm{expl}}(\tau)=Q_{B\mid W}(\tau)-Q_W(\tau)
$$

$$
\Delta_{\mathrm{unexpl}}(\tau)=Q_B(\tau)-Q_{B\mid W}(\tau)
$$

- **Percent explained:**  
  $100\times\frac{\Delta_{\mathrm{expl}}(\tau)}{\Delta_{\mathrm{total}}(\tau)}$  (when $\Delta_{\mathrm{total}}(\tau)\neq 0$).

### R sketch (lower-tail BMI at $\tau=0.25$)

```r
# Outcome: BMI (continuous)
# Reference group: Whites
# Comparison group: Blacks
# Covariates: stage, access, adherence, x_d (replace as needed)

library(quantreg)

tau <- 0.25

# 1) Fit quantile regression in Whites at tau
fit_ref_tau <- rq(BMI ~ stage + access + adherence + x_d,
                  tau = tau,
                  data = subset(dat, race == "White"))

# 2) Observed marginal quantiles (lower tail) in each group
QW_tau <- quantile(subset(dat, race == "White")$BMI, probs = tau, na.rm = TRUE)
QB_tau <- quantile(subset(dat, race == "Black")$BMI,  probs = tau, na.rm = TRUE)

# 3) Counterfactual: predict conditional tau-quantiles for Blacks using White coefficients
pred_B_tau <- predict(fit_ref_tau, newdata = subset(dat, race == "Black"))

# Take the tau-quantile of those predictions (approximate marginal counterfactual)
QB_givenW_tau <- quantile(pred_B_tau, probs = tau, na.rm = TRUE)

# 4) Decomposition at the lower-tail quantile
delta_total_tau <- as.numeric(QB_tau - QW_tau)
delta_expl_tau  <- as.numeric(QB_givenW_tau - QW_tau)
delta_unex_tau  <- as.numeric(QB_tau - QB_givenW_tau)
pct_expl_tau    <- ifelse(delta_total_tau != 0, 100 * delta_expl_tau / delta_total_tau, NA)

list(QW_tau = QW_tau,
     QB_tau = QB_tau,
     QB_givenW_tau = QB_givenW_tau,
     delta_total_tau = delta_total_tau,
     delta_expl_tau  = delta_expl_tau,
     delta_unexpl_tau= delta_unex_tau,
     percent_explained = pct_expl_tau)
```
> Hong, G., Graubard, B., Gastwirth, J., & Kim, M. (2024) *Quantile Regression Decomposition Analysis of Disparity Research Using Complex Survey Data: Application to Disparities in BMI and Telomere Length Between U.S. Minority and White Population Groups.* [PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12456447/)


---
<a id="dd-time"></a>
### 3. Time-Dependent Disparity Decomposition 

**Goal:** Show how the **explained** and **unexplained** portions of a disparity change over time.  
Example: fetal growth (e.g., estimated fetal weight or abdominal circumference) for **White vs Black** mothers at **30** and **35** weeks’ gestation.

> Lee, S., Kim, S., Kim, M., & Hong, G. (2025) *Decomposition of Longitudinal Disparities: an Application to the Fetal Growth-Singletons Study*. [Biostatistics, in press]




---
