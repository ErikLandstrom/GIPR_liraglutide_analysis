---
title: "GIPR statistical analysis 5vv"
author: "Erik Ländström"
date: "8 March 2019"
output: 
  html_document: 
    keep_md: yes
---



Read libraries.


```r
library(tidyverse)
library(modelr)
library(broom)
library(ggfortify)
library(sva)
library(qvalue)
```

Read functions.



# Analysis

Perform statistical analysis of imputed data with at least 5 valid values per 
group.

Read data.


```r
data <- read_tsv("GIPR_liraglutide_liver_5vv_imputed.tsv")
```

Separate conditions into separate columns.


```r
factored_data <- separate_condition_from_DEP_output(data, c("treatment", "gender", "series", "litter"), remove = FALSE)
```


## Linear model

Nest data.


```r
nested_data <- factored_data %>% 
  group_by(name) %>% 
  nest()
```

Linear model function.


```r
protein_model <- function(tb) {
  lm(LFQ ~ gender * litter, data = tb)
}
```

Perform linear regression on each protein and calculate/extract residuals.


```r
nested_model_data <- nested_data %>% 
  mutate(model = map(data, protein_model))

nested_resids_data <- nested_model_data %>% 
  mutate(resids = map2(data, model, add_residuals))

resids <- unnest(nested_resids_data, resids)
```

Calculate variance and compare with original data.


```r
variance <- resids %>% 
  group_by(name) %>% 
  summarise(var_lfq = var(LFQ),
            var_resid = var(resid))

variance %>% 
  summarise(mean_lfq = mean(var_lfq),
            mean_resid = mean(var_resid))
```

```
## # A tibble: 1 x 2
##   mean_lfq mean_resid
##      <dbl>      <dbl>
## 1    0.273      0.144
```

```r
# Test if the regression reduced the variance for all proteins
variance %>% 
  mutate(smaller = if_else(var_resid < var_lfq, 1, 0)) %>% 
  summarise(n = sum(smaller),
            all = n())
```

```
## # A tibble: 1 x 2
##       n   all
##   <dbl> <int>
## 1  2611  2611
```

```r
resids %>% 
  group_by(name, treatment) %>% 
  summarise(var_lfq = var(LFQ),
            var_resid = var(resid))
```

```
## # A tibble: 5,222 x 4
## # Groups:   name [?]
##    name  treatment var_lfq var_resid
##    <chr> <chr>       <dbl>     <dbl>
##  1 A1BG  L         0.477     0.198  
##  2 A1BG  P         0.255     0.0872 
##  3 A1CF  L         0.00643   0.00311
##  4 A1CF  P         0.0186    0.00439
##  5 A2M   L         0.453     0.304  
##  6 A2M   P         0.321     0.134  
##  7 AACS  L         0.120     0.122  
##  8 AACS  P         0.218     0.0594 
##  9 AADAC L         0.0302    0.0624 
## 10 AADAC P         0.0592    0.0445 
## # ... with 5,212 more rows
```

The linear model decreased the variance for every protein!

## PCA

Unite replicate (ID) with condition (grouping), and spread to intensities into
separate columns for each sample.


```r
lfq_for_pca <- resids %>% 
  dplyr::select(name, replicate, condition, LFQ) %>% 
  unite(replicate, replicate, condition) %>% 
  spread(key = "replicate", value = "LFQ")

resids_for_pca <- resids %>% 
  dplyr::select(name, replicate, condition, resid) %>% 
  unite(replicate, replicate, condition) %>% 
  spread(key = "replicate", value = "resid")
```


Perform PCA for both LFQ and residuals.


```r
lfq_pca <- format_data_for_PCA(lfq_for_pca, 17, c("sample", "treatment", "gender", "series", "litter"))
```

```
## Warning: `as_tibble.matrix()` requires a matrix with column names or a `.name_repair` argument. Using compatibility `.name_repair`.
## This warning is displayed once per session.
```

```r
resid_pca <- format_data_for_PCA(resids_for_pca, 17, c("sample", "treatment", "gender", "series", "litter"))
```


Calculate the variance explained by each principal component for original and
data processed with linear regression.


```r
# Calculate the variance that each PC explains
lfq_var_exp <- calculate_variance_explained(lfq_pca)
make_var_exp_kable(lfq_var_exp)
```

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Proportion of the variance that each principal component explains</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> PC </th>
   <th style="text-align:center;"> Variance </th>
   <th style="text-align:center;"> Proportion explained </th>
   <th style="text-align:center;"> Cumulative proportion explained </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 450.4394 </td>
   <td style="text-align:center;"> 0.1725 </td>
   <td style="text-align:center;"> 0.1725 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 302.6057 </td>
   <td style="text-align:center;"> 0.1159 </td>
   <td style="text-align:center;"> 0.2884 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 210.7004 </td>
   <td style="text-align:center;"> 0.0807 </td>
   <td style="text-align:center;"> 0.3691 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 193.2114 </td>
   <td style="text-align:center;"> 0.0740 </td>
   <td style="text-align:center;"> 0.4431 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 167.4022 </td>
   <td style="text-align:center;"> 0.0641 </td>
   <td style="text-align:center;"> 0.5072 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 155.7387 </td>
   <td style="text-align:center;"> 0.0596 </td>
   <td style="text-align:center;"> 0.5669 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 145.5238 </td>
   <td style="text-align:center;"> 0.0557 </td>
   <td style="text-align:center;"> 0.6226 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 133.9327 </td>
   <td style="text-align:center;"> 0.0513 </td>
   <td style="text-align:center;"> 0.6739 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 129.1350 </td>
   <td style="text-align:center;"> 0.0495 </td>
   <td style="text-align:center;"> 0.7234 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 119.7712 </td>
   <td style="text-align:center;"> 0.0459 </td>
   <td style="text-align:center;"> 0.7692 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 115.7007 </td>
   <td style="text-align:center;"> 0.0443 </td>
   <td style="text-align:center;"> 0.8135 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 108.5180 </td>
   <td style="text-align:center;"> 0.0416 </td>
   <td style="text-align:center;"> 0.8551 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 100.9678 </td>
   <td style="text-align:center;"> 0.0387 </td>
   <td style="text-align:center;"> 0.8938 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 99.0919 </td>
   <td style="text-align:center;"> 0.0380 </td>
   <td style="text-align:center;"> 0.9317 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 90.0665 </td>
   <td style="text-align:center;"> 0.0345 </td>
   <td style="text-align:center;"> 0.9662 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 88.1947 </td>
   <td style="text-align:center;"> 0.0338 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
</tbody>
</table>

```r
resid_var_exp <- calculate_variance_explained(resid_pca)
make_var_exp_kable(resid_var_exp)
```

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Proportion of the variance that each principal component explains</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> PC </th>
   <th style="text-align:center;"> Variance </th>
   <th style="text-align:center;"> Proportion explained </th>
   <th style="text-align:center;"> Cumulative proportion explained </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 514.2609 </td>
   <td style="text-align:center;"> 0.1970 </td>
   <td style="text-align:center;"> 0.1970 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 437.5550 </td>
   <td style="text-align:center;"> 0.1676 </td>
   <td style="text-align:center;"> 0.3645 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 292.9518 </td>
   <td style="text-align:center;"> 0.1122 </td>
   <td style="text-align:center;"> 0.4767 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 286.6988 </td>
   <td style="text-align:center;"> 0.1098 </td>
   <td style="text-align:center;"> 0.5865 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 244.5556 </td>
   <td style="text-align:center;"> 0.0937 </td>
   <td style="text-align:center;"> 0.6802 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 241.0802 </td>
   <td style="text-align:center;"> 0.0923 </td>
   <td style="text-align:center;"> 0.7725 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 209.5177 </td>
   <td style="text-align:center;"> 0.0802 </td>
   <td style="text-align:center;"> 0.8528 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 196.7997 </td>
   <td style="text-align:center;"> 0.0754 </td>
   <td style="text-align:center;"> 0.9282 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 187.5803 </td>
   <td style="text-align:center;"> 0.0718 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
</tbody>
</table>

```r
# Plot variance explained plots
plot_var_exp(lfq_var_exp)
```

![](GIPR_analysis_with_sva_files/figure-html/var_exp-1.png)<!-- -->

```r
plot_var_exp(resid_var_exp)
```

![](GIPR_analysis_with_sva_files/figure-html/var_exp-2.png)<!-- -->

The principal components from the linear model explains considerably more of the
variance than the LFQ intensities.


Plot the first two PCs of the PCA.


```r
plot_pca(lfq_pca, "LFQ", "treatment")
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pca-1.png)<!-- -->

```r
plot_pca(resid_pca, "resids", "treatment")
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pca-2.png)<!-- -->



Plot PC1 to PC4.


```r
plot_multiple_PCAs(lfq_pca, 1, 4, "treatment")
```

```
## [[1]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-1.png)<!-- -->

```
## 
## [[2]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-2.png)<!-- -->

```
## 
## [[3]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-3.png)<!-- -->

```r
plot_multiple_PCAs(resid_pca, 1, 4, "treatment")
```

```
## [[1]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-4.png)<!-- -->

```
## 
## [[2]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-5.png)<!-- -->

```
## 
## [[3]]
```

![](GIPR_analysis_with_sva_files/figure-html/plot_pcas-6.png)<!-- -->


# Surrogate variable analysis

Format data for sva.

Create expression data.


```r
exp_data <- resids %>% 
  dplyr::select(name, replicate, LFQ) %>% 
  spread(key = "replicate", value = "LFQ")

row.names(exp_data) <- exp_data$name
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
exp_data <- exp_data %>% 
  dplyr::select(-name)
```

Create phenotypic data.


```r
pheno_data <- resids %>% 
  dplyr::select(treatment, gender, series, litter, replicate) %>% 
  distinct() %>% 
  mutate(sample = row_number())

rownames(pheno_data) <- pheno_data$replicate
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
pheno_data <- pheno_data %>% 
  dplyr::select(sample, treatment, gender, series, litter)
```

Create linear models, full and null model.


```r
mod = model.matrix(~ as.factor(treatment), data = pheno_data)

mod0 = model.matrix(~ 1,data = pheno_data)
```

Sva using LFQ values


```r
svobj = sva(as.matrix(exp_data), mod, mod0, n.sv= 2)
```

```
## Number of significant surrogate variables is:  2 
## Iteration (out of 5 ):1  2  3  4  5
```

```r
svobj[[1]]
```

```
##              [,1]        [,2]
##  [1,]  0.28121621 -0.32240858
##  [2,]  0.12085909  0.08103956
##  [3,]  0.32955982 -0.33477472
##  [4,]  0.31803672  0.13682510
##  [5,]  0.20249532  0.19872057
##  [6,]  0.31024047  0.38598444
##  [7,]  0.17506837 -0.03399454
##  [8,] -0.14610477  0.03295110
##  [9,] -0.08279752  0.17625184
## [10,] -0.18480657 -0.10680021
## [11,]  0.07074298 -0.07001630
## [12,]  0.05278217  0.20590248
## [13,] -0.22228838 -0.26804769
## [14,] -0.33120371  0.05219279
## [15,] -0.33161042 -0.23091697
## [16,] -0.40640025  0.46423540
## [17,] -0.15578952 -0.36714426
```

Sva using LFQ plus the known batch effects: gender and series.


```r
mod = model.matrix(~ as.factor(treatment) * as.factor(gender) * as.factor(series), data = pheno_data)

mod0 = model.matrix(~ as.factor(gender) * as.factor(series),data = pheno_data)

svobj = sva(as.matrix(exp_data), mod, mod0, n.sv= 2)
```

```
## Number of significant surrogate variables is:  2 
## Iteration (out of 5 ):1  2  3  4  5
```

```r
svobj[[1]]
```

```
##             [,1]        [,2]
##  [1,]  0.3205847 -0.25132052
##  [2,]  0.1688711  0.11124770
##  [3,] -0.1278212 -0.59875472
##  [4,] -0.2325753 -0.27956608
##  [5,] -0.2342868 -0.13889789
##  [6,]  0.1903921  0.18744818
##  [7,] -0.2605272 -0.29652097
##  [8,]  0.1921298  0.20590660
##  [9,] -0.2666725  0.19462981
## [10,]  0.1491838  0.09577549
## [11,]  0.1481339 -0.03091534
## [12,]  0.1605083  0.12540432
## [13,]  0.2184897  0.08806001
## [14,] -0.3128468  0.12628874
## [15,]  0.1190174  0.12595411
## [16,] -0.4867301  0.43553185
## [17,]  0.2541492 -0.10027129
```

Sva analysis using residuals from linear model.


```r
exp_resids <- resids %>% 
  dplyr::select(name, replicate, resid) %>% 
  spread(key = "replicate", value = "resid")

row.names(exp_resids) <- exp_resids$name
```

```
## Warning: Setting row names on a tibble is deprecated.
```

```r
exp_resids <- exp_resids %>% 
  dplyr::select(-name)

mod = model.matrix(~ as.factor(treatment), data = pheno_data)

mod0 = model.matrix(~ 1,data = pheno_data)

svobj = sva(as.matrix(exp_resids), mod, mod0, n.sv= 2)
```

```
## Number of significant surrogate variables is:  2 
## Iteration (out of 5 ):1  2  3  4  5
```

```r
svobj[[1]]
```

```
##                [,1]          [,2]
##  [1,]  4.342562e-16 -1.345525e-15
##  [2,] -1.450184e-01 -2.865308e-02
##  [3,]  3.107872e-01 -4.632327e-01
##  [4,]  5.423766e-02 -1.905271e-01
##  [5,]  5.666220e-03 -3.006328e-02
##  [6,] -2.256727e-01  7.124762e-01
##  [7,]  2.624369e-01 -7.559128e-02
##  [8,] -2.624369e-01  7.559128e-02
##  [9,] -2.542256e-16  3.062080e-16
## [10,] -1.033075e-01  8.761065e-02
## [11,]  1.033075e-01 -8.761065e-02
## [12,] -8.395223e-17  1.885715e-16
## [13,]  3.888464e-01  2.244531e-01
## [14,] -3.937661e-01 -2.817511e-01
## [15,]  4.919751e-03  5.729805e-02
## [16,] -4.260929e-01 -1.999281e-01
## [17,]  4.260929e-01  1.999281e-01
```

```r
sv_factor1 <- tibble(replicate = as.numeric(rownames(pheno_data)),
       sv1 = svobj[[1]][, 1])
```

Join sva factor sv1 with previous results.


```r
resids_sv <- full_join(resids, sv_factor1, by = "replicate")
```

Nest data.


```r
nested_data <- resids_sv %>% 
  group_by(name) %>% 
  rename(resid_old = resid) %>% 
  nest()
```

Linear model including the sv1 factor.


```r
protein_model <- function(tb) {
  lm(resid_old ~ sv1, data = tb)
}
```

Perform linear model.


```r
nested_model_data <- nested_data %>% 
  mutate(model = map(data, protein_model))

nested_resids_data <- nested_model_data %>% 
  mutate(resids = map2(data, model, add_residuals))

residuals_sv <- unnest(nested_resids_data, resids)
```

Calculate variance.


```r
variance <- resids %>% 
  group_by(name) %>% 
  summarise(var_lfq = var(LFQ),
            var_resid = var(resid))

variance %>% 
  summarise(mean_var_lfq = mean(var_lfq),
            mean_var_resid = mean(var_resid))
```

```
## # A tibble: 1 x 2
##   mean_var_lfq mean_var_resid
##          <dbl>          <dbl>
## 1        0.273          0.144
```

```r
variance %>% 
  mutate(smaller = if_else(var_resid < var_lfq, 1, 0)) %>% 
  summarise(n = sum(smaller))
```

```
## # A tibble: 1 x 1
##       n
##   <dbl>
## 1  2611
```

```r
resids %>% 
  group_by(name, treatment) %>% 
  summarise(var_lfq = var(LFQ),
            var_resid = var(resid))
```

```
## # A tibble: 5,222 x 4
## # Groups:   name [?]
##    name  treatment var_lfq var_resid
##    <chr> <chr>       <dbl>     <dbl>
##  1 A1BG  L         0.477     0.198  
##  2 A1BG  P         0.255     0.0872 
##  3 A1CF  L         0.00643   0.00311
##  4 A1CF  P         0.0186    0.00439
##  5 A2M   L         0.453     0.304  
##  6 A2M   P         0.321     0.134  
##  7 AACS  L         0.120     0.122  
##  8 AACS  P         0.218     0.0594 
##  9 AADAC L         0.0302    0.0624 
## 10 AADAC P         0.0592    0.0445 
## # ... with 5,212 more rows
```

Calculate variance explained.


```r
resids_sv_for_pca <- residuals_sv %>% 
  dplyr::select(name, replicate, condition, resid) %>% 
  unite(replicate, replicate, condition) %>% 
  spread(key = "replicate", value = "resid")

resid_sv_pca <- format_data_for_PCA(resids_sv_for_pca, 17, c("sample", "treatment", "gender", "series", "litter"))

resid_sv_var_exp <- calculate_variance_explained(resid_sv_pca)
make_var_exp_kable(resid_sv_var_exp)
```

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Proportion of the variance that each principal component explains</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> PC </th>
   <th style="text-align:center;"> Variance </th>
   <th style="text-align:center;"> Proportion explained </th>
   <th style="text-align:center;"> Cumulative proportion explained </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 579.1862 </td>
   <td style="text-align:center;"> 0.2218 </td>
   <td style="text-align:center;"> 0.2218 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 387.5179 </td>
   <td style="text-align:center;"> 0.1484 </td>
   <td style="text-align:center;"> 0.3702 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 350.5730 </td>
   <td style="text-align:center;"> 0.1343 </td>
   <td style="text-align:center;"> 0.5045 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 297.1915 </td>
   <td style="text-align:center;"> 0.1138 </td>
   <td style="text-align:center;"> 0.6183 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 283.6084 </td>
   <td style="text-align:center;"> 0.1086 </td>
   <td style="text-align:center;"> 0.7270 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 258.3997 </td>
   <td style="text-align:center;"> 0.0990 </td>
   <td style="text-align:center;"> 0.8259 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 236.5208 </td>
   <td style="text-align:center;"> 0.0906 </td>
   <td style="text-align:center;"> 0.9165 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 218.0025 </td>
   <td style="text-align:center;"> 0.0835 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
</tbody>
</table>

```r
plot_var_exp(resid_sv_var_exp)
```

![](GIPR_analysis_with_sva_files/figure-html/var_exp2-1.png)<!-- -->


```r
plot_pca(resid_sv_pca, "resids", "treatment")
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
plot_multiple_PCAs(resid_pca, 1, 4, "treatment")
```

```
## [[1]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

```
## 
## [[2]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-3.png)<!-- -->

```
## 
## [[3]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-4.png)<!-- -->

```r
plot_multiple_PCAs(resid_sv_pca, 1, 4, "treatment")
```

```
## [[1]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-5.png)<!-- -->

```
## 
## [[2]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-6.png)<!-- -->

```
## 
## [[3]]
```

![](GIPR_analysis_with_sva_files/figure-html/unnamed-chunk-3-7.png)<!-- -->

# t-test

Perform multiple t-test.


```r
residuals_sv_ttest <- multiple_ttests(residuals_sv, treatment,
                group_1 = L,
                group_2 = P,
                mean_1 = mean_L,
                mean_2 = mean_P,
                values = resid,
                equal_var = TRUE)
```

Plot p-value histogram.


```r
plot_pvalue_histogram(residuals_sv_ttest, p_value, limit_max = 600)
```

![](GIPR_analysis_with_sva_files/figure-html/phist-1.png)<!-- -->


```r
fdr_correction_with_qvalue <- function(tb, p_value_col) {
  
  # Quote
  p_value_col <- enquo(p_value_col)
  
  # fdr correction
  qvalues <- qvalue(as_vector(tb %>%
                                ungroup() %>% # To make sure to only get one column
                                dplyr::select(!!p_value_col)))
  
  return(qvalues)
}
```

FDR correction

```r
q_vals <- fdr_correction_with_qvalue(residuals_sv_ttest, p_value) 
```


```r
resids_ttest <- multiple_ttests(resids_sv, treatment,
                group_1 = L,
                group_2 = P,
                mean_1 = mean_L,
                mean_2 = mean_P,
                values = resid,
                equal_var = TRUE)
```


```
## R version 3.5.2 (2018-12-20)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United Kingdom.1252 
## [2] LC_CTYPE=English_United Kingdom.1252   
## [3] LC_MONETARY=English_United Kingdom.1252
## [4] LC_NUMERIC=C                           
## [5] LC_TIME=English_United Kingdom.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] kableExtra_1.0.1    knitr_1.21          bindrcpp_0.2.2     
##  [4] qvalue_2.12.0       sva_3.28.0          BiocParallel_1.14.2
##  [7] genefilter_1.62.0   mgcv_1.8-26         nlme_3.1-137       
## [10] ggfortify_0.4.5     broom_0.5.1         modelr_0.1.3       
## [13] forcats_0.3.0       stringr_1.4.0       dplyr_0.7.8        
## [16] purrr_0.3.0         readr_1.3.1         tidyr_0.8.2        
## [19] tibble_2.0.1        ggplot2_3.1.0       tidyverse_1.2.1    
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-6         matrixStats_0.54.0   lubridate_1.7.4     
##  [4] bit64_0.9-7          webshot_0.5.1        httr_1.4.0          
##  [7] tools_3.5.2          backports_1.1.3      utf8_1.1.4          
## [10] R6_2.3.0             DBI_1.0.0            lazyeval_0.2.1      
## [13] BiocGenerics_0.26.0  colorspace_1.4-0     withr_2.1.2         
## [16] tidyselect_0.2.5     gridExtra_2.3        bit_1.1-14          
## [19] compiler_3.5.2       cli_1.0.1            rvest_0.3.2         
## [22] Biobase_2.40.0       xml2_1.2.0           labeling_0.3        
## [25] scales_1.0.0         digest_0.6.18        rmarkdown_1.11.3    
## [28] pkgconfig_2.0.2      htmltools_0.3.6      highr_0.7           
## [31] limma_3.36.5         rlang_0.3.1          readxl_1.2.0        
## [34] rstudioapi_0.9.0     RSQLite_2.1.1        bindr_0.1.1         
## [37] generics_0.0.2       jsonlite_1.6         RCurl_1.95-4.11     
## [40] magrittr_1.5         Matrix_1.2-15        Rcpp_1.0.0          
## [43] munsell_0.5.0        S4Vectors_0.18.3     fansi_0.4.0         
## [46] stringi_1.2.4        yaml_2.2.0           plyr_1.8.4          
## [49] grid_3.5.2           blob_1.1.1           ggrepel_0.8.0       
## [52] parallel_3.5.2       crayon_1.3.4         lattice_0.20-38     
## [55] haven_2.0.0          splines_3.5.2        annotate_1.58.0     
## [58] hms_0.4.2            pillar_1.3.1         reshape2_1.4.3      
## [61] stats4_3.5.2         XML_3.98-1.17        glue_1.3.0          
## [64] evaluate_0.12        cellranger_1.1.0     gtable_0.2.0        
## [67] assertthat_0.2.0     xfun_0.4             xtable_1.8-3        
## [70] survival_2.43-3      viridisLite_0.3.0    AnnotationDbi_1.42.1
## [73] memoise_1.1.0        IRanges_2.14.12
```


