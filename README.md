# MIC
Multilevel Integrative Clustering, packaged by Qian Li

# Installation:
```r
devtools::install_github("Qian-Li/MIC")
```

See our manuscript for more details!

# Quick Example

Time series simulation:
```r
ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)
```

MIC modle fitting:
```r
# Running time approx. 3 mins
output <- MIC(X = list_data, K = 4, NumRun = 5000)
```

Clustering result:
```r
# Accuracy relative to the true group labels
sum(clust_align(ts_sim$C, output$Cbest, type = 'vec') == ts_sim$C) / 40
```
