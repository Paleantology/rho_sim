---
title: "Simulate Data Sets"
author: "April Wright"
date: "1/17/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Simulate the FBD trees


```{r}
lambda = 1
mu = 0.2
tips = 100
trees = c()
num = 1

while (num < 10){
t = FossilSim::sim.fbd.age(numbsim = 1, lambda5 = .1, mu = .05, psi = .05, age = 100)
t
trees = c(trees, t)
num = num + 1
}
```

## Simulate fossils
```{r}
beta = 0.5
lambda.a = 1.2
num=1
for (t in trees) {
s = sim.taxonomy(tree = t, beta = beta, lambda.a = lambda.a)
f = sim.fossils.poisson(rate = rate, taxonomy = s)
file_string <- num + "/" + "num/ages.tsv"
write.csv(filestring, sep = "\t")
}
```

