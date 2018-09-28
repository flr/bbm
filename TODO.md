
# bbm() method design

- CREATE function to coerce FLStock + FLIndices to tmb 'model' object
  - Could this be a generic method for FLR?
- Returned model object can then be called with fit() or tmbstan()

```{r}
library(tmbstan)

library(bbm)

data(ane)

run <- BBM(catch=catch, indices=indices, control=control, inits=inits, debug=TRUE)

# browser() @ BBM:83

sfit <- tmbstan(model, chains=1, debug=TRUE, iter=50000, thin=100, warmup=1000)

pairs(sfit, pars=names(model$par)[1:5])
traceplot(sfit, pars=names(model$par), inc_warmup=FALSE)
```
