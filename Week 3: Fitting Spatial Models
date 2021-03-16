# In this session, we will:
# ??? Fit variograms (kriging);
# ??? Fit Gaussian processes with maximum likelihood and Bayesian methods;
# ??? Use the different models to predict/interpolate on a spatial grid.


# Kriging
library(geoR)
data(parana)
#summary(parana)
par(mar=c(4,2,2,2))
points(parana)

plot(parana)

sample_vario <- variog(parana, option = 'bin')
plot(sample_vario, pch=19)

sample_vario$n

sample_vario <- variog(parana, option='bin', max.dist=400)
plot(variog(parana, option='cloud', max.dist=400), pch = 19)

plot(variog(parana, option='bin', max.dist=400, bin.cloud = TRUE), bin.cloud = TRUE)

par(mar=c(4,4,2,2), mfrow=c(1,3))
plot(variog(parana, option='cloud', max.dist=400, estimator.type='modulus'), pch = 19)
sample.vario <- variog(parana, option='bin', max.dist=400, estimator.type='modulus', bin.cloud = TRUE)
plot(sample.vario, pch = 19)
plot(sample.vario, bin.cloud = TRUE)

vari.default <- variofit(sample.vario)
summary(vari.default)

vari.mat1.5 <- variofit(sample.vario, kappa=1.5, fix.kappa=TRUE)
summary(vari.mat1.5)
par(mar=c(4,4,2,2))
plot(sample.vario, pch = 19)
lines(vari.default)
lines(vari.mat1.5, lty = 2)

#variofit(sample.vario, kappa=1.5, fix.kappa=T, weights='cressie')
pred_grid <- expand.grid(seq(100, 800, by = 5), seq(0, 550, by = 5))
preds <- krige.conv(parana, loc=pred_grid, krige=krige.control(obj.model=vari.mat1.5))
image(preds, col = viridis::viridis(100), zlim = c(0,max(c(preds$predict))),
      coords.data = parana[1]$coords, main = 'Mean', xlab = 'x', ylab = 'y',
      x.leg = c(700, 900), y.leg = c(20, 70))
image(preds, values = preds$krige.var, col = heat.colors(100)[100:1],
      zlim = c(0,max(c(preds$krige.var))), coords.data = parana[1]$coords,
      main = 'Variance', xlab = 'x', ylab = 'y', x.leg = c(700, 900), y.leg = c(20, 70))


# Fitting Gaussian Processes
require(geoR)
set.seed(766)
ex.data <- grf(100, cov.pars=c(10, 0.2), cov.model="matern", kappa=1.5) # generating synthetic data

par(mar=c(4,2,2,2)); points(ex.data, pch = 17)

ml <- likfit(ex.data, ini=c(0.5, 0.5), fix.nug = TRUE, cov.model='matern', kappa=1.5)
summary(ml)

reml <- likfit(ex.data, ini=c(0.5, 0.5), fix.nug = TRUE, cov.model='matern', kappa=1.5,
               lik.met = "REML")
summary(reml)

par(mar=c(4,4,2,2))
plot(variog(ex.data), pch = 19)
lines.variomodel(ex.data, lwd=2, col = 'blue')
lines(ml)
lines(reml, lty = 2)

xv.ml <- xvalid(ex.data, model = ml)

par(mfcol = c(5,2), mar = c(4,4,1,1))
plot(xv.ml)

xvR.ml <- xvalid(ex.data, model = ml, reest = TRUE)

pred.grid <- expand.grid(seq(0, 1, l = 51), seq(0, 1, l = 51))
kc <- krige.conv(ex.data, loc = pred.grid, krige = krige.control(obj.model = ml))

par(mfcol=c(1,1),mar=c(4,4,2,2))
image(kc, col = viridis::viridis(100), zlim = c(-6,6), main = 'Mean - ML')
image(kc, values = kc$krige.var, zlim = c(0,1.65), coords.data = ex.data[1]$coords,
      main = 'Variance - ML')
kc_reml <- krige.conv(ex.data, loc = pred.grid, krige = krige.control(obj.model = reml))
image(kc_reml, col = viridis::viridis(100), zlim = c(-6,6), main = 'Mean - REML')
image(kc_reml, values = kc_reml$krige.var, zlim = c(0,1.65), coords.data = ex.data[1]$coords,
      main = 'Variance - REML')


# Bayesian Estimation
ex.grid <- as.matrix(expand.grid(seq(0,1,l=21), seq(0,1,l=21)))
ex.bayes <- krige.bayes(geodata = ex.data, loc=ex.grid,
                        model = model.control(cov.m="matern", kappa=1.5),
                        prior = prior.control(phi.discrete=seq(0, 1, l=51),
                                              phi.prior="reciprocal"))

par(mar=c(4,4,2,2))
plot(ex.bayes, type="h", tausq.rel = FALSE, col=c("red", "blue"))

par(mfrow=c(1,3), mar=c(4,4,2,2))
hist(ex.bayes, breaks = 15)

image(ex.bayes, col = viridis::viridis(100), zlim = c(-6,6), main="Mean", xlab = 'x', ylab = 'y')
image(ex.bayes, val="variance", main="Variance", xlab = 'x', ylab = 'y')
image(ex.bayes, val= "simulation", number.col=1, col = viridis::viridis(100), zlim = c(-6,6),
      main="Simulation 1", xlab = 'x', ylab = 'y')
image(ex.bayes, val= "simulation", number.col=2, col = viridis::viridis(100), zlim = c(-6,6),
      main="Simulation 2", xlab = 'x', ylab = 'y')

plot(variog(ex.data), pch = 19)
lines.variomodel(ex.data, lwd=2, col = 'blue')
lines(ml)
lines(reml, lty = 2)
lines(ex.bayes, summ = mean, lty=1, lwd=2, col = 'green')
lines(ex.bayes, summ = median, lty=2, lwd=2, col = 'green')
lines(ex.bayes, summary = "mode", post = "parameters", lty=3, lwd=2, col = 'green')








