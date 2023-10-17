# install packages
library(gamair)
library(mgcv)
library(knitr)
library(kableExtra)
library(dplyr)

# load the CanWeather dataset
data(canWeather)

# regions
regions <- levels(CanWeather$region)

# EDA
## summary statistics
kbl(summary(CanWeather)[-7,-5], booktabs=T, linesep='') %>%
  kable_styling(latex_options = "hold_position")
## plotting
par(mfrow=c(2,2))
# Scatter plot: Temperature vs. Time
plot(CanWeather$time, CanWeather$T,
     main = "Temperature vs. Time", 
     col=adjustcolor('black',alpha=0.1),
     xlab = "Time", ylab = "Temperature (T)")
# Scatter plot: Temperature vs. Latitude
plot(CanWeather$latitude, CanWeather$T,
     main = "Temperature vs. Latitude",
     col=adjustcolor('black',alpha=0.1),
     xlab = "Latitude", ylab = "Temperature (T)")
# Box plot: Temperature vs. Region
boxplot(CanWeather$T ~ CanWeather$region,
        main = "Temperature vs. Region",
        xlab = "Region", ylab = "Temperature (T)")
# Box plot: Temperature vs. Place
boxplot(CanWeather$T ~ CanWeather$place, 
        main = "Temperature vs. Place",
        xlab = "Place", ylab = "Temperature (T)")

# fit spline regression model
gamFit <- function(location, i, p) {
  mod <- gam(T~s(time,bs='bs'), data=location)
  if (p) {
    mgcvfit <- predict(mod,data.frame(time=1:365),se.fit=TRUE)
    plot(location$T~location$time, xlab=paste(regions[i], 'Day in a Year'),
         ylab='Temperature in Celsius', col=adjustcolor('black', alpha=0.2),
         main = paste0("Spline Reg. fit on ", regions[i]))
    with(mgcvfit,lines(fit, col='red'))
    with(mgcvfit,lines(fit-2*se.fit,lty='dashed', col='orange'))
    with(mgcvfit,lines(fit+2*se.fit,lty='dashed', col='orange'))
  }
  return(mod)
}
# fit polynomial regression model
polyFit <- function(location, i, degree, p) {
  mod <- lm(location$T ~ poly(location$time,deg=degree))
  if (p) {
    plot(location$T~location$time, xlab=paste(regions[i], 'Day in a Year'),
         ylab='Temperature in Celsius', col=adjustcolor('black', alpha=0.2),
         main = paste0("Poly. Reg. fit on ", regions[i]))
    predint <- predict(mod, data.frame(rep(1:365, length(unique(location$place)))),
                       interval='confidence', level = 0.95)
    lines(predint[,1], col='red')
    lines(predint[,2],lty='dashed', col='orange')
    lines(predint[,3],lty='dashed', col='orange')
  }
  # compute gcv
  RSS <- sum(resid(mod)^2)
  n <- length(resid(mod))
  X <- model.matrix(mod)
  H <- X %*% solve(t(X) %*% X) %*% t(X)
  tr_H <- sum(diag(H))
  GCV <- RSS / (n * (1 - tr_H/n)^2)
  return(list(model = mod, gcv = GCV))
}

# fit global additive model
gam_glob <- gam(T ~ s(time, bs='bs') + s(latitude, bs='cr'), data=CanWeather)

# plot polynomial fits
for (i in 1:length(regions)) {
  polyFit(CanWeather[which(CanWeather$region == regions[i]), ], i, 6, T)
}
# plot spline fits
for (i in 1:length(regions)) {
  gamFit(CanWeather[which(CanWeather$region == regions[i]), ], i, T)
}
# plot global additive fits
for (i in 1:length(regions)) {
  region <- regions[i]
  latitudes_in_region <- CanWeather$latitude[CanWeather$region == region]
  location <- CanWeather[CanWeather$region == region,]
  mean_latitude <- mean(latitudes_in_region, na.rm = TRUE)
  
  mgcvfit <- predict(gam_glob,data.frame(time=1:365, latitude=mean_latitude), se.fit=TRUE)
  plot(location$T~location$time, xlab=paste(regions[i], 'Day in a Year'),
       ylab='Temperature in Celsius', col=adjustcolor('black', alpha=0.2),
       main = paste0("GAM: Temp vs. ", regions[i]))
  with(mgcvfit,lines(fit, col='red'))
  with(mgcvfit,lines(fit-2*se.fit,lty='dashed', col='orange', lwd=1.5))
  with(mgcvfit,lines(fit+2*se.fit,lty='dashed', col='orange', lwd=1.5))
}

# polynomial gcv and r2
gcv_scores_poly <- list()
r.sq_scores_poly <- list()
for (i in 1:length(regions)) {
  fit <- polyFit(CanWeather[which(CanWeather$region == regions[i]), ], i, 6, F)
  gcv_scores_poly[[regions[i]]] <- fit$gcv
  r.sq_scores_poly[[regions[i]]] <- summary(fit$model)$r.squared
}
# spline gcv and r2
gcv_scores_gam <- list()
r.sq_scores_gam <- list()
for (i in 1:length(regions)) {
  fit <- gamFit(CanWeather[which(CanWeather$region == regions[i]), ], i, F)
  gcv_scores_gam[[regions[i]]] <- summary(fit)$sp.criterion
  r.sq_scores_gam[[regions[i]]] <- summary(fit)$r.sq
}
# global additive model gcv and adj r2
r.sq_scores_gam_glob <- summary(gam_glob)$r.sq
gcv_gam_glob <- gam_glob$gcv.ubre

# organize in data frame
df <- data.frame(
  'Polynomial'=c(mean(unlist(r.sq_scores_poly)), mean(unlist(gcv_scores_poly))),
  'Spline'=c(mean(unlist(r.sq_scores_gam)), mean(unlist(gcv_scores_gam))),
  'GAM'=c(r.sq_scores_gam_glob, gcv_gam_glob))
rownames(df) <- c('R squared', 'GCV estimate')
kbl(df, booktabs=T, linesep='') %>%
  kable_styling("hold_position")
