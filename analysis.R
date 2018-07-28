strikes<-read.csv("strikes_clean.csv")
# Country excluded in interaction terms
M0<-lm(log(Strike+1)~1,data=strikes)
Mfull<-lm(log(Strike+1)~(.-Centr-Country)^2+Country,data=strikes)
Mstart<-lm(log(Strike+1)~Year + Infl + Unemp+Dens+Demo+Country,data=strikes)
Mfwd<-step(object = M0, scope = list(lower = M0, upper = Mfull),
           direction = "forward", trace = FALSE)
Mbwd<-step(object = Mfull, scope = list(lower = M0, upper = Mfull),
           direction = "backward", trace = FALSE)
Mstep <- step(object = Mstart, scope = list(lower = M0, upper = Mfull),
              direction = "both", trace = FALSE)

#Country included in interaction terms
M01<-lm(log(Strike+1)~1,data=strikes)
Mfull1<-lm(log(Strike+1)~(.-Centr)^2,data=strikes)
Mstart1<-lm(log(Strike+1)~Year + Infl + Unemp,data=strikes)
Mfwd1<-step(object = M01, scope = list(lower = M01, upper = Mfull1),
            direction = "forward", trace = FALSE)

#Residuals
y.hat<-predict(M)
h<-hatvalues(M)
res<-resid(M)
Msum <- summary(M)
sigma.hat <- Msum$sigma
#studentized residuals
stan.res <- res/sigma.hat
stud.res<-studres(M)
#press residuals
press<-res/(1-h)
#DIFFTS
dfts<-dffits(M)

p <- length(coef(M))
n <- nobs(M)
hbar <- p/n # average leverage
stud.res <- stud.res*sqrt(1-hbar) # at h = hbar, stud.res = stan.res
press <- press*(1-hbar)/sigma.hat # at h = hbar, press = stan.res
dfts <- dfts*(1-hbar)/sqrt(hbar) # at h = hbar, dfts = stan.res
#plot all residuals against predicted values
par(mfrow = c(1,1), mar = c(4,4,1.5,1))
cex <- .8
plot(y.hat, rep(0, length(y.hat)), type = "n", # empty plot to get the axis range
     ylim = range(stan.res, stud.res, press, dfts), cex.axis = cex,
     xlab = "Predicted Values", ylab = "Residuals",main="M2")
# dotted line connecting each observations residuals for better visibility
segments(x0 = y.hat,
         y0 = pmin(stan.res, stud.res, press, dfts),
         y1 = pmax(stan.res, stud.res, press, dfts),
         lty = 2)
points(y.hat, stan.res, pch = 21, bg = "black", cex = cex)
points(y.hat, stud.res, pch = 21, bg = "blue", cex = cex)
points(y.hat, press, pch = 21, bg = "red", cex = cex)
points(y.hat, dfts, pch = 21, bg = "orange", cex = cex)
abline(h=0,lty=2)

# Plot residuals against leverages
plot(h, rep(0, length(y.hat)), type = "n", cex.axis = cex,
     ylim = range(stan.res, stud.res, press, dfts),
     xlab = "Leverages", ylab = "Residuals",main="M2")
segments(x0 = h,
         y0 = pmin(stan.res, stud.res, press, dfts),
         y1 = pmax(stan.res, stud.res, press, dfts),
         lty = 2)
points(h, stan.res, pch = 21, bg = "black", cex = cex)
points(h, stud.res, pch = 21, bg = "blue", cex = cex)
points(h, press, pch = 21, bg = "red", cex = cex)
points(h, dfts, pch = 21, bg = "orange", cex = cex)
abline(v = hbar, col = "grey60", lty = 2)
legend("topright", legend = c("Standardized", "Studentized", "PRESS", "DFFITS"),
       pch = 21, pt.bg = c("black", "blue", "red", "orange"), title = "Residual Type:",
       cex = 0.7, pt.cex = 1)

#residual histograms
hist(stud.res,prob=TRUE,main ="M2")
lines(seq(-4,4,by=0.1),dnorm(seq(-4,4,by=0.1)),col="red")
#Residual QQ-plots
qqnorm(studres(M),main="M2")
qqline(studres(M))
# cook's distance vs. leverage
D <- cooks.distance(M)
# flag some of the points
infl.ind <- which.max(D) # top influence point

lev.ind <- h > 2*hbar # leverage more than 2x the average
clrs <- rep("black", len = n)
clrs[lev.ind] <- "blue"
clrs[infl.ind] <- "red"
par(mfrow = c(1,1), mar = c(4,4,2,2))
cex <- 1
plot(h, D, xlab = "Leverage", ylab = "Cook's Influence Measure",
     pch = 21, bg = clrs, cex = cex, cex.axis = cex,main="M2")
p <- length(coef(M))
n <- 625
hbar <- p/n # average leverage
abline(v = 2*hbar, col = "grey60", lty = 2) # 2x average leverage
legend("topright", legend = c("High Leverage", "High Influence"), pch = 21,
       pt.bg = c("blue", "red"), cex = 0.75, pt.cex = cex)

#Cross-validation
M1 <- Mbwd
M2 <- Mbwd1
Mnames <- expression(M[1], M[2])

AIC1 <- AIC(M1)
AIC2 <- AIC(M2) # for M2
c(AIC1 = AIC1, AIC2 = AIC2) 

# PRESS
press1 <- resid(M1)/(1-hatvalues(M1)) # M1
press2 <- resid(M2)/(1-hatvalues(M2)) # M2
c(PRESS1 = sum(press1^2), PRESS2 = sum(press2^2)) 

# plot PRESS statistics
par(mar = c(3, 6, 1, 1))
boxplot(x = list(press1^2, press2^2), names = Mnames,
        ylab = expression(PRESS[i]^2), col = c("yellow", "orange"))

nreps <- 2e3 # number of replications
ntot <- nrow(strikes) # total number of observations
ntrain <- 500 # size of training set
ntest <- ntot-ntrain # size of test set
sse1 <- rep(NA, nreps) # sum-of-square errors for each CV replication
sse2 <- rep(NA, nreps)
Lambda <- rep(NA, nreps) # likelihod ratio statistic for each replication
system.time({
  for(ii in 1:nreps) {
    if(ii%%400 == 0) message("ii = ", ii)
    # randomly select training observations
    train.ind <- sample(ntot, ntrain) # training observations
    
    M1.cv <- update(M1, subset = train.ind)
    M2.cv <- update(M2, subset = train.ind)
    # testing residuals for both models
    # that is, testing data - predictions with training parameters
    M1.res <- log(Strike[-train.ind]+1) - predict(M1.cv, newdata = strikes[-train.ind,])
    M2.res <- log(Strike[-train.ind]+1) - predict(M2.cv, newdata = strikes[-train.ind,])
    # total sum of square errors
    sse1[ii] <- sum((M1.res)^2)
    sse2[ii] <- sum((M2.res)^2)
    # testing likelihood ratio
    M1.sigma <- sqrt(sum(resid(M1.cv)^2)/ntrain) # MLE of sigma
    M2.sigma <- sqrt(sum(resid(M2.cv)^2)/ntrain)
    Lambda[ii] <- sum(dnorm(M1.res, mean = 0, sd = M1.sigma, log = TRUE))
    Lambda[ii] <- Lambda[ii] - sum(dnorm(M2.res, mean = 0, sd = M2.sigma, log = TRUE))
  }
})

# plot cross-validation SSE and Lambda
par(mfrow = c(1,2))
par(mar = c(4.5, 4.5, .1, .1))
boxplot(x = list(sse1, sse2), names = Mnames, cex = .7,
        ylab = expression(SS[err]^{test}), col = c("yellow", "orange"))
hist(Lambda, breaks = 50, freq = FALSE, xlab = expression(Lambda^{test}),
     main = "", cex = .7)
abline(v = mean(Lambda), col = "red") # average value
