rm(list=ls())                

library(quantmod)
library(quadprog)
library(zoo)
library(tseries)
library(ggplot2)
library(PerformanceAnalytics)
library(psych)

###############################################################
# This part we will download ETF data set from Yahoo Finance.
# We use library(quantmod) and library(tseries) to download data
# for selected ETF into data environment. In order to ease the 
# impact of stock split and dividends, we use adjusted close 
# price of ETFs. Moreover, to test the change of the Market
# environment we use VIX, which is the volatility of
# SPY, to be the indicator of the change of market environment.
###############################################################
# select ETFs
sym_bols <- c("FXE", "EWJ","GLD","QQQ","SPY",
              "SHV","DBA","USO", "XBI","ILF","GAF","EPP","FEZ")

data_env <- new.env()  # new environment for data
# download data for sym_bols into data_env
getSymbols(sym_bols, env=data_env, adjust=TRUE,
           from="2006-07-12", to="2015-03-01")
# combine daily adjusted close prices
# If close prices are needed, changing Ad into Cl at next function.
etf_series_ad <- do.call(merge,
                         eapply(data_env, Ad)[sym_bols])
colnames(etf_series_ad) <-
  sapply(colnames(etf_series_ad),
         function(col_name)
           strsplit(col_name, split="[.]")[[1]])[1, ]
# get ^VIX
VIX <- get.hist.quote(instrument = "^VIX", start = "2006-07-12", 
                      end = "2015-03-01", quote = "Close")
VIX <- matrix(VIX, nrow=2173,ncol=1)
# save(data_env, etf_series_ad, sym_bols, file='etf_data.Rdata')
# save(VIX, file='vix.Rdata')

###############################################################
# Next, we select 60 days for short term time period and 120
# days for long term time period to be the sample estimators.
# Then we use selected ETF adjusted close prices to calculate
# daily log returns and annually expected log returns for 60 and
# 120 days and covariance matrixes for 60 and 120 days. 
###############################################################
# load("etf_data.Rdata")
# load("vix.Rdata")

# change the time series into matrix for further calculation
data<-matrix(etf_series_ad,nrow=2173,ncol=13)
# Scrub NA
# Since some ETFs have no data in early time, we set NA into 0.
# The missing data, which are very little compared with large historical
# data, would not impact our strategies.
for(i in 1:13){
  for(j in 1:2173){
    if(is.na(data[j,i])) data[j,i]=0
  }
}

# Calculate daily log returns and 60 and 120 days annually expected log returns
log_r<-matrix(0,nrow=2173,ncol=13)
log_r_60<-matrix(0,nrow=2173,ncol=13)
log_r_120<-matrix(0,nrow=2173,ncol=13)
for (i in 1:13){
  for (j in 2:2173){
    if(data[j-1,i]==0) log_r[j,i]=0 else log_r[j,i]<-(data[j,i]/data[(j-1),i]-1)
    if (j>=61) log_r_60[j,i]<-mean(log_r[(j-59):j,i])*250
    if (j>=121) log_r_120[j,i]<-mean(log_r[(j-119):j,i])*250 
  }
}

# Calculate covariance matrixes for 60 and 120 days
cov_60<-array(0,dim=c(2173,13,13))
cov_120<-array(0,dim=c(2173,13,13))
# j is the number of the date
for (j in 61:2173){
  cov_60[j,,]<-cov(log_r[(j-59):j,1:13])*250
  if(j>120) cov_120[j,,]<-cov(log_r[(j-119):j,1:13])*250
}

#########################################################################
# In this part we build our first optimized portfolio for trend
# estimators. we choose R60 and C120, which are expected returns
# for 60 days and covariance for 120 days.
# Later on we will change the trend estimators from R60 and C120
# into the trend estimators for R120 and C120 and the trend estimators
# for R120 and C60 and the trend estimators for R60 and C60.
# Therefore, we have four combinations of estimators to optimize portfolios,
# which are (R60, C120), (R120, C120), (R120, C60) and (R60, C60).
# The optimized portfolio is based on Markowitz Optimization, minimizing
# the portfolio variance.

# Market indicator:
# As we previous explanation, we choose VIX index to be an indicator
# of the change of market environment.
# When VIX >= 25, we consider the Market is turbulence. Although we do not
# test whether VIX get into bubble or not, we can see the Market gets into
# an unstable situation. Then we adjust the frequency of rebalance weekly
# for 5 days.
# Otherwise, when VIX < 25, the Market is relative stable. We adjust
# the frequency of rebalance monthly for 20 days.
# Under the constraint, we build long term and short term portfolios
# using a combination of estimators to see the behavior of the strategies
# during different historical periods, which are before the 2008 crisis,
# during the 2008 crisis and after the crisis.

# Assumption:
# We do neglect transaction cost and do not consider SHV as a Riskless Security.
#########################################################################
# Sol1 is the weight set of R60 and C120
Sol1<-array(0,dim=c(2193,13))
A<-diag(1,nrow=13,ncol=13)
C<-matrix(c(A,-A),nrow=13,ncol=26)
B<-array(1,dim=c(13,1))
dvec<-array(0,dim=c(13,1))
# set 15% annual return Target and restricted weights to the interval
# between [-200%, 200%]
bvec<-c(1,0.15,-2*B,-2*B)
# use to test whether the sum of weights is 1
temp<-0
# the realized performance is from 178th row to 2173th row
i=178
while (i <= 2173){
  if(VIX[i]>=25){
    Amat<-matrix(c(B,log_r_60[i,],C),nrow=13,ncol=28)
    Dmat<-cov_120[i,,]
    for(j in 0:4)
      Sol1[i+j,]<-solve.QP(Dmat,dvec,Amat,bvec,meq=2)$solution
    #check the solution
    if(sum(Sol1[i,])-1>0.0001) temp<-temp+1
    i=i+5
    }
  else{
  Amat<-matrix(c(B,log_r_60[i,],C),nrow=13,ncol=28)
  Dmat<-cov_120[i,,]
  for(j in 0:19)
    Sol1[i+j,]<-solve.QP(Dmat,dvec,Amat,bvec,meq=2)$solution
  #check the solution
  if(sum(Sol1[i,])-1>0.0001) temp<-temp+1
  i=i+20
  }
}
##
# change (R60, C120) into (R120, C120) and (R120, C60) and (R60, C60) 
# to get weight sets, which are Sol2 with (R120, C120), Sol3 with
# (R120,C60) and Sol4 with (R60, C60).
# save(Sol1, Sol2, Sol3, Sol4, log_r, file='etf_weights.Rdata')

#########################################################################
# Finally we calculate Key Indicators for each ETF and our Optimal Strategies
# over the whole period from 03/23/2007 to 02/27/2015.
# The following indicators are Cumulated PnL, Daily Mean Geometric Return,
# Daily Min Return, Daily Max Return, Max Drawdown, Volatility, Sharpe Ratio,
# Sortino Ratio, Skewness, Kurtosis, Modified VaR and CVaR.
# We calculate the indicators based on the first strategy with the trend
# estimators (R60, C120). Then we calculate the indicators based on other
# three strategies with the trend estimators (R120, C120), (R120, C60) and
# (R60, C60) separately.
#########################################################################
# load("etf_data.Rdata")
# load("vix.Rdata")

# Weights: Sol1(R60, C120), Sol2(R120, C120), Sol3(R120, C60), Sol4(R60, C60)
load("etf_weights.Rdata")
Sol1 <- Sol1[1:2173,1:13]
Sol2 <- Sol2[1:2173,1:13]
Sol3 <- Sol3[1:2173,1:13]
Sol4 <- Sol4[1:2173,1:13]

# Plot one ETF to view the behavior of the underlying strategy.
# Later on we will plot all assets of the strategy with trend 
# estimators (R60, C120) as an example to see the trading behaviors.
FXE <- xts(Sol1[,1], order.by = index(etf_series_ad))
plot(FXE, type = "l")

# Performance and Risk Reporting
# Cumulated PnL or Return
PL <- vector("numeric",2173)
PL[178] <- 1
for (i in 179:2173){
  PL[i] <- PL[i-1]*(1+sum(Sol1[i,]*log_r[i,]))
}
PL <- xts(PL, order.by=index(etf_series_ad))
plot(PL[179:2173],type='l',main="Cumulated daily PnL for R60 and C120")

# Calculate optimal portfolio daily log returns
Return_O <- vector("numeric",2173)
for (i in 179:2173){
  Return_O[i] <- sum(Sol1[i,]*log_r[i,])
}

# Daily Mean Geometric Return
Geometric <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Geometric[1,i] <- geometric.mean(1+log_r[178:2173,i])-1
  #  Geometric[1,i] <- (Geometric[1,i]+1)^250-1 # annualized return
}
Geometric[1,14] <- geometric.mean(1+Return_O[178:2173])-1
#Geometric[1,14] <- (Geometric[1,14]+1)^250-1 # annualized return

# Daily Min Return
MinR <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  MinR[1,i] <- min(log_r[178:2173,i]) 
}
MinR[1,14] <- min(Return_O[178:2173])

# Daily Max Return
MaxR <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  MaxR[1,i] <- max(log_r[,i]) 
}
MaxR[1,14] <- max(Return_O[178:2173])

# Max Drawdown
MaxDown <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  MaxDown[i] <- maxDrawdown(log_r[178:2173,i])
}
MaxDown[1,14] <- maxDrawdown(Return_O[178:2173])

# Volatility
Volt <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Volt[1,i] <- sd(log_r[178:2173,i])*sqrt(250)
}
Volt[1,14] <- sd(Return_O[178:2173])*sqrt(250)

# Sharpe Ratio
Sharpe <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Sharpe[1,i] <- sharpe(log_r[178:2173,i],geometric.mean(1+log_r[178:2173,6])-1,scale=sqrt(250))
}
Sharpe[1,14] <- sharpe(Return_O[178:2173],geometric.mean(1+log_r[178:2173,6])-1,scale=sqrt(250))

# Sortino Ratio
Sortino <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Sortino[1,i] <- SortinoRatio(log_r[178:2173,i],0)
}
Sortino[1,14] <- SortinoRatio(Return_O[178:2173],0)

# Skewness
Skewness <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Skewness[1,i] <- skewness(log_r[178:2173,i]) 
}
Skewness[1,14] <- skewness(Return_O[178:2173])

# Kurtosis
Kurtosis <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  Kurtosis[1,i] <- kurtosis(log_r[178:2173,i]) 
}
Kurtosis[1,14] <- kurtosis(Return_O[178:2173])

# Modified VaR
VaR <- matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  VaR[1,i] <- VaR(log_r[178:2173,i], p = 0.95, method = c("modified")) 
}
VaR[1,14] <- VaR(Return_O[178:2173], p = 0.95, method = c("modified"))

# CVaR
CVaR <-  matrix(0,nrow=1,ncol=14)
for (i in 1:13){
  CVaR[1,i] <- ES(log_r[178:2173,i], p = 0.95, method = c("modified"))
}
CVaR[1,14] <- ES(Return_O[178:2173], p = 0.95, method = c("modified"))

# Combine all the solution
Analytics <- matrix(rbind(Geometric,MinR,MaxR,MaxDown,Volt,Sharpe,Sortino,Skewness,Kurtosis,
                          VaR,CVaR),nrow=11,ncol=14)
rownames(Analytics) <- c("Daily Mean Geometric Return", "Daily Min Return", "Daily Max Return", 
                         "Max Drawdown", "Volatility", "Sharpe Ratio","Sortino Ratio","Skewness",
                         "Kurtosis","Modified VaR","CVaR")
colnames(Analytics) <- c(sym_bols,"Optimal Portfolio")

Analytics1 <- Analytics
########
# Analytics1 <- Analytics
# Analytics2 <- Analytics
# Analytics3 <- Analytics
# Analytics4 <- Analytics
# write.zoo(Analytics1, file='Analytics1.csv', sep=",")
# write.zoo(Analytics2, file='Analytics2.csv', sep=",")
# write.zoo(Analytics3, file='Analytics3.csv', sep=",")
# write.zoo(Analytics4, file='Analytics4.csv', sep=",")
# save(Analytics1, Analytics2, Analytics3, Analytics4, file='etf_report.Rdata')

###########################
# Plot graphs for analytics
###################################################################################
# In order to test the sensitivity of our optimal portfolio, we separately
# calculate P&L based on the strategy with trend estimators (R60, C120) rebalancing
# once for each fixed time periods, which are 20 days and 5 days indicating 
# long term and short term. And comparing with the previous strategy with VIX
# constraint.

# Moreover, we change the target return to see the change of the behavior of
# our optimal portfolio. We use cumulated daily P&L as an indicator of the
# performance of our optimal portfolio.

# Therefore, we calculate cumulated daily P&L with the strategy of 20 days rebalance,
# 5 days rebalance and the strategy with 20% target return.
# save(PL_20, PL_5, PL_vix, PL_0.2, file='etf_PL.Rdata')
###################################################################################
# Plot The Behaviors Under Different Strategies
load("etf_PL.Rdata")
# PL_20 for 20 days rebalance, PL_5 for 5 days rebalance
# PL_vix rebalance depends on the index of VIX
# PL_0.2 change the target return from 15% into 20%
par(mfrow=c(2,1))
par(mar=c(3, 2, 1, 1), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
F <- cbind(PL_20,PL_vix, PL_5)
G <- cbind(PL_vix,PL_0.2)
matplot(F[179:2173], type="l", main="The Behaviors of (R60, C120) for 15% Target")
matplot(G[179:2173], type="l", main="The Behaviors of (R60, C120) for 15% And 20% Targets")
# write.zoo(F, file='F.csv', sep=",")
# write.zoo(G, file='G.csv', sep=",")
# We plot the combinations of cumulated daily P&L in EXCEL file.

#########################################################################
# Plot the behaviors of 13 assets under the strategy with the trend
# estimators of (R60, C120).
# get separate data
FXE <- xts(Sol1[,1], order.by = index(etf_series_ad))
EWJ <- xts(Sol1[,2], order.by = index(etf_series_ad))
GLD <- xts(Sol1[,3], order.by = index(etf_series_ad))
QQQ <- xts(Sol1[,4], order.by = index(etf_series_ad))
SPY <- xts(Sol1[,5], order.by = index(etf_series_ad))
SHV <- xts(Sol1[,6], order.by = index(etf_series_ad))
DBA <- xts(Sol1[,7], order.by = index(etf_series_ad))
USO <- xts(Sol1[,8], order.by = index(etf_series_ad))
XBI <- xts(Sol1[,9], order.by = index(etf_series_ad))
ILF <- xts(Sol1[,10], order.by = index(etf_series_ad))
GAF <- xts(Sol1[,11], order.by = index(etf_series_ad))
EPP <- xts(Sol1[,12], order.by = index(etf_series_ad))
FEZ <- xts(Sol1[,13], order.by = index(etf_series_ad))

# plot one ETF to view the behavior of the underlying strategy
par(mfrow=c(4,3))
par(mar=c(1, 1, 1, 1), mgp=c(0, 0, 0),oma=c(0.2,0.2,0.2,0.2), cex.lab=0.8, cex.axis=0.8,
    cex.main=0.8, cex.sub=0.8, xaxs="i")
# Add vertical lines indicating 2008 crisis period
cross_es1 <- as.POSIXct("2008-01-02")
cross_es2 <- as.POSIXct("2009-12-31")
plot(FXE, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(EWJ, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(GLD, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(QQQ, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(SPY, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(DBA, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(USO, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(XBI, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(ILF, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(GAF, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(EPP, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(FEZ, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")

# consider riskless security or not
# We plot the weights set of SHV under our four optimal strategies
SHV1 <- SHV
SHV2 <- xts(Sol2[,6], order.by = index(etf_series_ad))
SHV3 <- xts(Sol3[,6], order.by = index(etf_series_ad))
SHV4 <- xts(Sol4[,6], order.by = index(etf_series_ad))

par(mfrow=c(2,2))
par(mar=c(3, 2, 1, 1), mgp=c(2, 1, 0), cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.5)
plot(SHV1, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(SHV2, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(SHV3, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
plot(SHV4, type = "l")
abline(v=cross_es1, col="red", lty="dashed")
abline(v=cross_es2, col="red", lty="dashed")
# save(SHV1, SHV2, SHV3, SHV4, file='etf_shv.Rdata')
#####################################################