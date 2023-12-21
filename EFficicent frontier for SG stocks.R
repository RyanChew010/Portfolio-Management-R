library(DEoptim)
library(Rglpk)
library(pso)
library(GenSA)
library(corpcor)
library(testthat)
library(nloptr)
library(MASS)
library(robustbase)
library(fPortfolio)
library(timeSeries)
library(quantmod)
library(dplyr)
library(PerformanceAnalytics)
library(ggplot2)
library(TTR)

# Genting, ST, YZJ, ShengShiong, SATS
s1 <- c("G13.SI", "S63.SI", "BS6.SI", "OV8.SI", "S58.SI")

# Download prices and create returns from adjusted prices 
data1 <- lapply(s1, FUN = function(x){
  ROC(Ad(getSymbols(x, from = "2022-12-20", to = "2023-12-21", auto.assign = F)),
      type = "discrete") * 100
})
data1
# convert to data frame 
ret1 <- as.data.frame(do.call(merge, data1))
#change column names 
colnames(ret1) <- gsub(".SI.Adjusted", "", colnames(ret1))
#remove the first row of missing values 
ret1 <- ret1[-1, ]
# add dates column
ret1 <- data.frame(Date = as.Date(row.names(ret1)), ret1)
row.names(ret1) <- NULL

#Saving the dataframe (Only for reproducibility)
#saveRDS(ret1, file = "data/port_ret.Rds")

#---------------------------------------------------------------------
# Plotting the data 
library(ggplot2)
library(tidyr)
library(pander)
library(ggthemes)

# overview
pander(tail(ret1), split.table = Inf)

#Convert to long 
ret_long <- pivot_longer(ret1, cols = -c(Date), values_to = "Return", names_to = 
                           "stock")

# Plotting 
port_p1 <- ggplot(ret_long, aes(Date, Return, colour = stock)) + 
  geom_path(stat = "identity") + facet_grid(stock ~.) + theme_minimal() + 
  labs(x = "Date", y = "Returns")
port_p1

np1 <- 200 
ret2 <- ret1[, -1]
mu1 <- colMeans(ret2)
na1 <- ncol(ret2)
varc1<- cov(ret2)

riskp1 <- NULL
retp1 <- NULL 

for (i in 1:np1) {
  w = diff(c(0, sort(runif(na1 - 1)), 1)) #Random weights
  r1 = t(w) %*% mu1 #Matrix multiplication
  sd1 = t(w) %*% varc1 %*% w 
  retp1 = rbind(retp1, r1)
  riskp1 = rbind(riskp1, sd1)
  
}

d_p1 <- data.frame(Ret = retp1, Risk = riskp1)
d_p1

plot(d_p1$Risk, d_p1$Ret, xlab = "risk", ylab = "return", 
     main = "Frontier Portfolios", col = "black")

#GGPLOT
#First Layer
p1 <- ggplot(d_p1, aes(Risk, Ret, colour = Ret))
#Scatterplot 
p1 <- p1 + geom_point()
# Scatter with density and identified portfolio risk return 
# Highest, lowest return and min risk 
p1 + geom_point() + geom_hline(yintercept = c(max(d_p1$Ret), median(d_p1$Ret), 
                                              min(d_p1$Ret)), colour = 
                                 c("darkgreen", "darkgray", "darkred"), linewidth = 1) + 
  geom_vline(xintercept = d_p1[(d_p1$Risk == min(d_p1$Risk)), ][, 2]) + 
  labs(colour = "Portfolio Return", x = "Portfolio Risk", y = "Portfolio Return", 
       title = "Random Feasible Portfolios") + theme_bw()

library(ROI)
library(ROI.plugin.glpk)
library(ROI.plugin.symphony)
library(ROI.plugin.quadprog)
library(PortfolioAnalytics)
#Initialize with asset names uses time series data 
data_p2 <- zoo(ret1[, -1], order.by = as.Date(ret1$Date))
# create specification 
port <- portfolio.spec(assets = c(colnames(data_p2)))
#Add long only constraint 
port <- add.constraint(portfolio = port, type = "long_only")
# add full investment constraint 
port <- add.constraint(portfolio = port, type = "weight_sum", min_sum = 0.9,
                       max_sum = 1.01)
# add objective: minimize risk 
port <- add.objective(portfolio = port, type = "risk", name = "StdDev")
# add objective: maximize return 
port <- add.objective(portfolio = port, type = "return", name = "mean")

#Optimize portfolios
rand_p <- optimize.portfolio(R = data_p2, portfolio = port, 
                             optimize_method = "random", trace = TRUE, 
                             search_size = 1000)
# plot (Also plots equally weighted portfolio) 
chart.RiskReward(rand_p, risk.col = "StdDev", return.col = "mean", 
                 chart.assets = T) 
chart.Weights(rand_p)

# Efficient frontier
port_msd <- add.objective(portfolio = port, type = "risk", name = "StdDev")
minvar1 <- optimize.portfolio(R = data_p2, portfolio = port_msd,
                              optimize_method = "ROI", trace = TRUE) #Error??
minvar1

#fPortfolio 
library(fPortfolio)
data_p2 <- as.timeSeries(data_p2)
pspec <- portfolio.spec(assets = c(colnames(data_p2)))
# Random portfolios for efficient frontier

eff_front2 <- portfolioFrontier(data = data_p2, constraints = "LongOnly")
plot(eff_front2, c(1,2,4,5,6))

#Another function (allows to add different points and lines from other ptf)
tailoredFrontierPlot(eff_front2, sharpeRatio = F, risk = "Sigma")

#Weights Efficient frontier 
weightsPlot(eff_front2)

# Min Variance and Portfolio for given level of return (L0ngOnly)
(minvar2 = minvariancePortfolio(data_p2))

#Sharpe Ratio 
pspec3 = portfolioSpec()
setRiskFreeRate(pspec3) = 0.0272
setNFrontierPoints(pspec3) = 100
# create efficient frontier
eff_front4 = portfolioFrontier(data_p2, spec = pspec3, constraints = "LongOnly")

# find the tangency port
tgport1 = tangencyPortfolio(data = data_p2, spec = pspec3, constraints = "LongOnly")
tgport1

# create frontier plot
frontierPlot(eff_front4, pch = 1, auto = F, xlim = c(0, 2.0), ylim = c(0, 0.1))
# add min variance point 
minvariancePoints(object = eff_front4, auto = F, col = "red", pch = 20)
#add tangency point 
tangencyPoints(object = eff_front4, col = "blue", pch = 20)
tangencyLines(object = tgport1, col = "darkblue")







