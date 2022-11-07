## Defining % of sample to withhold for model verification
train_percent = 0.95

# red is for training data, blue is for full data sets

#Library Function for ARMA Forecasting

library(nlme)
library(lmtest)
library(stats)
library(strucchange)
library(ggplot2)
library(TSSS)
library(forecast)

#Data Entry 
data <- read.csv("C:/Users/msrini9/Desktop/IE 594/SNP 500.csv",header=TRUE)
data_title <- "S&P 500 Price, Monthly"
log_title <- "log(S&P 500 Price, Monthly)"
start_entry<-c(0,1); dataFreq = 12
data_ts <- ts(log(data$Close[1:(nrow(data)*train_percent)]), start=start_entry, freq=dataFreq)
data_ts_full <-ts(log(data$Close), start=start_entry, freq=dataFreq)

## Start Analysis
# Plot the raw data
par(mfrow=c(2,2))
plot(exp(data_ts), ylab="Closing Price", main=paste(c(data_title, "Training Data")),col="red")
plot(data_ts, ylab="Log Closing Price", main=paste(c(log_title, "Training Data")), col="red")
plot(exp(data_ts_full), ylab="Closing Price", main=paste(c(data_title, "Full Data Set")), col="blue")
plot(data_ts_full, ylab="Log Closing Price", main=paste(c(log_title,"Full Data Set")), col="blue")

## Plot Seasonal/Trend/Remainder to and see % of variance in each

plot(kitty <-stl(data_ts, s.window=12), main=paste(c("STL Decomposition", data_title)),col="red")
summary(kitty)

## Determine Polynomial Trend Order Trend
trend=time(data_ts)
trendf=time(data_ts_full)
sin_cos = cbind(sin(2*pi*trend/12), cos(2*pi*trend/12),sin(2*pi*trend/6),cos(2*pi*trend/6),
                sin(2*pi*trend/4), cos(2*pi*trend/4),sin(2*pi*trend/3), cos(2*pi*trend/3))
trend2 = trend^2; trend3 = trend^3; trend4 = trend^4; trend5 = trend^5;
trend6 = trend^6; trend7 = trend^7; trend8 = trend^8; trend9 = trend^9;
trend10 = trend^10

trend_model_1 = gls(data_ts ~ trend + sin_cos)
trend_model_2 = gls(data_ts ~ trend + trend2 + sin_cos)
trend_model_3 = gls(data_ts ~ trend + trend2 + trend3 + sin_cos)
trend_model_4 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + sin_cos)
trend_model_5 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 + sin_cos)
trend_model_6 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 +
                      trend6 + sin_cos)
trend_model_7 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 +
                      trend6 + trend7 + sin_cos)
trend_model_8 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 +
                      trend6 + trend7 + trend8 + sin_cos)
trend_model_9 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 +
                      trend6 + trend7 + trend8 + trend9 + sin_cos)
trend_model_10 = gls(data_ts ~ trend + trend2 + trend3 + trend4 + trend5 +
                       trend6 + trend7 + trend8 + trend9 + trend10 + sin_cos)

#Finding F Values of the Polynomial order trend

F_2_1 = ((sum(trend_model_1$residuals^2)-sum(trend_model_2$residuals^2))/1) /
  (sum(trend_model_2$residuals^2)/(length(data_ts)-length(trend_model_2$coefficients)))
F_3_2 = ((sum(trend_model_2$residuals^2)-sum(trend_model_3$residuals^2))/1) /
  (sum(trend_model_3$residuals^2)/(length(data_ts)-length(trend_model_3$coefficients)))
F_4_3 = ((sum(trend_model_3$residuals^2)-sum(trend_model_4$residuals^2))/1) /
  (sum(trend_model_4$residuals^2)/(length(data_ts)-length(trend_model_4$coefficients)))
F_5_4 = ((sum(trend_model_4$residuals^2)-sum(trend_model_5$residuals^2))/1) /
  (sum(trend_model_5$residuals^2)/(length(data_ts)-length(trend_model_5$coefficients)))
F_6_5 = ((sum(trend_model_5$residuals^2)-sum(trend_model_6$residuals^2))/1) /
  (sum(trend_model_6$residuals^2)/(length(data_ts)-length(trend_model_6$coefficients)))
F_7_6 = ((sum(trend_model_6$residuals^2)-sum(trend_model_7$residuals^2))/1) /
  (sum(trend_model_7$residuals^2)/(length(data_ts)-length(trend_model_7$coefficients)))
F_8_7 = ((sum(trend_model_7$residuals^2)-sum(trend_model_8$residuals^2))/1) /
  (sum(trend_model_8$residuals^2)/(length(data_ts)-length(trend_model_8$coefficients)))
F_9_8 = ((sum(trend_model_8$residuals^2)-sum(trend_model_9$residuals^2))/1) /
  (sum(trend_model_9$residuals^2)/(length(data_ts)-length(trend_model_9$coefficients)))
F_10_9 = ((sum(trend_model_9$residuals^2)-sum(trend_model_10$residuals^2))/1) /
  (sum(trend_model_10$residuals^2)/(length(data_ts)-length(trend_model_10$coefficients)))


F_2_1
1-pf(F_2_1,1,length(trend)-length(trend_model_2$coefficients))

F_3_2
1-pf(F_3_2,1,length(trend)-length(trend_model_3$coefficients))

F_4_3
1-pf(F_4_3,1,length(trend)-length(trend_model_4$coefficients))

F_5_4
1-pf(F_5_4,1,length(trend)-length(trend_model_5$coefficients))

F_6_5
1-pf(F_6_5,1,length(trend)-length(trend_model_6$coefficients))

F_7_6
1-pf(F_7_6,1,length(trend)-length(trend_model_7$coefficients))

F_8_7
1-pf(F_8_7,1,length(trend)-length(trend_model_8$coefficients))

F_9_8
1-pf(F_9_8,1,length(trend)-length(trend_model_9$coefficients))

F_10_9
1-pf(F_10_9,1,length(trend)-length(trend_model_10$coefficients))

# Reduce the model size to the smallest appropriate model till trend 3

trend.m = gls(data_ts ~ trend + trend2 + trend3  + sin_cos)
summary(trend.m)

#removing sin and cosine periodicity with high p-values

remove.1 = gls(data_ts ~ trend + trend2 + trend3 +
                 sin(2*pi*trend/12) +cos(2*pi*trend/4) +sin(2*pi*trend/6) +
                 cos(2*pi*trend/6)  + sin(2*pi*trend/3) )
summary(remove.1)

# remove.1 is the best model

trend_residuals = remove.1$residuals
trend_model = remove.1
par(mfrow=c(1,1))
ts.plot(data_ts, trend_model$fitted, ylab="Log Closing Price",
        main=paste(c(log_title, "vs Forecast")), col=c("red","blue"),lwd=1)

#Add additional regressors as needed for model

trend_regressors = cbind(trend, trend2, trend3,
                         sin(2*pi*trend/12), 
                         sin(2*pi*trend/6), cos(2*pi*trend/6),
                         cos(2*pi*trend/4),sin(2*pi*trend/3)
)
RSS_trend = sum(trend_residuals^2)
test_range = c((length(data_ts)+1): length(data_ts_full))

test_xreg = cbind(trendf[test_range], trendf[test_range]^2, trendf[test_range]^3,
                  sin(2*pi*trendf[test_range]/12),
                  sin(2*pi*trendf[test_range]/6), cos(2*pi*trendf[test_range]/6), cos(2*pi*trendf[test_range]/4),
                  sin(2*pi*trendf[test_range]/3))


# Get model details and trend adjusted residual plots

summary(trend_model)
par(mfrow=c(1,1))
plot(trend_residuals, pch=20, col="red", main=paste(data_title,": Trend Adjusted Residuals"),lwd=1)
abline(h=0, col="black", lty="dotted")


## Plot acf and pcf of residuals

par(mfrow=c(1,2))
acf(trend_residuals,(length(data_ts))/20, col="red") #use 5% of data for autocovariance 
pacf(trend_residuals,(length(data_ts))/20, col="red")#use 5% of data for autocovariance


## Modeling Procedure ARMA (2n,2n-1)

# AR(1)
arma_1_0 = arima(trend_residuals, order=c(1,0,0))
arma_1_0
RSS_1_0 = arma_1_0$sigma2*length(arma_1_0$residuals)

# ARMA(2,1)
n=1; p=2*n; d=0; q=2*n-1
arma_2_1 = arima(trend_residuals, order=c(p,d,q))
arma_2_1
RSS_2_1 = arma_2_1$sigma2*length(arma_2_1$residuals)
RSS_2_1
F_stat = ((RSS_1_0 - RSS_2_1)/2)/(RSS_2_1/(length(data_ts)-(4*n)))
1-pf(F_stat,4,length(trend_residuals)-4*n)

# ARMA(4,3)
n=2; p=2*n; d=0; q=2*n-1
arma_4_3 = arima(trend_residuals, order=c(p,d,q), method="CSS")
arma_4_3
RSS_4_3 = arma_4_3$sigma2*length(arma_4_3$residuals)
RSS_4_3
F_stat_1 = ((RSS_2_1 - RSS_4_3)/4)/(RSS_4_3/(length(data_ts)-(4*n)))
F_stat_1
1-pf(F_stat_1,4,length(trend_residuals)-4*n)

# ARMA(6,5)
n=3; p=2*n; d=0; q=2*n-1
arma_6_5 = arima(trend_residuals, order=c(p,d,q), method="CSS")
arma_6_5
RSS_6_5 = arma_6_5$sigma2*length(arma_6_5$residuals)
RSS_6_5
F_stat_2 = ((RSS_4_3 - RSS_6_5)/4)/(RSS_6_5/(length(data_ts)-(4*n)))
F_stat_2
1-pf(F_stat_2,4,length(trend_residuals)-4*n)

#integrated model

trend_regressors = cbind(trend, trend2, trend3,
                         sin(2*pi*trend/12), sin(2*pi*trend/6), cos(2*pi*trend/6),
                         cos(2*pi*trend/4), sin(2*pi*trend/3))

test_xreg = cbind(trendf[test_range], trendf[test_range]^2,trendf[test_range]^3, 
                  sin(2*pi*trendf[test_range]/12),
                  sin(2*pi*trendf[test_range]/6), cos(2*pi*trendf[test_range]/6),
                  cos(2*pi*trendf[test_range]/4),
                  sin(2*pi*trendf[test_range]/3))

integrated_model = arima(data_ts, order=c(4,0,3), xreg = trend_regressors)
integrated_model

## Print coefs, std errors and calculate p-values for each of the coefficients
variable_matrix = data.frame(estimate=integrated_model$coef,
                             std_err=rep(0,length(integrated_model$coef)), p_value=rep(0,length(integrated_model$coef)))

for(i in 1:length(integrated_model$coef)){
  variable_matrix[i,2] = sqrt(integrated_model$var.coef[i,i])
  variable_matrix[i,3] = dt(integrated_model$coef[i]/sqrt(integrated_model$var.coef[i,i]),
                            length(data_ts)-length(integrated_model$coef))
}
variable_matrix

# integrated_model residual plots
par(mfrow=c(2,1))
plot(integrated_model$residuals, col="red", main="Integrated Model Residuals", )
pacf(integrated_model$residuals, main="PACF: Integrated Model Residuals", col="red")

## Forecasting graph
joint_pred = predict(integrated_model, n.ahead=ceiling(length(data_ts_full)*(1-train_percent)), newxreg=test_xreg)
U = joint_pred$pred + 2*joint_pred$se
L = joint_pred$pred - 2*joint_pred$se

par(mfrow=c(1,1))
ts.plot(data_ts_full, joint_pred$pred, col=c("red","blue"),
        main=paste(c(log_title, "vs Forecast with 95% Confidence Interval of", train_percent*100 ,"percent Training dataset")))
lines(U, col="Black", lty="dashed")
lines(L, col="Black", lty="dashed")


## Diagnostics on residuals
Box.test(integrated_model$residuals, 13, type="Box-Pierce") #correlated residuals? 
Box.test(integrated_model$residuals, 13, type="Ljung")  #correlated residuals?
shapiro.test(integrated_model$residuals) # normality of residuals?


#close-up on forecast period
ts.plot(data_ts_full, joint_pred$pred, col=c("red","blue"),
        main=paste(c(log_title, "vs  Forecast with 95% Confidence Interval of ",train_percent*100 ,"percent Training dataset")),
        xlim=c(length(data_ts_full)*train_percent/12,length(data_ts_full)/12))
lines(U, col="Black", lty="dashed")
lines(L, col="Black", lty="dashed")
abline(v=(length(data_ts)/12), col="black", lty="dotted")

#### MSE Calculation ###
mse_data = ts(data_ts_full[test_range], start=c((length(data_ts))/12), freq=12)
mse_data
joint_pred$pred
error = mse_data - joint_pred$pred
mean(error^2)

## Integrated_model residuals
par(mfrow=c(1,1))
plot(integrated_model$residuals, pch=20, col="red", main=paste(data_title,": Int. Model Resids"))
abline(h=0, col="black", lty="dotted")

## calculate roots. model is in increasing exponent order. Must hand input.
polyroot(c(1,-0.5939757880,1.1289269396,0.6167066446, -0.2646756250))

## Future Forecasting graph
joint_pred_Future = predict(integrated_model, n.ahead=ceiling(length(data_ts_full)*(1-(train_percent-.3))), newxreg=test_xreg)
U = joint_pred_Future$pred + 2*joint_pred_Future$se
L = joint_pred_Future$pred - 2*joint_pred_Future$se

par(mfrow=c(1,1))
ts.plot(data_ts_full, joint_pred_Future$pred, col=c("red","blue"),
        main=paste(c(log_title, "vs Forecast with 95% Confidence Interval of", train_percent*100 ,"percent Training dataset")))
lines(U, col="Black", lty="dashed")
lines(L, col="Black", lty="dashed")

#close-up on forecast period
ts.plot(data_ts_full, joint_pred_Future$pred, col=c("red","blue"),
        main=paste(c(log_title, "vs Future Forecast with 95% Confidence Interval of ",train_percent*100 ,"percent Training dataset")),
        xlim=c(length(data_ts_full)*train_percent/11.5,length(data_ts_full)/11.5))
lines(U, col="Black", lty="dashed")
lines(L, col="Black", lty="dashed")
abline(v=(length(data_ts)/12), col="black", lty="dotted")
exp(joint_pred_Future$pred) #Predicted value for the future
