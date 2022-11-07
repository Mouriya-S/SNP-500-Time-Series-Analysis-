
library("forecast") #load forecast library

###################Load data#####################
#Load the data using the Import DataSet dropdown menu in the environment tab (top right)
data <-read.csv("C:/Users/msrini9/Desktop/IE 594/SNP 500.csv",header=TRUE)
  

###########Extract time series and plot#########

dataFreq= 12 #Data frequency of time series. Set to 12 for monthly  data
startEntry= c(1992,1) #Time stamp of first entry in time series e.g. c(1992,1) implies first month of 1992 if data frequency equals 12

#create a time series  
dataTS <- ts(data$Close, frequency=dataFreq, 
                 start=startEntry)

#plot time series.
plot(dataTS,main = "SNP 500 CLosing Price Monthly",
    xlab="Monthly",ylab="Closing Price") 

###########Decompose time series and plot#########
#Decompose methods return an object containing

tsDecomp <- decompose(dataTS, type="additive") #classical decomposition. Can be set to additive
plot(tsDecomp) #plot decomposed time series
summary(tsDecomp)


###########Prepare time series for forecasting#########
###We partition the time series into a training set for forecasting and a test set to evaluate accuracy####
trainSetStart= c(1992,1) #training set start location in time series (typically the first entry)
trainSetEnd= c(2019,4) #training set end location in time series (typically covers 90% of time series)
testSetStart= c(2019,5) #test set start location in time series (typically location of entry after training set ends)
testSetEnd= c(2022,4) #test set end location in time series (typically end of time series)

#extract training set  
Price_Train <- window(dataTS,start=trainSetStart,end=trainSetEnd) 

#extract test set
Price_Test <- window(dataTS,start=testSetStart,end=testSetEnd) #extract test set

###########Forecast#########
numForcPeriods = 36 #number of periods to forecast in the future 
    
  
HWForcModel <- HoltWinters(Price_Train,seasonal="additive") #Train Holt-Winters forecasting model. for additive 
HWForecast <- forecast(HWForcModel, h=numForcPeriods) #Forecast using Holt-Winters model trained in previous step

HWForecast
plot(HWForecast, main="Plot of training data, testing data, and forecast with 80% and 95% prediction intervals",xlab="Year.Montly", 
     ylab="Closing Price") #plot the training data, and forecast with prediction intervals
lines(Price_Test,col=2) #add the testing data line to plot
legend("topleft", lty=1, col=c(1,4,2), 
legend=c("Training Data","Forecast","Testing Data")) #create plot legend

###########Analyze forecasting error#########

error = HWForecast$mean-Price_Test   #difference between forecast and actual demand
summary(error)
mean(error^2)
AD=abs(error) #absolute value of error
AD


#Creating empty vectors to store errors
MSE <- matrix(AD, nrow = numForcPeriods, ncol = 1)
MAD <- matrix(AD, nrow = numForcPeriods, ncol = 1)
MAPE <- matrix(AD, nrow = numForcPeriods, ncol = 1)
bias <- matrix(AD, nrow = numForcPeriods, ncol = 1)
TS <- matrix(AD, nrow = numForcPeriods, ncol = 1)

#Labeling  columns of matrices using name of error

colnames(MSE) <- "MSE"
colnames(MAD) <- "MAD"
colnames(MAPE) <- "MAPE"
colnames(bias) <- "bias"
colnames(TS) <- "TS"

#compute errors 

for(t in 1:numForcPeriods){
  MSE[t] <- mean(error[1:t]*error[1:t])
  MAD[t] <- mean(AD[1:t])
  MAPE[t] <- mean(100*abs(error[1:t]/Price_Test[1:t]))
  bias[t] <- sum(error[1:t])
  TS[t]= bias[t]/MAD[t]
}

#combined vectors into a dataframe, also appending year and month information in the first two columns

error_Meas <- data.frame(floor(time(error)),cycle(error),Price_Test,HWForecast$mean,error,AD,MSE,MAD,MAPE,bias,TS)
colnames(error_Meas)[1] <- "Year"
colnames(error_Meas)[2] <- "Month"
colnames(error_Meas)[3] <- "Actual Price"
colnames(error_Meas)[4] <- "Forecast"
error_Meas
