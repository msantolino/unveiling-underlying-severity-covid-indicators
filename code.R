library(forecast)
library(odpc)
library(gdpc)
library(tvReg)

#### Read the daily positive cases, hospitalizations, ICU admissions and deaths for each group age
indicador <- read.csv("datos_edades.csv") #database with information from 2020-05-11 to 2022-03-27
indicador <- indicador[,-1]
indicador$fecha<-as.Date(indicador$fecha,format="%Y-%m-%d") #date
indicador$grupo<-as.factor(indicador$grupo) #age group

#### Select the group with all the ages (Todos = All)
BDcom<-indicador[indicador$grupo=="Todos",]
# Spanish names of the variables date, positive cases, hospital, ICU, deaths and group
names(BDcom)<-c("fecha","positivo","hospital","uci","muerte","grupo")

#### Remove the weekly seasonality
# positive cases
positivo <- ts(BDcom$positivo, BDcom$fecha, frequency=7)
positivoC<-seasadj(stl(positivo, s.window = 7, t.window = 50, t.jump = 1,robust=T)) 

# hospital admissions
hospital <- ts(BDcom$hospital, BDcom$fecha, frequency=7)
hospitalC<-seasadj(stl(hospital, s.window = 7, t.window = 50, t.jump = 1,robust=T))

# ICU admissions
uci <- ts(BDcom$uci, BDcom$fecha, frequency=7)
uciC<-seasadj(stl(uci, s.window=35,robust=T))

# deaths
muerte <- ts(BDcom$muerte, BDcom$fecha, frequency=7)
muerteC<-seasadj(stl(muerte, s.window = 7, t.window = 50, t.jump = 1,robust=T))


#### Kernel smoothing

positivoK<-ksmooth(BDcom$fecha,as.vector(positivoC), "normal", bandwidth = 12)$y # positive cases
hospitalK<-ksmooth(BDcom$fecha,as.vector(hospitalC), "normal", bandwidth = 12)$y # hospitalizations
uciK<-ksmooth(BDcom$fecha,as.vector(uciC), "normal", bandwidth = 12)$y # ICU admissions
muerteK<-ksmooth(BDcom$fecha,as.vector(muerteC), "normal", bandwidth = 12)$y # deaths


#### Plot Figure 1

f1 <- seq(from = BDcom$fecha[1], to = BDcom$fecha[686], length.out = 20)
plot(as.vector(muerte)~BDcom$fecha,type="l", ylim=c(-1,2700), ylab="Value", xlab="", xaxt = "n")
lines(as.vector(positivo)/100~BDcom$fecha, col="darkgrey")
lines(as.vector(hospital)~BDcom$fecha, col="blue3")
lines(as.vector(uci)~BDcom$fecha, col="green")
lines(muerteK~BDcom$fecha, col="red", lwd = 0.7)
lines(hospitalK~BDcom$fecha, col="red", lwd = 0.7)
lines(uciK~BDcom$fecha, col="red", lwd = 0.7)
lines(positivoK/100~BDcom$fecha, col="red", lwd = 0.7)
axis(1, f1, format(f1, "%Y-%m-%d"), cex.axis = .7, las = 2)
mtext(side=1, text="Date", line=4)
legend("topleft", legend=c("Positives/100","Hospitalizations","ICU admissions","Deaths"),col=c("darkgrey","blue3","green","black"),lty=1, cex=0.8, inset = c(0.03, 0.075))

#### Stationarity analysis of the corrected time series (after removing seasonality and Kernel smoothing)

adf.test(positivoK,k=1) # stationary positive cases
adf.test(hospitalK,k=1) # stationary hospitalizations
adf.test(uciK,k=1) # stationary ICU admissions
adf.test(muerteK,k=1) # stationary deaths

#### Correlations between positive cases, hospitalizations, ICU admissions and deaths
cor(cbind(as.vector(positivoK),as.vector(hospitalK), as.vector(uciK),as.vector(muerteK)))

#### One-sided dynamic principal components

# Standardize the data (series with mean = 1)
expST<-(positivoK/mean(positivoK)) # positive cases
exhST<-(hospitalK/mean(hospitalK)) # hospitalizations
exuST<-(uciK/mean(uciK)) # ICU
exmST<-(muerteK/mean(muerteK)) # deaths

# We include the deaths, hospitalizations and ICU admissions to compute its principal components
covidst3<-cbind(exmST,exhST,exuST)

# We check the number of lags to include in the model (using two different criteria and imposing 1 component)
autosel <- cv.odpc(covidst3, h=1, k_list = c(1:10), max_num_comp = 1, window_size = 2, ncores_k = 2, ncores_w = 2)
autosel
autosel2 <- crit.odpc(covidst3,k_list = c(1:10),  max_num_comp = 1)
autosel2

# We fit the model with 1 lag
fit <- odpc(covidst3, ks = c(1))
fit
fit[[1]]$B # deaths, hospitalizatiions and ICU
fit[[1]]$a

# explained variance by the component
expart <- 1 - fit[[1]]$mse/mean(apply(covidst3, 2, var))
expart

# Compare the results to the ones without using lags (static principal components)
fitCst<- gdpc(covidst3, k = 0)
fitCst

#### Plot Figure 2

fecha2 <- BDcom$fecha[2:686]
f2 <- seq(from = BDcom$fecha[2], to = BDcom$fecha[686], length.out = 20)
plot(fit[[1]]$f~fecha2,type="l",ylim=c(0,5), ylab="Value", xlab="", xaxt = "n",col="red",lwd=2)
lines(exmST[2:686]~fecha2,col="black")
lines(exhST[2:686]~fecha2,col="blue")
lines(exuST[2:686]~fecha2,col="green")
axis(1, f2, format(f2, "%Y-%m-%d"), cex.axis = .7, las = 2)
mtext(side=1, text="Date", line=4)
legend("topleft", legend=c("Hospitalizations","ICU admissions","Deaths","ODPC indicator"),col=c("blue","green","black","red"),lty=1, cex=0.8, inset = c(0.03, 0.075),lwd=c(1,1,1,2))

#### Reconstruction of the series (deaths, hospitalizations and ICU)
recons <- fitted(fit, num_comp = 1)

# Mean Squared Error of the reconstruction
mean((covidst3[3:686,1]-recons[,1])**2) # deaths
mean((covidst3[3:686,2]-recons[,2])**2) # hospitalizations
mean((covidst3[3:686,3]-recons[,3])**2) # ICU

#### One day ahead forecasts of the series

predict<-function(data,m){
  n<-nrow(data)
  A<-matrix(0, (n-m), 3)
  for(i in m:(n-1)){
    datai<-data[1:i,]
    fit <- odpc(datai, ks = c(1))
    k<-i+1-m
    A[k,]<-forecast.odpcs(fit, h = 1,Z=datai,add_residuals="TRUE")
  }
  return(A)
}

# at least we use 100 observations to fit the odpc and make the predictions
prediccionR <- predict(covidst3,m=100)


# Mean Squared Error of the predictions
mean((covidst3[101:686,1]-prediccionR[,1])**2) # deaths
mean((covidst3[101:686,2]-prediccionR[,2])**2) # hospitalizations
mean((covidst3[101:686,3]-prediccionR[,3])**2) # ICU

#### Plot Figure 3
par(mfcol = c(2, 3), oma=c(2.5,0,0,0))

fecha4 <- BDcom$fecha[3:686]
f4 <- seq(from = BDcom$fecha[3], to = BDcom$fecha[686], length.out = 15)
fecha5 <- BDcom$fecha[101:686]
f5 <- seq(from = BDcom$fecha[101], to = BDcom$fecha[686], length.out = 15)


plot(exhST[3:686]~fecha4,type="l", ylab = "Hospitalizations", xlab = "", xaxt = "n")
lines(recons[,2]~fecha4,col="blue", lty = 2)
axis(1, f4, format(f4, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(a)", line=6)

plot(covidst3[101:686,2]~fecha5,type="l", ylab = "Hospitalizations", xlab = "", xaxt = "n")
lines(prediccionR[,2]~fecha5,col="green", lty = 2)
axis(1, f5, format(f5, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(d)", line=6)

plot(exuST[3:686]~fecha4,type="l", ylab = "ICU admissions", xlab = "", xaxt = "n")
lines(recons[,3]~fecha4,col="blue", lty = 2)
axis(1, f4, format(f4, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(b)", line=6)

plot(covidst3[101:686,3]~fecha5,type="l", ylab = "ICU admissions", xlab = "", xaxt = "n")
lines(prediccionR[,3]~fecha5,col="green", lty = 2)
axis(1, f5, format(f5, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(e)", line=6)

plot(exmST[3:686]~fecha4,type="l", ylab = "Deaths", xlab = "", xaxt = "n")
lines(recons[,1]~fecha4,col="blue", lty = 2)
axis(1, f4, format(f4, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(c)", line=6)

plot(covidst3[101:686,1]~fecha5,type="l", ylab = "Deaths", xlab = "", xaxt = "n")
lines(prediccionR[,1]~fecha5,col="green", lty = 2)
axis(1, f5, format(f5, "%Y-%m-%d"), cex.axis = .9, las = 2)
mtext(side=1, text="(f)", line=6)

#### TVLM

# Correlations between the principal component (severity) and lagged values of the positive cases 
par(mfcol = c(1,1), oma=c(0.0,0,0,0))
print(ccf(expST[1:685],fit[[1]]$f,lag.max=20)) # maximum value at -8 lags (8 days earlier)

positivo<-expST[1:677] # positive cases 
gravedad<-fit[[1]]$f[9:685] # 8 days lagged principal component (severity indicator)

# Fit the model
model.tvLM <- tvLM(gravedad ~ 0+positivo, bw = 0.25)
summary(model.tvLM)

#### Plot Figure 4
fecha3 <- fecha2[9:685]
f3 <- seq(from = fecha3[1], to = fecha3[677], length.out = 20)
plot(model.tvLM$coeff~fecha3,type="l", ylim=c(0,2), ylab="Coefficient", xlab="", xaxt = "n")
axis(1, f3, format(f3, "%Y-%m-%d"), cex.axis = .7, las = 2)
mtext(side=1, text="Date", line=4)
abline(v=fecha3[300], col = "red",lwd=2, lty=2) # 15 March 2021, 4,6% fully-vaccination rate
abline(v=fecha3[500], col = "red",lwd=2, lty=2) # 1 October 2021, 82,9% fully-vaccination rate


#### Same procedure for the different group ages

#### 20-49

#### Select the group 20-49 ages 
BDcom20<-indicador[indicador$grupo=="20-49",]
# names of the variables in English are date, positive cases, hospital, ICU, deaths and group
names(BDcom20)<-c("fecha","positivo","hospital","uci","muerte","grupo")

#### Remove the weekly seasonality
# positive cases
positivo20 <- ts(BDcom20$positivo, BDcom20$fecha, frequency=7)
positivoC20<-seasadj(stl(positivo20, s.window = 7, t.window = 50, t.jump = 1,robust=T)) 

# hospital admissions
hospital20 <- ts(BDcom20$hospital, BDcom20$fecha, frequency=7)
hospitalC20<-seasadj(stl(hospital20, s.window = 7, t.window = 50, t.jump = 1,robust=T))

# ICU admissions
uci20 <- ts(BDcom20$uci, BDcom20$fecha, frequency=7)
uciC20<-seasadj(stl(uci20, s.window=35,robust=T))

# deaths
muerte20 <- ts(BDcom20$muerte, BDcom20$fecha, frequency=7)
muerteC20<-seasadj(stl(muerte20, s.window = 7, t.window = 50, t.jump = 1,robust=T))


#### Kernel smoothing

positivoK20<-ksmooth(BDcom20$fecha,as.vector(positivoC20), "normal", bandwidth = 12)$y # positive cases
hospitalK20<-ksmooth(BDcom20$fecha,as.vector(hospitalC20), "normal", bandwidth = 12)$y # hospitalizations
uciK20<-ksmooth(BDcom20$fecha,as.vector(uciC20), "normal", bandwidth = 12)$y # ICU admissions
muerteK20<-ksmooth(BDcom20$fecha,as.vector(muerteC20), "normal", bandwidth = 12)$y # deaths

#### One-sided dynamic principal components

# Standardize the data (divided by the mean of the whole population)
expST20<-(positivoK20/mean(positivoK)) # positive cases
exhST20<-(hospitalK20/mean(hospitalK)) # hospitalizations
exuST20<-(uciK20/mean(uciK)) # ICU
exmST20<-(muerteK20/mean(muerteK)) # deaths

covidst20<-cbind(exmST20,exhST20,exuST20)

# We fit the model with 1 lag
fit20 <- odpc(covidst20, ks = c(1))
fit20

#### TVLM

positivo20<-expST20[1:677] # positive cases 
gravedad20<-fit20[[1]]$f[9:685] # 8 days lagged principal component (severity indicator)

# Fit the model
model.tvLM20 <- tvLM(gravedad20 ~ 0+positivo20, bw = 0.25)
summary(model.tvLM20)

coef20 <- model.tvLM20$coeff

#### 50-69

#### Select the group 50-69 ages 
BDcom50<-indicador[indicador$grupo=="50-69",]
# names of the variables in English are date, positive cases, hospital, ICU, deaths and group
names(BDcom50)<-c("fecha","positivo","hospital","uci","muerte","grupo")

#### Remove the weekly seasonality
# positive cases
positivo50 <- ts(BDcom50$positivo, BDcom50$fecha, frequency=7)
positivoC50<-seasadj(stl(positivo50, s.window = 7, t.window = 50, t.jump = 1,robust=T)) 

# hospital admissions
hospital50 <- ts(BDcom50$hospital, BDcom50$fecha, frequency=7)
hospitalC50<-seasadj(stl(hospital50, s.window = 7, t.window = 50, t.jump = 1,robust=T))

# ICU admissions
uci50 <- ts(BDcom50$uci, BDcom50$fecha, frequency=7)
uciC50<-seasadj(stl(uci50, s.window=35,robust=T))

# deaths
muerte50 <- ts(BDcom50$muerte, BDcom50$fecha, frequency=7)
muerteC50<-seasadj(stl(muerte50, s.window = 7, t.window = 50, t.jump = 1,robust=T))


#### Kernel smoothing

positivoK50<-ksmooth(BDcom50$fecha,as.vector(positivoC50), "normal", bandwidth = 12)$y # positive cases
hospitalK50<-ksmooth(BDcom50$fecha,as.vector(hospitalC50), "normal", bandwidth = 12)$y # hospitalizations
uciK50<-ksmooth(BDcom50$fecha,as.vector(uciC50), "normal", bandwidth = 12)$y # ICU admissions
muerteK50<-ksmooth(BDcom50$fecha,as.vector(muerteC50), "normal", bandwidth = 12)$y # deaths

#### One-sided dynamic principal components

# Standardize the data (divided by the mean of the whole population)
expST50<-(positivoK50/mean(positivoK)) # positive cases
exhST50<-(hospitalK50/mean(hospitalK)) # hospitalizations
exuST50<-(uciK50/mean(uciK)) # ICU
exmST50<-(muerteK50/mean(muerteK)) # deaths

covidst50<-cbind(exmST50,exhST50,exuST50)

# We fit the model with 1 lag
fit50 <- odpc(covidst50, ks = c(1))
fit50

#### TVLM

positivo50<-expST50[1:677] # positive cases 
gravedad50<-fit50[[1]]$f[9:685] # 8 days lagged principal component (severity indicator)

# Fit the model
model.tvLM50 <- tvLM(gravedad50 ~ 0+positivo50, bw = 0.25)
summary(model.tvLM50)

coef50 <- model.tvLM50$coeff

#### 70-79

#### Select the group 70-79 ages 
BDcom70<-indicador[indicador$grupo=="70-79",]
# names of the variables in English are date, positive cases, hospital, ICU, deaths and group
names(BDcom70)<-c("fecha","positivo","hospital","uci","muerte","grupo")

#### Remove the weekly seasonality
# positive cases
positivo70 <- ts(BDcom70$positivo, BDcom70$fecha, frequency=7)
positivoC70<-seasadj(stl(positivo70, s.window = 7, t.window = 50, t.jump = 1,robust=T)) 

# hospital admissions
hospital70 <- ts(BDcom70$hospital, BDcom70$fecha, frequency=7)
hospitalC70<-seasadj(stl(hospital70, s.window = 7, t.window = 50, t.jump = 1,robust=T))

# ICU admissions
uci70 <- ts(BDcom70$uci, BDcom70$fecha, frequency=7)
uciC70<-seasadj(stl(uci70, s.window=35,robust=T))

# deaths
muerte70 <- ts(BDcom70$muerte, BDcom70$fecha, frequency=7)
muerteC70<-seasadj(stl(muerte70, s.window = 7, t.window = 50, t.jump = 1,robust=T))


#### Kernel smoothing

positivoK70<-ksmooth(BDcom70$fecha,as.vector(positivoC70), "normal", bandwidth = 12)$y # positive cases
hospitalK70<-ksmooth(BDcom70$fecha,as.vector(hospitalC70), "normal", bandwidth = 12)$y # hospitalizations
uciK70<-ksmooth(BDcom70$fecha,as.vector(uciC70), "normal", bandwidth = 12)$y # ICU admissions
muerteK70<-ksmooth(BDcom70$fecha,as.vector(muerteC70), "normal", bandwidth = 12)$y # deaths

#### One-sided dynamic principal components

# Standardize the data (divided by the mean of the whole population)
expST70<-(positivoK70/mean(positivoK)) # positive cases
exhST70<-(hospitalK70/mean(hospitalK)) # hospitalizations
exuST70<-(uciK70/mean(uciK)) # ICU
exmST70<-(muerteK70/mean(muerteK)) # deaths

covidst70<-cbind(exmST70,exhST70,exuST70)

# We fit the model with 1 lag
fit70 <- odpc(covidst70, ks = c(1))
fit70

#### TVLM

positivo70<-expST70[1:677] # positive cases 
gravedad70<-fit70[[1]]$f[9:685] # 8 days lagged principal component (severity indicator)

# Fit the model
model.tvLM70 <- tvLM(gravedad70 ~ 0+positivo70, bw = 0.25)
summary(model.tvLM70)

coef70 <- model.tvLM70$coeff

#### 80*

#### Select the group 80+ ages 
BDcom80<-indicador[indicador$grupo=="80+",]
# names of the variables in English are date, positive cases, hospital, ICU, deaths and group
names(BDcom80)<-c("fecha","positivo","hospital","uci","muerte","grupo")

#### Remove the weekly seasonality
# positive cases
positivo80 <- ts(BDcom80$positivo, BDcom80$fecha, frequency=7)
positivoC80<-seasadj(stl(positivo80, s.window = 7, t.window = 50, t.jump = 1,robust=T)) 

# hospital admissions
hospital80 <- ts(BDcom80$hospital, BDcom80$fecha, frequency=7)
hospitalC80<-seasadj(stl(hospital80, s.window = 7, t.window = 50, t.jump = 1,robust=T))

# ICU admissions
uci80 <- ts(BDcom80$uci, BDcom80$fecha, frequency=7)
uciC80<-seasadj(stl(uci80, s.window=35,robust=T))

# deaths
muerte80 <- ts(BDcom80$muerte, BDcom80$fecha, frequency=7)
muerteC80<-seasadj(stl(muerte80, s.window = 7, t.window = 50, t.jump = 1,robust=T))


#### Kernel smoothing

positivoK80<-ksmooth(BDcom80$fecha,as.vector(positivoC80), "normal", bandwidth = 12)$y # positive cases
hospitalK80<-ksmooth(BDcom80$fecha,as.vector(hospitalC80), "normal", bandwidth = 12)$y # hospitalizations
uciK80<-ksmooth(BDcom80$fecha,as.vector(uciC80), "normal", bandwidth = 12)$y # ICU admissions
muerteK80<-ksmooth(BDcom80$fecha,as.vector(muerteC80), "normal", bandwidth = 12)$y # deaths

#### One-sided dynamic principal components

# Standardize the data (divided by the mean of the whole population)
expST80<-(positivoK80/mean(positivoK)) # positive cases
exhST80<-(hospitalK80/mean(hospitalK)) # hospitalizations
exuST80<-(uciK80/mean(uciK)) # ICU
exmST80<-(muerteK80/mean(muerteK)) # deaths

covidst80<-cbind(exmST80,exhST80,exuST80)

# We fit the model with 1 lag
fit80 <- odpc(covidst80, ks = c(1))
fit80

#### TVLM

positivo80<-expST80[1:677] # positive cases 
gravedad80<-fit80[[1]]$f[9:685] # 8 days lagged principal component (severity indicator)

# Fit the model
model.tvLM80 <- tvLM(gravedad80 ~ 0+positivo80, bw = 0.25)
summary(model.tvLM80)

coef80 <- model.tvLM80$coeff


#### Plot Figure 5

plot(model.tvLM$coeff~fecha3,type="l", ylim=c(0,21), ylab="Coefficient", xlab="", xaxt = "n")
lines(coef20~fecha3,col="blue")
lines(coef50~fecha3,col="green")
lines(coef70~fecha3,col="yellow")
lines(coef80~fecha3,col="brown")
axis(1, f3, format(f3, "%Y-%m-%d"), cex.axis = .7, las = 2)
mtext(side=1, text="Date", line=4)
abline(v=fecha3[300], col = "brown",lwd=2, lty=2) #
abline(v=fecha3[353], col = "yellow",lwd=2, lty=2) #
abline(v=fecha3[390], col = "green",lwd=2, lty=2) #
abline(v=fecha3[427], col = "blue",lwd=2, lty=2) #
legend("topright", legend=c("All","20-49","50-69","70-79","80+"),col=c("black","blue","green","yellow","brown"),lty=1, cex=0.8, inset = c(0.03, 0.075))


# decrease in the value of the coefficients
# 80+
100*(1-(coef80[520]/coef80[300]))
# ratio 70-79
100*(1-(coef70[520]/coef70[353]))
# ratio 50-69
100*(1-(coef50[520]/coef50[390]))
# ratio 20-49
100*(1-(coef20[540]/coef20[427]))
