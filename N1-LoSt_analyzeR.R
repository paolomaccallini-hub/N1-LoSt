# file name: N1-LoSt_analyze
#
# Rome, 13th October 2026
#
# It performs all the analyses on the daily variables using "all_daily_data.csv"
#
#---- Directories and subroutines ------------------------------------------------
#
current_dir<-getwd()
source("N1-LoSt_func.R")
#
for (sub in c("Histograms","Tables","Time_series","Score_analysis","Correlations",
              "Regressions","LeanTest")) {
  folder_path<-file.path(current_dir,"Output",sub)
  if (!dir.exists(folder_path)) dir.create(folder_path,recursive=TRUE)
}
#
#---- Build a table with some metadata and plot a few histograms ---------------
#
# all data
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
#
# Plot a few histograms 
#
folder_path<-file.path("Output","Histograms") 
Plot_Hist(mydata,folder_path) 
#
# Generate metadata on some variables
#
folder_path<-file.path("Output","Tables") 
results<-Table_results(mydata,folder_path)
#
#---- Study of Score (main exit variable) --------------------------------------
#
##--- Plot normalized time series-----------------------------------------------
#
folder_path<-file.path("Output","Time_series") 
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
mydata<-mydata[1:(nrow(mydata)-1),] # remove last day
#
variables<-c("date","Score_T","notafk","Score_m","Steps.C") # (date always first!)
start<-"2021-03-01"
end<-0
type<-"scores"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
##--- Autocorrelation ----------------------------------------------------------
#
folder_path<-file.path(current_dir,"Output","Score_analysis") 
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
mydata<-subset.data.frame(mydata,date>=as.integer(epoch("2021-03-01","YMD")))
mydata<-mydata[1:(nrow(mydata)-1),] # remove last day
#
# Impute missing values for Score_m
#
mydata$Score_m<-zoo::na.approx(mydata$Score_m,maxgap=10,na.rm=FALSE)
#
# Autocorrelation
#
Score_m.ts<-ts(mydata$Score_m,frequency=1) # transform outcome in a ts object
for (LM in c(100,200,300,400)) {
  file.name<-paste0("autocorrelation_Score_m_lag",LM,".jpeg")  
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=500)
  par(mfrow=c(1,2))
  acf(Score_m.ts,lwd=2,lag.max=LM) # autocorrelation
  pacf(Score_m.ts,lwd=2,lag.max=round(LM/2)) # autocorrelation corrected for previous lags (partial acf)
  dev.off()
}
#
##--- Fast Fourier Transform and Decomposition ---------------------------------
#
folder_path<-file.path(current_dir,"Output","Score_analysis") 
if(!dir.exists(folder_path)) {
  dir.create(folder_path) # create the folder for histograms if it is not present
}
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
mydata<-subset.data.frame(mydata,date>=as.integer(epoch("2021-03-01","YMD")))
mydata<-mydata[1:(nrow(mydata)-1),] # remove last day
#
# Impute missing values for Score_m
#
x<-zoo::na.approx(mydata$Score_m,maxgap=10,na.rm=FALSE)
#
# Remove trend
#
t<-seq(1,(length(x)),by=1)
window<-3*365 # trend window in days
x<-residuals(loess(x~t,span=window/length(x))) # remove trend over three years 
#
# FFT
#
X<-fft(x)
m<-length(x)
fs<-365 # frequency: measures/year
fr<-c(0:(m-1))*fs/m # frequency range
power<-(abs(X)^2)/m # power
#
xr<-c() # reconstructed signal
vector<-power[1:round(length(power)/2)]
positions<-order(-vector)[1:10] # take indexes with first 10 more powerful harmonics 
for (i in 1:length(t)) { 
  xr[i]<-0
  for (k in positions) {
    kf<-k-1                                          # 0-based frequency index
    fac<-if (kf==0 || (m%%2==0 && kf==m/2)) 1 else 2 # DC/Nyquist not doubled
    xr[i]<-xr[i]+(fac/m)*abs(X[k])*cos(2*pi*(t[i]-1)*kf/m+Arg(X[k]))
  }
}
#
# Plot both the reconstructed signal form main harmonics and periodogram
#
file.name<-"fourier.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
par(mfrow=c(2,1))
y.lim<-c(min(x),max(x))
plot(t[1:length(x)],x,type="l",lwd=1,
     xlim=c(1,length(t)),ylim=y.lim,xlab="time (days)",ylab="signal",col="black")
par(new=T)
plot(t,xr,type="l",lwd=3,
     xlim=c(1,length(t)),ylim=y.lim,xlab="",ylab="",col="red")
max.x<-max(positions)
max.y<-max(vector[-1])
max.y<-max.y+max.y/20
plot(fr,power,type="l",lwd=2,xlab="frequency (cycle/year)",xlim=c(fr[2],fr[max.x]),
     ylim=c(0,max.y)) 
abline(v=c(1:fr[max.x]),col="black",lty=2,lwd=1)
top<-positions[positions>1][1] # strongest harmonic, excluding frequency 0
legend("topright",bty="n",legend=paste0("top frequency: ",round(fr[top],3)," cycles/year",
                                        " (period ",round(fs/fr[top])," days)"))
dev.off()  
#
##--- Decomposition of score (you need to run the FTT code first!) -------------
#
# Impute missing values for Score_m
#
x<-zoo::na.approx(mydata$Score_m,maxgap=10,na.rm=FALSE)
#
# Generate time series: assume as period the one of the strongest harmonic 
#
period<-round(fs/fr[top])
score.ts<-ts(data=x,start=c(1,1),frequency=period)
str<-paste("Periods of",period,"days from",HRD(mydata$date[1]))
#
# Additive decomposition
#
x<-stats::decompose(score.ts,type="additive")
file.name<-"Additive decomposition.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
plot(x,xlab=str)
dev.off()
#
# Multiplicative decomposition
#
x<-stats::decompose(score.ts,type="multiplicative")
file.name<-"Multiplicative decomposition.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
plot(x,xlab=str)
dev.off()
#
# Seasonal and Trend decomposition using Loess (LOcally Estimated Scatterplot Smoothing)
#
x<-stl(score.ts,s.window="periodic")
file.name<-"STL_decomposition.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
plot(x)
title(str)
dev.off()
#
# Holt-Winters decomposition using exponential smoothing and forecast
#
x<-HoltWinters(score.ts)
file.name<-"holt_winters.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=1000,height=500)
plot(x,predict(x,period),xlab=str)
dev.off()
#
# Holt-Winters forecast validation
#
h<-period # days to forecast
n<-length(score.ts)
train<-window(score.ts,end=time(score.ts)[n-h])
test<-as.numeric(window(score.ts,start=time(score.ts)[n-h+1]))
hw<-HoltWinters(train)
pred.hw<-as.numeric(predict(hw,h))
model<-lm(test~pred.hw)
#
file.name<-"holt_winters_validation.jpeg"
file.name<-file.path(folder_path,file.name)
jpeg(file.name,quality=100,res=100,width=800,height=800)
plot(test,pred.hw,xlab="predicted values",ylab="values",pch=19)
abline(model,lwd=3,col="red")
r2adj<-summary(model)$adj.r.squared
TS<-summary(model)$fstatistic[1]
df1<-summary(model)$fstatistic[2]
df2<-summary(model)$fstatistic[3]
pval<-pf(TS,df1=df1,df2=df2,lower.tail=F)
str_r2adj<-sprintf(paste0("%.",5,"f"),r2adj)
str_pval<-sprintf("%.2e",pval)
title(main=paste("R2_adj=",str_r2adj,",p val=",str_pval))
dev.off()
#
##--- Search for differences in monthly scores ---------------------------------
#
folder_path<-file.path(current_dir,"Output","Score_analysis") 
if(!dir.exists(folder_path)) {
  dir.create(folder_path) # create the folder for histograms if it is not present
}
#
for (score in c("Score_m","Score_T","Steps","notafk")) {
  #
  file.name<-"Output/all_daily_data.csv"
  mydata<-read.csv(file.name,header=T,sep=",")
  mydata<-mydata[1:(nrow(mydata)-10),] # remove last days
  #
  if (score=="Score_m"|score=="Score_T") {
    mydata<-subset.data.frame(mydata,date>=as.integer(epoch("2021-03-01","YMD")))   
  } else if (score=="Steps"|score=="Steps.C") {
    mydata<-subset.data.frame(mydata,date>=as.integer(epoch("2023-01-01","YMD")))
  } else if (score=="notafk") {
    mydata<-subset.data.frame(mydata,date>=as.integer(epoch("2025-03-03","YMD")))
  }
  x<-zoo::na.approx(mydata[[score]],maxgap=15,na.rm=FALSE)
  t<-seq_along(x)
  x<-residuals(lm(x~t))
  month<-as.numeric(format(as.Date(as.POSIXct(mydata$date,origin="1970-01-01",tz="UTC")),"%m"))
  #
  file.name<-paste0(score,"_best_month.jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1200,height=800)
  boxplot(x~month,xlab="month",ylab=paste("detrended",score))
  model<-aov(x~month) # are all the means equal?
  pval<-summary(model)[[1]][["Pr(>F)"]][1]
  str_pval<-sprintf("%.2e",pval)
  title(main=paste("one-way ANOVA p-val =",str_pval))
  dev.off()
}
#
##--- Good days bye year -------------------------------------------------------
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
#
# Percentage of good days from 2021 by year
# 
mydata$Year<-as.numeric(mydata$Year)
mydata<-subset.data.frame(mydata,Year>2020)
#
# Build data frame for output
#
years<-seq(2021,max(mydata$Year,na.rm=T))
myresults<-data.frame(year=years)
#
r<-1
for (year in years) {
  df<-subset.data.frame(mydata,Year==year)
  den<-sum(!is.na(df$Score_m))
  #
  # good days
  #
  num<-0
  for (i in 1:nrow(df)) {
    if(!is.na(df$Score_1[i])&df$Score_1[i]==6) num<-num+1
    if(!is.na(df$Score_2[i])&df$Score_2[i]==6) num<-num+1
    if(!is.na(df$Score_3[i])&df$Score_3[i]==6) num<-num+1
  }
  num<-num/3
  myresults$good.days[r]<-round(num/den,2)
  #
  # middle days
  #
  num<-0
  for (i in 1:nrow(df)) {
    if(!is.na(df$Score_1[i])&df$Score_1[i]==5) num<-num+1
    if(!is.na(df$Score_2[i])&df$Score_2[i]==5) num<-num+1
    if(!is.na(df$Score_3[i])&df$Score_3[i]==5) num<-num+1
  }
  num<-num/3
  myresults$middle.days[r]<-round(num/den,2)
  #
  # bad days
  #
  num<-0
  for (i in 1:nrow(df)) {
    if(!is.na(df$Score_1[i])&df$Score_1[i]==4) num<-num+1
    if(!is.na(df$Score_2[i])&df$Score_2[i]==4) num<-num+1
    if(!is.na(df$Score_3[i])&df$Score_3[i]==4) num<-num+1
  }
  num<-num/3
  myresults$bad.days[r]<-round(num/den,2)
  #
  # mean score and standard error
  #
  myresults$Score.mean[r]<-round(mean(df$Score_m,na.rm=T),2)
  myresults$Score.sd[r]<-round(sd(df$Score_m,na.rm=T),2)
  #
  # count
  #
  num<-0
  for (i in 1:nrow(df)) {
    if(!is.na(df$Score_m[i])) num<-num+1
  }
  myresults$tot.days[r]<-round(num)
  #
  # Steps
  #
  myresults$Steps.mean[r]<-round(mean(df$Steps,na.rm=T))
  #
  # Work
  #
  myresults$Work.mean[r]<-round(mean(df$Work,na.rm=T))
  myresults$Work.a.mean[r]<-round(mean(df$Work.a,na.rm=T))
  myresults$Work.an.mean[r]<-round(mean(df$Work.an,na.rm=T))
  #
  # TA.out
  #
  myresults$TA.out.mean[r]<-round(mean(df$TA.out,na.rm=T))
  #
  # TA.out.sd
  #
  myresults$TA.out.sd.mean[r]<-round(mean(df$TA.out.sd,na.rm=T),2)
  #
  # TA.in
  #
  myresults$TA.in.mean[r]<-round(mean(df$TA.in,na.rm=T))
  #
  # TAmean
  #
  myresults$TAmean.mean[r]<-round(mean(df$TAmean,na.rm=T))
  #
  # AH.out
  #
  myresults$AH.out.mean[r]<-round(mean(df$AH.out,na.rm=T),2)
  #
  # AQI
  #
  myresults$AQI.mean[r]<-round(mean(df$AQI,na.rm=T),2)
  #
  # Provigil
  #
  myresults$Provigil.mean[r]<-round(mean(df$Provigil,na.rm=T))
  #
  # Bupropion
  #
  myresults$Bupropion.mean[r]<-round(mean(df$Bupropion,na.rm=T))
  #
  # VitE
  #
  myresults$VitE.mean[r]<-round(mean(df$VitE,na.rm=T))
  #
  # VitB1
  #
  myresults$VitB1.mean[r]<-round(mean(df$Thiamin,na.rm=T),2)
  #
  # VitD3
  #
  myresults$VitD3.mean[r]<-round(mean(df$VitD3,na.rm=T),2)
  #
  r<-r+1
}
#
# Save
#
folder_path<-file.path(current_dir,"Output","Tables") # where to put the histograms
if(!dir.exists(folder_path)) {
  dir.create(folder_path) # create the folder for histograms if it is not present
} 
write.csv(myresults,file.path(folder_path,"good_days.csv"),row.names=F,quote=F)
#
#---- Study of blood work-------------------------------------------------------
#
# Read all data
#
file.name<-"Output/all_daily_data.csv"
myexit0<-read.csv(file.name,header=T,sep=",")
#
# Prepare exit.v as a mean around blood draw
#
nd<-5 # use odd number to have blood drown at center of interval
exit.v<-"Score_m" # select exit variable (Score_m, Score_T, Sys, etc.)
myexit0<-subset.data.frame(myexit0,select=c(exit.v,"date","blood.draw","TAmean","Humidity"))
myexit<-myexit0
for (i in 1:nrow(myexit0)) {
  if (myexit0$blood.draw[i]==1) {
    myexit[[exit.v]][i]<-mean(myexit0[[exit.v]][(i-round(nd/2)):(i+round(nd/2))],na.rm=T)
  }
}
myexit<-subset.data.frame(myexit,blood.draw==1)
myexit$blood.draw<-NULL # we don't need it anymore
#
# Read blood analyses
#
myBW<-Blood.work() # Retrieve all blood work
#
# Prepare data for statistics and save it
#
mydata<-merge(myBW,myexit,by="date",all.x=T)
file.name<-paste0("Output/Tables/BloodWork_",nd,".csv")
write.table(mydata,file=file.name,quote=F,row.names=F,col.names=T,sep=",")
t<-mydata$date # for detrending
mydata$date<-NULL # we don't need it anymore
mydata$blood.draw<-NULL # we don't need it
for (j in 1:ncol(mydata)) {
  model<-lm(mydata[,j]~t,na.action=na.exclude)
  mydata[,j]<-residuals(model)
}
#
# Correlations
#
folder_path<-"Output/Correlations"
data.select<-paste0(exit.v,"_blood_work")
n.cex<-15/ncol(mydata)
mtd<-"pearson"
p.cap<-16
alpha<-0.05
Correlations(mydata,nd,data.select,folder_path,n.cex,alpha,mtd,p.cap)
#
# Regressions
#
v<-colnames(mydata)
index<-c()
for (j in v) {
  if (j!=exit.v) {
    df<-subset.data.frame(mydata,select=c(j,exit.v))
    df<-na.omit(df)
    if (nrow(df)>=5) index<-c(index,j)
  }
}
mydata<-subset.data.frame(mydata,select=c(index,exit.v))
#
folder_path<-"Output/Regressions"
data.select<-"blood_work"
Regressions(mydata,exit.v,nd,data.select,folder_path)
#
# Time series
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
mydata<-mydata[1:(nrow(mydata)-1),] # remove last day
variables<-c("date","Score_m","NSE","ALT","Neutrophils.per","Lymphocytes.per")
start<-"2021-03-01"
end<-0
type<-"blood work"
folder_path<-"Output/Time_series"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
#---- Study of lean test -------------------------------------------------------
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata0<-read.csv(file.name,header=T,sep=",")
#
# Plotting: HR
#
variables<-c("HR","HR.1","HR.2","HR.3","HR.4","HR.5","HR.6","HR.7","HR.8","HR.9","HR.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Heart Rate (bmp)"
deg<-2
LeanTest_plot(mydata,data.select,deg)             
#
# Plotting: Delta HR
#
variables<-c("Delta.HR.1","Delta.HR.2","Delta.HR.3","Delta.HR.4","Delta.HR.5",
             "Delta.HR.6","Delta.HR.7","Delta.HR.8","Delta.HR.9","Delta.HR.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Delta Heart Rate (bmp)"
deg<-2
LeanTest_plot(mydata,data.select,deg)             
#
# Plotting: Sys
#
variables<-c("Sys","Sys.1","Sys.2","Sys.3","Sys.4","Sys.5","Sys.6","Sys.7","Sys.8","Sys.9","Sys.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Systolic (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)  
#
# Plotting: Delta Sys
#
variables<-c("Delta.Sys.1","Delta.Sys.2","Delta.Sys.3","Delta.Sys.4","Delta.Sys.5",
             "Delta.Sys.6","Delta.Sys.7","Delta.Sys.8","Delta.Sys.9","Delta.Sys.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Delta Systolic (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)   
#
# Plotting: Dia
#
variables<-c("Dia","Dia.1","Dia.2","Dia.3","Dia.4","Dia.5","Dia.6","Dia.7","Dia.8","Dia.9","Dia.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Diastolic (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)  
#
# Plotting: Delta Dia
#
variables<-c("Delta.Dia.1","Delta.Dia.2","Delta.Dia.3","Delta.Dia.4","Delta.Dia.5",
             "Delta.Dia.6","Delta.Dia.7","Delta.Dia.8","Delta.Dia.9","Delta.Dia.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Delta Diastolic (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)   
#
# Plotting: MAP
#
variables<-c("MAP","MAP.1","MAP.2","MAP.3","MAP.4","MAP.5","MAP.6","MAP.7","MAP.8","MAP.9","MAP.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Mean Arterial Pressure (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
# Plotting: SVI
#
variables<-c("SVI","SVI.1","SVI.2","SVI.3","SVI.4","SVI.5","SVI.6","SVI.7","SVI.8","SVI.9","SVI.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Stroke Volume Index (ml on m2)"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
# Plotting: CPI
#
variables<-c("CPI","CPI.1","CPI.2","CPI.3","CPI.4","CPI.5","CPI.6","CPI.7","CPI.8","CPI.9","CPI.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Cardiac Power Index"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
# Plotting: PP
#
variables<-c("PP","PP.1","PP.2","PP.3","PP.4","PP.5","PP.6","PP.7","PP.8","PP.9","PP.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Sys - Dia (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg) 
#
# Time series
#
variables<-c("date","Score_m","Delta.HR.6","Delta.Sys.6","Steps.C")
start<-"2025-12-01"
end<-0
type<-"Lean Test"
folder_path<-"Output/Time_series"
nd<-7
Data.plot(mydata0,nd,variables,start,end,folder_path,type)
#
#---- Correlations, Regressions, Radar, PCA ------------------------------------
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata0<-read.csv(file.name,header=T,sep=",")
#
# Select covariates
#
data.select<-"Climate_Bresser"
exit.v<-"Score_m"
variables<-c(exit.v,c("Temp","TA.out","TA.in","TAmean","AH.out","RH.out","Pressure","Light",
                      "Clouds","Precipitation","PM10","PM2.5","SO2","NO","AQI","CO"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-5
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Climate_BG"
exit.v<-"Score_m"
variables<-c(exit.v,c("Temp","MRT","TG","TA","AH","RH","TAmean","Pressure","Light",
                      "Clouds","Precipitation","PM10","PM2.5","SO2","NO","AQI","CO"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-5
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain_L"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF",
                      "Delta_L","Theta_L","Alpha_L","Beta_L","Gamma_L",
                      "L.Outer.HbDiff","L.Outer.HbTot"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-1
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain_R"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF",
                      "Delta_R","Theta_R","Alpha_R","Beta_R","Gamma_R",
                      "R.Outer.HbDiff","R.Outer.HbTot"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-1
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF",
                      "Delta","Theta","Alpha","Beta","Gamma",
                      "Outer.HbDiff","Outer.HbTot"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-1
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Lean Test"
exit.v<-"Score_m"
variables<-c(exit.v,c("HR","Delta.HR.1","Delta.HR.2","Delta.HR.3","Delta.HR.4","Delta.HR.5",
                      "Delta.HR.6","Delta.HR.7","Delta.HR.8","Delta.HR.9","Delta.HR.10",
                      "Sys","Sys.1","Sys.2","Sys.3","Sys.4","Sys.5","Sys.6","Sys.7","Sys.8","Sys.9","Sys.10",
                      "Dia","Dia.1","Dia.2","Dia.3","Dia.4","Dia.5","Dia.6","Dia.7","Dia.8","Dia.9","Dia.10"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-1
CoReRaPCA(mydata,exit.v,data.select,nd)
#
#---- Marginal dynamic regressions ---------------------------------------------
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata0<-read.csv(file.name,header=T,sep=",")
#
# Select exit variable
#
exit.vs<-c("Score_m","Score_T","notafk") 
#
# Select covariates 
#
covariates<-c("Provigil","Bupropion","Sulbutiamine",
              "Aniracetam","D_serine","Sarcosine",
              "Beta.Lactam","Corticosteroid","Nasonex",
              "Pregnenolone","DHEA","Finasteride.Topic","Melatonin","T4","T3",
              "Destrometorfano","Tachipirine","Ibuprofene",
              "Compression.suit","Compression.socks",
              "Ketotifen","Desloratadina",
              "Pioglitazone","Deursil","Bismuth.subsalicylate",
              "Mestinon","Etilefrine","Acetazolamide","Nimodipine",
              "Creatine","ALCar","B.complex",
              "VitE","VitD3","HydroxycobalaminTot","B12.Tot",
              "Thiamin","Riboflavin","Niacin","VitB6","Folic.acid",          
              "Folinic.acid","Cyanocobalamin","Biotin",               
              "Pantothenic.acid","Choline.bitartrate","Myo.Inositol",  
              "Omega3","EPA","DHA",
              "Work.an","Work.a","Steps","Steps.B","Steps.C","Bike","Run","ST",
              "Assistant","Cleaning","Relatives","FLU","Travel","blood.draw",
              "kcal","protein.g","carbo.g","fat.g",
              "TA.in","TA.in.out","AH.in","RH.in",
              "TA.out","AH.out","RH.out","Light",
              "TA","TG","MRT","AH","RH","TG.night","TG.day.night","TG.day",
              "PM2.5","PM10","NO","NO2","CO","SO2","AQI",
              "TAmean","Pressure","Humidity","Clouds","Abs.Humidity",
              "Lamp.power","Lamp.peak_lambda","Lamp.min.on",
              "SpO2","Temp","glucose","ketones",
              "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF","TotPower",
              "Delta","Theta","Alpha","Beta","Gamma",
              "Delta_L","Theta_L","Alpha_L","Beta_L","Gamma_L",
              "Delta_R","Theta_R","Alpha_R","Beta_R","Gamma_R",
              "L.Outer.HbDiff","L.Outer.HbTot",
              "R.Outer.HbDiff","R.Outer.HbTot",
              "L.Outer.Hb","L.Outer.Hb02",
              "R.Outer.Hb","R.Outer.Hb02",
              "L.Outer.660","R.Outer.660",
              "L.Inner.HbDiff","L.Inner.HbTot",
              "R.Inner.HbDiff","R.Inner.HbTot",
              "L.Inner.Hb","L.Inner.Hb02",
              "R.Inner.Hb","R.Inner.Hb02",
              "L.Inner.660","R.Inner.660",
              "HR","Delta.HR.1","Delta.HR.2","Delta.HR.3","Delta.HR.4","Delta.HR.5",
              "Delta.HR.6","Delta.HR.7","Delta.HR.8","Delta.HR.9","Delta.HR.10",
              "Sys","Sys.1","Sys.2","Sys.3","Sys.4","Sys.5","Sys.6","Sys.7","Sys.8","Sys.9","Sys.10",
              "Dia","Dia.1","Dia.2","Dia.3","Dia.4","Dia.5","Dia.6","Dia.7","Dia.8","Dia.9","Dia.10",
              "SVI","SVI.1","SVI.2","SVI.3","SVI.4","SVI.5","SVI.6","SVI.7","SVI.8","SVI.9","SVI.10",
              "CPI","CPI.1","CPI.2","CPI.3","CPI.4","CPI.5","CPI.6","CPI.7","CPI.8","CPI.9","CPI.10",
              "MAP","MAP.1","MAP.2","MAP.3","MAP.4","MAP.5","MAP.6","MAP.7","MAP.8","MAP.9","MAP.10",
              "PP","PP.1","PP.2","PP.3","PP.4","PP.5","PP.6","PP.7","PP.8","PP.9","PP.10")
#
# Marginal regression
#
type<-"marginal"
alpha<-0.05
DynaReg(mydata0,exit.vs,covariates,type,alpha)
#
#---- Multinomial dynamic regressions ------------------------------------------
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata0<-read.csv(file.name,header=T,sep=",")
#
# Select exit variable
#
exit.vs<-c("Score_m","Score_T","max.Delta.HR") 
#
# Select covariates 
#
data.select<-"Climate_BG"
covariates<-c("TG","AH","TAmean","Pressure","Light","Clouds","Precipitation","AQI")
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select)
#
# Select covariates 
#
data.select<-"Brain Left"
covariates<-c("Delta_L","Theta_L","Alpha_L","Beta_L",
              "L.Outer.Hb02","L.Outer.HbTot",
              "L.Inner.Hb02","L.Inner.HbTot",
              "MRT","Temp","SpO2","Steps.B")
days<-c(0,1,2,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
# Select covariates 
#
data.select<-"Brain Right"
covariates<-c("Delta_R","Theta_R","Alpha_R","Beta_R",
              "R.Outer.Hb02","R.Outer.HbTot",
              "R.Inner.Hb02","R.Inner.HbTot",
              "MRT","Temp","SpO2","Steps.B")
days<-c(0,1,2,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
# Select covariates 
#
data.select<-"Brain"
covariates<-c("Theta","Alpha","Beta","Gamma",
              "Outer.Hb02","Outer.HbTot",
              "Inner.Hb02","Inner.HbTot",
              "MRT","Temp","SpO2","Steps.B")
days<-c(0,1,2,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
# Select covariates 
#
data.select<-"HRV"
covariates<-c("SDNN","pNN50","rMSSD","LFonHF","Steps.C","TG","Provigil","Mestinon","Etilefrine")
days<-c(0,1,2,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
# Select covariates 
#
data.select<-"Drugs"
covariates<-c("Provigil","Bupropion",
              "Aniracetam","D_serine","Sarcosine","Destrometorfano",
              "Beta.Lactam",
              "Tachipirine","Ibuprofene",
              "Corticosteroid","Pregnenolone","DHEA","Melatonin","T4","T3",
              "Ketotifen","Desloratadina",
              "Pioglitazone","Deursil",
              "Mestinon","Etilefrine","Acetazolamide",
              "Nimodipine","Compression.suit",
              "Creatine","ALCar",
              "VitE","VitD3","B12.Tot",
              "Thiamin","Riboflavin","Niacin","VitB6","Folic.acid",          
              "Folinic.acid","Biotin",   
              "Pantothenic.acid","Choline.bitartrate","Myo.Inositol",  
              "DHA","EPA")
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select)
#
# Select covariates 
#
data.select<-"Exercise"
covariates<-c("Steps.C","WL.bin","Bike","TA.out","ST")
days<-c(0,1,2,4,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
# Select covariates 
#
data.select<-"Lean test"
covariates<-c("max.Delta.HR","min.PP","min.Delta.Sys","ST","ST.hour","Steps.C","TA.out")
days<-c(0,1,2,4,7)
DynaReg(mydata0,exit.vs,covariates,type="multinomial",alpha=0.05,data.select,days)
#
#---- Normalized diagrams of time series ---------------------------------------
#
folder_path<-file.path(current_dir,"Output","Time_series") 
#
# Read
#
file.name<-"Output/all_daily_data.csv"
mydata<-read.csv(file.name,header=T,sep=",")
mydata<-mydata[1:(nrow(mydata)-1),] # remove last day
#
# normalized diagrams (date always first!)
#
variables<-c("date","Score_m","L.Outer.HbTot","L.Outer.HbDiff","L.Inner.HbTot","L.Inner.HbDiff")
start<-"2025-09-19"
end<-0
type<-"Left_fNIRS"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type,detrend=T)
#
variables<-c("date","Score_m","R.Outer.HbTot","R.Outer.HbDiff","R.Inner.HbTot","R.Inner.HbDiff")
start<-"2025-09-19"
end<-0
type<-"Right_fNIRS"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type,detrend=T)
#
variables<-c("date","Score_m","Outer.HbTot","Outer.HbDiff","Inner.HbTot","Inner.HbDiff")
start<-"2025-09-19"
end<-0
type<-"fNIRS"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type,detrend=T)
#
variables<-c("date","Score_m","LF","LFnu","HF","LFonHF")
start<-"2022-07-09"
end<-0
type<-"HRV_freq"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","Sys","Dia","HR")
start<-"2023-06-29"
end<-0
type<-"cardiac"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","Steps","Steps.B","Steps.C")
start<-"2022-12-01"
end<-0
type<-"steps"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","Work","Work.a","Work.an")
start<-"2023-01-01"
end<-0
type<-"work"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","TA.out","RH.out","AH.out")
start<-"2021-12-06"
end<-0
type<-"TA"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","Weight","Skeletal.muscle","Body.Fat")
start<-"2021-05-03"
end<-0
type<-"BIA"
nd<-30
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","glucose","ketones")
start<-"2024-02-06"
end<-"2024-10-06"
type<-"metabolism"
nd<-7
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","TG","MRT","TA")
start<-"2024-01-29"
end<-0
type<-"TG"
nd<-15
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","TP9_delta","TP9_theta","TP9_alpha","TP9_beta","TP9_gamma")
start<-"2025-06-17"
end<-0
type<-"TP9"
nd<-10
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","TP10_delta","TP10_theta","TP10_alpha","TP10_beta","TP10_gamma")
start<-"2025-06-17"
end<-0
type<-"TP10"
nd<-10
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","AF7_delta","AF7_theta","AF7_alpha","AF7_beta","AF7_gamma")
start<-"2025-06-17"
end<-0
type<-"AF7"
nd<-10
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","AF8_delta","AF8_theta","AF8_alpha","AF8_beta","AF8_gamma")
start<-"2025-06-17"
end<-0
type<-"AF8"
nd<-10
Data.plot(mydata,nd,variables,start,end,folder_path,type)
#
variables<-c("date","Score_m","Delta","Theta","Alpha","Beta","Gamma")
start<-"2025-06-17"
end<-0
type<-"EEG"
nd<-10
Data.plot(mydata,nd,variables,start,end,folder_path,type)
