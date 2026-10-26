# file name: N1-LoSt_func
#
# year: 2024
# month: 06
# day: 19
#
#---- Packages, parameters, API keys -------------------------------------------
#
library(RHRV) # for HRv data analysis
library(leaps) # for model reduction
library(mixR) # for bimodal distributons
library(cluster) # needed for pam
library(corrplot) # for corrplot
library(runner) # for moving mean
library(scatterplot3d) # needed for 3D visualization
library(MASS) # for discriminant analysis
library(e1071) # for naiveBayes()
library(boot) # for cv.glm()
library(car) # for vif()
library(xts)
library(zoo)
library(pracma) # for numerical integration
library(ellipse) # for generating ellipse
library(readxl) # for reading .xlsx format
library(fmsb) # for the radar chart
library(ppcor) # partial correlation
library(readr)
library(httr) # for API
library(jsonlite) # for API
library(dplyr) 
library(purrr)
library(lubridate)
library(rgl) # for 3D plotting 
library(data.table) # for fread
library(eegkit) # for EEG Fourier analysis
library(pracma) # for numeric integration
library(tidyr)
library(ggplot2)
library(scales)
#
# parameters
#
R<-8.31 # constant for ideal gasses (J/(K*mole))
mu_v<-18.01528 # molecular weight of water (g/mole)
mu_d<-28.97 # molecular weight of dry air (g/mole)
LR<-16.7 # Lewis ratio (K/kPa)
set.seed(12) # seed
#
# API keys
#
OWM_key<-c()
OWM_key<-c(Sys.getenv("OWM_key1"),Sys.getenv("OWM_key2"))
if (any(OWM_key=="")) warning("OpenWeather API keys not found in .Renviron")
#
#---- TimeStamp handling ------------------------------------------------------
#
# Convert human readable date (HRD) into Unix epoch (or Unix time or POSIX time or 
# Unix timestamp). It works on multiple HRD formats. 
#
epoch<-function(date,format) {
  date<-as.character(date)
  #
  # Define possible date formats
  #
  if (format=="MDY") {
    formats <- c("%m/%d/%Y %H:%M:%S","%m/%d/%Y","%m-%d-%Y %H:%M:%S","%m-%d-%Y")
  } else if (format=="DMY") {
    formats <- c("%d/%m/%Y %H:%M:%S","%d/%m/%Y","%d-%m-%Y %H:%M:%S","%d-%m-%Y")
  } else if (format=="YDM") {
    formats <- c("%Y/%d/%m %H:%M:%S","%Y/%d/%m","%Y-%d-%m %H:%M:%S","%Y-%d-%m")
  } else if (format=="YMD") {  
    formats <- c("%Y/%m/%d %H:%M:%S","%Y/%m/%d","%Y-%m-%d %H:%M:%S","%Y-%m-%d")
  } else {
    stop("date format not supported")
  }
  #
  # Try each format until one works
  #
  for (fmt in formats) {
    EPC <- try(as.POSIXct(date,format=fmt,tz = "UTC"),silent=TRUE)
    if (!is.na(EPC)) return(as.numeric(EPC))
  }
}
#
# Convert Unix epoch (or Unix time or POSIX time or Unix timestamp) into HRD 
# (human readable date, YYYY-MM-DD)
#
HRD<-function(epoch) {
  date<-as.character(as.POSIXct(as.numeric(epoch),origin="1970-01-01", tz="GMT"))
  return(date)
}
#
#---- Equivalent drugs --------------------------------------------------------
#
# Add a fictitious antibiotic
#
Beta.Lactam<-function(mydata) {
  #
  # 1 g Amoxicillin = 1 Beta.Lactam
  # 833 UI Ospen = 1 Beta.Lactam
  # 500 mg Penilevel = 1 Beta.Lactam
  # 500 mg Cefuroxime = 1 Beta.Lactam
  #
  Beta.Lactam<-rep(NA,nrow(mydata))
  for (i in 1:nrow(mydata)) {
    if (!is.na(mydata$Amoxicillin[i])) Beta.Lactam[i]<-mydata$Amoxicillin[i]*1
    if (!is.na(mydata$Ospen[i])) {
      if (is.na(Beta.Lactam[i])) {
        Beta.Lactam[i]<-mydata$Ospen[i]*(1/833)
      } else {
        Beta.Lactam[i]<-Beta.Lactam[i]+mydata$Ospen[i]*(1/833)
      }
    }
    if (!is.na(mydata$Penilevel[i])) {
      if (is.na(Beta.Lactam[i])) {
        Beta.Lactam[i]<-mydata$Penilevel[i]*(1/500)
      } else {
        Beta.Lactam[i]<-Beta.Lactam[i]+mydata$Penilevel[i]*(1/500)
      }
    }
    if (!is.na(mydata$Cefuroxime[i])) {
      if (is.na(Beta.Lactam[i])) {
        Beta.Lactam[i]<-mydata$Cefuroxime[i]*(1/500)
      } else {
        Beta.Lactam[i]<-Beta.Lactam[i]+mydata$Cefuroxime[i]*(1/500)
      }
    }
  }
  mydata$Beta.Lactam<-Beta.Lactam
  return(mydata)
}
#
# Add a fictitious corticosteroid (prednisone-equivalent dose, mg)
#
Corticosteroid<-function(mydata) {
  #
  # 1 mg Prednisone = 1 mg Prednisone.eq
  # 6 mg Deflazacort = 5 mg Prednisone.eq
  #
  # Not included:
  # Fludrocortisone: mineralocorticoid, negligible glucocorticoid effect at 0.025-0.1 mg
  # Nasonex: nasal mometasone, negligible systemic absorption
  #
  P<-mydata$Prednisone
  D<-mydata$Deflazacort*5/6
  mydata$Corticosteroid<-ifelse(is.na(P)&is.na(D),NA,rowSums(cbind(P,D),na.rm=TRUE))
  return(mydata)
}
#
# Total B12 as intramuscular-equivalent dose (mcg)
#
B12.Tot<-function(mydata) {
  #
  # 1 mcg IM hydroxycobalamin = 1 mcg B12.Tot
  # 1 mcg PO hydroxycobalamin = 1/20 mcg B12.Tot
  # 1 mcg PO cyanocobalamin   = 1/20 mcg B12.Tot
  #
  B12<-cbind(mydata$HydroxycobalaminIM,
             mydata$HydroxycobalaminPO/20,
             mydata$Cyanocobalamin/20)
  mydata$B12.Tot<-ifelse(rowSums(!is.na(B12))==0,NA,rowSums(B12,na.rm=TRUE))
  return(mydata)
}
#
# Convert PO hydroxycobalamin to IM hydroxycobalamin
#
Hydroxycobalamin.Tot<-function(mydata) {
  #
  # 1 mg PO hydroxycobalamin = 1/20 mg IM hydroxycobalamin
  #
  HydroxycobalaminTot<-rep(NA,nrow(mydata))
  for (i in 1:nrow(mydata)) {
    if (!is.na(mydata$HydroxycobalaminPO[i])) {
      if (is.na(mydata$HydroxycobalaminIM[i])) {
        HydroxycobalaminTot[i]<-mydata$HydroxycobalaminPO[i]*1/20
      } else {
        HydroxycobalaminTot[i]<-mydata$HydroxycobalaminIM[i]+mydata$HydroxycobalaminPO[i]*1/20
      }
    }
  }
  mydata$HydroxycobalaminTot<-HydroxycobalaminTot
  return(mydata)
}
#
#---- Cardiac variables -------------------------------------------------------
#
# It calculates some derived cardiac parameters 
# 
# Linear models for SVI comes from MC Van Campen 2023. It refers to 52 HCs in
# supine position and end tilt. I also tried the model for 362 ME patients in supine position
# and it gives worse regressions with my symptoms, so I removed it from the subroutine.
#
Cardiac.Derived<-function(mydata) {
  #
  K<-2.21*1e-06 # conversion from mmHg*ml/min to watt
  #
  # Supine MAP, CPI, PP
  #
  mydata$MAP<-mydata$Dia+(mydata$Sys-mydata$Dia)/3 # mean arterial pressure (mmHg)
  mydata$SVI<-(99.18-mydata$HR)/0.9769 # supine stroke volume index healthy (ml/m^2)
  mydata$CPI<-mydata$SVI*mydata$HR*mydata$MAP # cardiac power index (ml*mmHg/(min*m^2))
  mydata$CPI<-K*mydata$CPI # cardiac power index (W/m^2)
  mydata$PP<-mydata$Sys-mydata$Dia # delta between max pressure - min pressure (pulse pressure) (mmHg)
  #
  # Standing MAP, CPI, PP
  #
  for (i in 1:10) {
    #
    mydata[[paste0("MAP",".",i)]]<-mydata[[paste0("Dia",".",i)]]+(mydata[[paste0("Sys",".",i)]]-mydata[[paste0("Dia",".",i)]])/3 
    mydata[[paste0("SVI",".",i)]]<-(130.7-mydata[[paste0("HR",".",i)]])/1.968
    mydata[[paste0("CPI",".",i)]]<-mydata[[paste0("SVI",".",i)]]*mydata[[paste0("HR",".",i)]]*mydata[[paste0("MAP",".",i)]]
    mydata[[paste0("CPI",".",i)]]<-K*mydata[[paste0("CPI",".",i)]] 
    mydata[[paste0("PP",".",i)]]<-mydata[[paste0("Sys",".",i)]]-mydata[[paste0("Dia",".",i)]] 
    #
  }
  return(mydata)
}
#
# Time domain and frequency domain analysis
#
HRV<-function(hrv.data,days) {
  #
  # filter the RR intervals and calculate the instantaneous non-interpolated HR
  #
  for (i in 1:days) {
    hrv.data[[i]]<-BuildNIHR(hrv.data[[i]])
    hrv.data[[i]]<-FilterNIHR(hrv.data[[i]])
  }
  #
  # interpolated HR
  #
  for (i in 1:days) {
    hrv.data[[i]]<-InterpolateNIHR(hrv.data[[i]],freqhr=4)
  }
  #
  # time domain analysis
  #
  for (i in 1:days) {
    hrv.data[[i]]<-CreateTimeAnalysis(hrv.data[[i]])
  }
  #
  # frequency domain analysis
  #
  for (i in 1:days) {
    hrv.data[[i]]<-CreateFreqAnalysis(hrv.data[[i]])
    # hrv.data[[i]]<-CalculatePSD(hrv.data[[i]],method="pgram",doPlot=F) # alternative method, needs a different downstream script
    hrv.data[[i]]<-CalculatePowerBand(hrv.data[[i]],indexFreqAnalysis=1,LFmin=0.05,
                                      LFmax=0.15,HFmin=0.15,HFmax=0.4,size=300,shift=300,
                                      type="fourier")
  }
  return(hrv.data)
}
#
# Derive deltas from Lean test
#
LeanTest_func<-function(Score){
  #
  # deltas for HR
  #
  variables<-c("HR","HR.1","HR.2","HR.3","HR.4","HR.5","HR.6","HR.7","HR.8","HR.9","HR.10","date")
  mydata<-subset.data.frame(Score,select=variables)
  for (j in 2:11) {
    mydata[,j]<-mydata[,j]-mydata[,1]
    colnames(mydata)[j]<-paste0("Delta.HR.",(j-1))
  }
  mydata.HR<-mydata[,-1] # remove first column
  #
  # deltas for Sys
  #
  variables<-c("Sys","Sys.1","Sys.2","Sys.3","Sys.4","Sys.5","Sys.6","Sys.7","Sys.8","Sys.9","Sys.10","date")
  mydata<-subset.data.frame(Score,select=variables)
  for (j in 2:11) {
    mydata[,j]<-mydata[,j]-mydata[,1]
    colnames(mydata)[j]<-paste0("Delta.Sys.",(j-1))
  }
  mydata.Sys<-mydata[,-1] # remove first column
  #
  # deltas for Dia
  #
  variables<-c("Dia","Dia.1","Dia.2","Dia.3","Dia.4","Dia.5","Dia.6","Dia.7","Dia.8","Dia.9","Dia.10","date")
  mydata<-subset.data.frame(Score,select=variables)
  for (j in 2:11) {
    mydata[,j]<-mydata[,j]-mydata[,1]
    colnames(mydata)[j]<-paste0("Delta.Dia.",(j-1))
  }
  mydata.Dia<-mydata[,-1] # remove first column
  #
  mydata<-merge(mydata.HR,mydata.Sys,by="date",all=T)
  mydata<-merge(mydata,mydata.Dia,by="date",all=T)
  #
  # Add overall measures for standing test and return results
  #
  mydata$max.Delta.HR<-do.call(pmax,c(mydata[paste0("Delta.HR.",1:10)],na.rm=TRUE))
  mydata$min.PP<-do.call(pmin,c(mydata[paste0("PP.",1:10)],na.rm=TRUE))
  mydata$min.Delta.Sys<-do.call(pmin,c(mydata[paste0("Delta.Sys.",1:10)],na.rm=TRUE))
  #
  return(mydata)
}
#
# Used to plot results of lean tests
#
LeanTest_plot<-function(mydata,data.select,deg) {
  #
  # Figure name and directory
  #
  file.name<-paste0(data.select,".jpeg")
  file.name<-file.path("Output/LeanTest",file.name)
  jpeg(file.name,quality=100,res=150,width=2000,height=1000)
  #
  # Time axis
  #
  if (ncol(mydata)==11) {
    time<-seq(0,10,1)
  } else {
    time<-seq(1,10,1)
  }
  #
  # Plot data
  #
  bp<-boxplot(mydata,names=time,xlab="Minutes standing",ylab=data.select,
              ylim=c(min(mydata,na.rm=T)-min(mydata,na.rm=T)/5,
                     max(mydata,na.rm=T)+max(mydata,na.rm=T)/5),
              main=data.select) # plot boxes
  #
  means<-sapply(mydata,mean,na.rm=T) # calculate means
  points(x=1:length(means),y=means,pch=19,col="red") # add means to plot
  means<-round(means,0)
  text(x=1:length(means),y=means+3,labels=means,cex=1,col="red")
  #
  ns<-sapply(mydata, function(x) sum(!is.na(x))) # calculate sample size
  text(x=1:length(ns),y=min(mydata,na.rm=T)-5,labels=ns,cex=1)
  #
  # Polynomial regression HR ~ time
  #
  long<-data.frame(time=rep(time,each=nrow(mydata)),HR=unlist(mydata,use.names=FALSE))
  long<-long[!is.na(long$HR),]
  fit<-lm(HR~stats::poly(time,deg),data=long) # stats:: avoids pracma::poly
  #
  new<-data.frame(time=seq(0,max(time),length.out=200))
  pr<-predict(fit,newdata=new,interval="confidence")
  #
  polygon(c(new$time,rev(new$time))+1,c(pr[,"lwr"],rev(pr[,"upr"])),
          col=adjustcolor("blue",alpha.f=0.2),border=NA)   # 95% confidence band
  lines(new$time+1,pr[,"fit"],col="blue",lwd=2)            # fitted curve
  legend("topleft",bty="n",text.col="blue",legend=paste("polynomial degree",deg,
                                                        ", R2dj =",round(summary(fit)$adj.r.squared,2)))
  #
  dev.off()
}
#
#---- Environmental variables --------------------------------------------------
#
# It calculates absolute humidity (g/m^3) and water vapour pressure
#
AHum<-function(mydata) {
  TA <- mydata$TA
  ep <- (mydata$RH / 100) * 6.112 * exp(17.62 * TA / (243.12 + TA))
  mydata$AH <- round(mu_v / (R * (273.15 + TA)) * 100 * ep, 2)
  mydata
}
#
# It calculates mean radiant temperature (°C)
#
MRT<-function(mydata, v = 0.3, epsilon = 0.95, D_g = 0.05) {
  TG <- mydata$TG
  TA <- mydata$TA
  k  <- (1.1e8 * v^0.6) / (epsilon * D_g^0.4)   # computed once
  mydata$MRT <- round(((TG + 273.15)^4 + k * (TG - TA))^0.25 - 273.15, 2)
  mydata
}
#
# Load and edit data from Tenmars black globe thermometer
# 
Tenmars<-function() {
  #
  # read all files 
  #
  folder_path<-file.path(getwd(),"Data","BlackGlobe")
  file.names<-list.files(folder_path,pattern="\\.csv$",ignore.case=TRUE,full.names=TRUE)
  array<-unlist(lapply(file.names,readLines,warn=FALSE))
  #
  # files are Latin-1: convert so that the degree sign is read correctly
  #
  array<-iconv(array,from="latin1",to="UTF-8")
  #
  # keep only lines starting with a date and time
  #
  array<-array[grepl("^[0-9]{2}/[0-9]{2}/[0-9]{4} [0-9]{2}:",array)]
  #
  # extract values: six numbers followed by " °C", one followed by " %"
  #
  to.num<-function(x) as.numeric(sub(",",".",x,fixed=TRUE))
  celsius<-regmatches(array,gregexpr("-?[0-9]+,[0-9]+(?= \u00b0C)",array,perl=TRUE))
  percent<-regmatches(array,gregexpr("[0-9]+,[0-9]+(?= %)",array,perl=TRUE))
  #
  ok<-lengths(celsius)==6 & lengths(percent)==1
  if (any(!ok)) warning(sum(!ok)," malformed lines skipped")
  array<-array[ok]
  celsius<-matrix(to.num(unlist(celsius[ok])),ncol=6,byrow=TRUE)
  #
  # build data frame
  #
  mydata<-data.frame(date=as.numeric(as.POSIXct(substr(array,1,10),format="%d/%m/%Y",tz="UTC")),
                     time=as.numeric(substr(array,12,13)),
                     WBGT.in=celsius[,1],WBGT.out=celsius[,2],TA=celsius[,3],
                     TG=celsius[,4],WET=celsius[,5],DEW=celsius[,6],
                     RH=to.num(unlist(percent[ok])))
  #
  # sort readings chronologically
  #
  ts<-as.numeric(as.POSIXct(substr(array,1,19),format="%d/%m/%Y %H:%M:%S",tz="UTC"))
  mydata<-mydata[order(ts),]
  #
  # add absolute humidity
  #
  mydata<-AHum(mydata)
  #
  # add mean radiant temperature
  #
  mydata<-MRT(mydata)
  #
  # calculate daily means and standard deviations
  #
  days<-unique(mydata$date)
  n.days<-length(days)
  mydata.daily<-data.frame(date=rep(NA,n.days),
                           WBGT.in=rep(NA,n.days),WBGT.out=rep(NA,n.days),WET=rep(NA,n.days),
                           TA=rep(NA,n.days),TA.sd=rep(NA,n.days),TG=rep(NA,n.days),TG.sd=rep(NA,n.days),
                           MRT=rep(NA,n.days),RH=rep(NA,n.days),AH=rep(NA,n.days),
                           HL=rep(NA,n.days),HS=rep(NA,n.days),CpR=rep(NA,n.days),
                           Cres=rep(NA,n.days),Eres=rep(NA,n.days),Esk=rep(NA,n.days),
                           TG.night=rep(NA,n.days),TG.day=rep(NA,n.days),TG.day.night=rep(NA,n.days))
  for (i in 1:n.days) {
    df<-subset.data.frame(mydata,date==days[i])
    mydata.daily$date[i]<-days[i]
    mydata.daily$WBGT.in[i]<-mean(df$WBGT.in,na.rm=T)
    mydata.daily$WBGT.out[i]<-mean(df$WBGT.out,na.rm=T)
    mydata.daily$TA[i]<-mean(df$TA,na.rm=T)
    mydata.daily$TA.sd[i]<-sd(df$TA,na.rm=T)
    mydata.daily$TG[i]<-mean(df$TG,na.rm=T)
    mydata.daily$TG.sd[i]<-sd(df$TG,na.rm=T)
    mydata.daily$MRT[i]<-mean(df$MRT,na.rm=T)
    mydata.daily$WET[i]<-mean(df$WET,na.rm=T)
    mydata.daily$RH[i]<-mean(df$RH,na.rm=T)
    mydata.daily$AH[i]<-mean(df$AH,na.rm=T)
    wake<-subset.data.frame(df,time>8&time<22)
    mydata.daily$TG.day[i]<-mean(wake$TG,na.rm=T)
    if (i>1) {
      yesterday<-subset.data.frame(mydata,date==days[i-1])
      night_1<-subset.data.frame(yesterday,time>=22)
      night_2<-subset.data.frame(df,time<=8)
      mydata.daily$TG.night[i]<-mean(c(night_1$TG,night_2$TG),na.rm=T)
      mydata.daily$TG.day.night[i]<-mydata.daily$TG.day[i]/mydata.daily$TG.night[i]
    }
  }
  #
  return(mydata.daily)
  #
}
#
# Load and edit data from Bresser thermometer
# 
Bresser<-function() {
  #
  # read
  #
  current_dir<-getwd()
  sub_dir<-"Bresser"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.names<-list.files(folder_path)
  #
  for (i in 1:length(file.names)) {
    if (substr(file.names[[i]],1,12)=="Weathercloud") {
      f.name<-file.path(folder_path,file.names[[i]])
      if (i==1) {
        Bresser<-fread(f.name,header=F,skip=1,fill=Inf)
        Bresser<-data.frame(lapply(Bresser,function(x) gsub(",", ".", x)))
      }
      if (i>1) {
        df<-fread(f.name,header=F,skip=1,fill=Inf)[,1:18]
        df<-data.frame(lapply(df,function(x) gsub(",", ".", x)))
        Bresser<-rbind(Bresser,df)
      }
    }
  }
  Bresser[,2:18]<-data.frame(lapply(Bresser[,2:18],function(x) as.numeric(x)))
  colnames(Bresser)<-c("Date","TA.in","TA.out","Wind.chill","Dew.in","Dew.out",
                       "Index.in","Index.out",
                       "RH.in","RH.out","Wind.max","Wind.min","Wind.dir","Pressure.B",
                       "Rain","Rain2","Light","UV")
  Bresser<-subset.data.frame(Bresser,select=c("Date","TA.in","TA.out","Dew.out",
                                              "Index.out","RH.in","RH.out","Light"))
  #
  # calculate AH.in
  #
  df<-data.frame(TA=Bresser$TA.in,RH=Bresser$RH.in)
  df<-AHum(df)
  Bresser$AH.in<-df$AH
  #
  # calculate AH.out
  #
  df<-data.frame(TA=Bresser$TA.out,RH=Bresser$RH.out)
  df<-AHum(df)
  Bresser$AH.out<-df$AH
  #
  # Convert date into epoch (DD/MM/YYYY or YYYY-MM-DD), once per distinct date
  #
  d<-substr(Bresser$Date,1,10)
  ud<-unique(d)
  ue<-vapply(ud,function(x) {
    if (substr(x,3,3)=="/") epoch(x,format="DMY") else stop("Date is in the wrong format!")
  },numeric(1),USE.NAMES=FALSE)
  Bresser$date<-ue[match(d,ud)]
  #
  # Calculate daily mean values
  #
  dates<-unique(Bresser$date)
  v<-rep(NA,length(dates))
  mydata<-data.frame(TA.in=v,RH.in=v,AH.in=v,TA.out=v,TA.in.out=v,TA.out.sd=v,RH.out=v,AH.out=v,Light=v,date=v)
  for (i in 1:length(dates)) {
    df<-subset.data.frame(Bresser,date==dates[i])
    mydata$TA.in[i]<-mean(df$TA.in,na.rm=T)
    mydata$TA.out[i]<-mean(df$TA.out,na.rm=T)
    mydata$RH.in[i]<-mean(df$RH.in,na.rm=T)
    mydata$AH.in[i]<-mean(df$AH.in,na.rm=T)
    mydata$TA.in.out[i]<-mean(c(mydata$TA.in[i],mydata$TA.out[i]),na.rm=T)
    mydata$TA.out.sd[i]<-sd(df$TA.out,na.rm=T)
    mydata$RH.out[i]<-mean(df$RH.out,na.rm=T)
    mydata$AH.out[i]<-mean(df$AH.out,na.rm=T)
    mydata$AH.out.sd[i]<-sd(df$AH.out,na.rm=T)
    mydata$Light[i]<-mean(df$Light,na.rm=T)
    mydata$date[i]=dates[i]
  }
  #
  return(mydata)
} 
#
# Derive coordinates from city name
# 
Get.coord<-function(city.name,country.name,api_key,limit) {
  #
  url<-paste0(
    "http://api.openweathermap.org/geo/1.0/direct?q=", 
    URLencode(city.name), 
    "&limit=",limit, 
    "&appid=",api_key
  )
  #
  # Make the API request
  #
  response<-GET(url)
  #
  # Check if the request was successful
  #
  if (status_code(response)!=200) {
    stop("Failed to retrieve data. Status code:",status_code(response))
  }
  #
  # Parse the JSON response
  #
  data<-fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  #
  # Check if data is empty
  #
  if (length(data)==0) {
    return(paste("No data found for city:",city_name))
  }
  #
  # Extract latitued
  #
  data<-subset.data.frame(data,country==country.name)
  data<-data[1,]
  #
  return(data)
}
#
# Retrieve daily environmental data from Open Weather Map API
# 
Get.weather<-function(timestamp,time.zone,lat,lon,place,country,api_key) {
  #
  url<-paste0(
    "https://api.openweathermap.org/data/3.0/onecall/day_summary?",
    "lat=", lat,
    "&lon=", lon,
    "&date=", HRD(timestamp),
    "&tz=", time.zone,
    "&appid=", api_key,
    "&units=", "metric"
  )
  #
  # Make the API request
  #
  response<-GET(url)
  #
  # Check if the request was successful
  #
  if (status_code(response)!=200) {
    print("Failed to retrieve data. Status code: ",status_code(response))
    return(data.frame())
  }
  print(status_code(response))
  #
  # Parse the JSON response
  #
  data<-fromJSON(httr::content(response,as="text",encoding="UTF-8"))
  #
  # Extract the components of interest (pollutants and weather details)
  #
  result<-data.frame(date=timestamp,
                     Place=place,
                     Country=country,
                     TAmin=data$temperature$min,
                     TAmax=data$temperature$max,
                     TAmean=mean(data$temperature$morning,data$temperature$afternoon,
                                 data$temperature$evening,data$temperature$night),
                     Pressure=data$pressure$afternoon,
                     Humidity=data$humidity$afternoon,
                     Clouds=data$cloud_cover$afternoon,
                     Precipitation=data$precipitation$total)
  #
  return(result)
  #
}
#
# Retrieve daily pollutants from Open Weather Map API
# 
Get.pollution<-function(timestamp,time.zone,lat,lon,place,country,api_key) {
  #
  url<-paste0(
    "http://api.openweathermap.org/data/2.5/air_pollution/history?",
    "lat=", lat,
    "&lon=", lon,
    "&start=", timestamp,
    "&end=", as.integer(timestamp+60*60*24),
    "&appid=", api_key
  )
  #
  # Make the API request
  #
  response<-GET(url)
  #
  # Check if the request was successful
  #
  if (status_code(response)!=200) {
    print("Failed to retrieve data. Status code: ",status_code(response))
    return(data.frame())
  }
  print(status_code(response))
  #
  # Parse the JSON response
  #
  data<-fromJSON(httr::content(response,as="text",encoding="UTF-8"))
  #
  # Extract the components of interest (pollutants and weather details)
  #
  result<-data.frame(date=timestamp,
                     Place=place,
                     Country=country,
                     AQI=mean(data$list$main$aqi,na.rm=T),
                     CO=mean(data$list$components$co,na.rm=T),
                     NO=mean(data$list$components$no,na.rm=T),
                     NO2=mean(data$list$components$no2,na.rm=T),
                     O3=mean(data$list$components$o3,na.rm=T),
                     SO2=mean(data$list$components$so2,na.rm=T),
                     PM2.5=mean(data$list$components$pm2_5,na.rm=T),
                     PM10=mean(data$list$components$pm10,na.rm=T),
                     NH3=mean(data$list$components$nh3,na.rm=T))
  #
  return(result)
  #
}
#
# Load and edit data from Bresser thermometer
#
ShellyTurbo_func<-function() {
  #
  # Physical Costants
  #
  C3 <- 2898.756 # Wien's displacement constant in um*K
  sigma <- 5.67e-8  # Stefan-Boltzmann constant in W/m^2/K^4
  emissivity <- 0.95  # Emissivity of the filament
  efficiency <- 0.90  # Efficiency of the filament
  #
  # Costants for the filament
  #
  N <- 25  # Number of coils
  D <- 0.01  # Diameter of the filament in meters
  p <- 0.004  # Pitch of the filament in meters
  w <- 0.0027  # width of the filament in meters
  #
  filament_temperature<-function(power){
    #
    # Calculate the filament temperature based on the power input.
    # The formula used is derived from the Stefan-Boltzmann law.
    #
    filament_surface <- N*w*((pi*D)^2 + p^2)^0.5
    temperature_K <- (power*efficiency/(emissivity*sigma*filament_surface))^0.25 
    temperature_C <- temperature_K - 273.15  # Convert to Celsius
    return (temperature_C)  
  }
  #
  filament_total_emissive_power<-function(power){
    #
    # Calculate the filament total emissive power based on the power input.
    #
    filament_surface <- N*w*((pi*D)^2 + p^2)^0.5
    total_emissive_powe <- power * efficiency/filament_surface
    return (total_emissive_powe)  
  }
  #
  filament_peak_lambda<-function(power) {
    #
    # Calculate wavelength at maximum power (Wien's displacement law).
    #
    peak_lambda <- C3 / (filament_temperature(power) + 273.15)  # in um
    return (peak_lambda)  
  }
  #
  # read
  #
  current_dir<-getwd()
  sub_dir<-"ShellyTurbo"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.names<-list.files(folder_path)
  #
  mydata<-read.csv(file.path(folder_path,file.names[1]),header=T,sep=",")
  if (length(file.names)>1) {
    for (i in 2:length(file.name)) mydata<-rbind(mydata,
                                                 read.csv(file.path(folder_path,file.names[i]),header=T,sep=","))
  }
  mydata<-unique(mydata)
  #
  # Editing
  #
  colnames(mydata)<-c("timestamp","Lamp.power","Lamp.voltage")
  #
  # Add epoch
  #
  mydata$timestamp<-substr(mydata$timestamp,1,10)
  unique_ts<-unique(mydata$timestamp)
  mydata$date<-NA
  for (i in 1:length(unique_ts)) {
    index<-which(mydata$timestamp==unique_ts[i])
    mydata$date[index]<-epoch(unique_ts[i],format="DMY")
  }
  #
  # Transform data into numeric variables
  #
  mydata$power<-as.numeric(gsub(",",".",mydata$Lamp.power))
  mydata$voltage<-as.numeric(gsub(",",".",mydata$Lamp.voltage))
  #
  # Calculate wave length at maximum emissive power (Wien law) 
  # and total emissive power (stefan-Boltzmann law)
  #
  mydata$peak_lambda    <- filament_peak_lambda(mydata$power)
  mydata$emissive_power <- filament_total_emissive_power(mydata$power)
  #
  # For each day, indicate mean power, mean voltage, hours of use.
  # Consider that sampling has a frequency of 1/min
  #
  days<-unique(mydata$date)
  vector<-rep(NA,length(days))
  mydata2<-data.frame(date=vector)
  for (i in 1:length(days)) {
    df<-subset.data.frame(mydata,date==days[i])
    df<-subset.data.frame(df,power>50)  
    if (nrow(df)>0) {
      mydata2$Lamp.power[i]<-mean(df$power)
      mydata2$Lamp.voltage[i]<-mean(df$voltage)
      mydata2$Lamp.peak_lambda[i]<-mean(df$peak_lambda)
      mydata2$Lamp.emissive_power[i]<-mean(df$emissive_power)
      mydata2$Lamp.min.on[i]<-nrow(df)
    } else {
      mydata2$Lamp.power[i]<-0
      mydata2$Lamp.voltage[i]<-0
      mydata2$Lamp.peak_lambda[i]<-NA
      mydata2$Lamp.emissive_power[i]<-NA
      mydata2$Lamp.min.on[i]<-nrow(df)
    }
    mydata2$date[i]<-days[i]
  }
  #
  return(mydata2)
} 
#
#---- Wearables ----------------------------------------------------------------
#
# Load and edit data from Bio Impedance Scale.
# Then remove the original file from the folder.
# 
Fitdays<-function() {
  #
  # Read the file exported by the app. Despite the .csv extension it is a
  # legacy Excel (.xls) workbook, so read_xls is used instead of read.csv.
  #
  current_dir <- getwd()
  sub_dir     <- "Fitdays"
  folder_path <- file.path(current_dir,"Data",sub_dir)
  file.name   <- list.files(folder_path)
  file.name   <- file.path(folder_path, file.name)
  BIA <- as.data.frame(readxl::read_xls(file.name[1], col_types = "text"))
  names(BIA) <- make.names(names(BIA))
  #
  # The app writes "Skeletal Muscle" and "Muscle Mass"; keep the old names
  #
  names(BIA)[names(BIA) == "Skeletal.Muscle"] <- "Skeletal.muscle"
  names(BIA)[names(BIA) == "Muscle.Mass"]     <- "Muscle.mass"
  #
  # Add epoch. Date looks like "11:33 Dec.29,2020"
  #
  month_abbrev <- c("Jan","Feb","Mar","Apr","May","Jun",
                    "Jul","Aug","Sep","Oct","Nov","Dec")
  #
  for (i in 1:nrow(BIA)) {
    Month <- substr(BIA$Date[i], 7, 9)
    Month <- match(Month, month_abbrev)
    Year  <- substr(BIA$Date[i], nchar(BIA$Date[i]) - 3, nchar(BIA$Date[i]))
    Day   <- substr(BIA$Date[i], 11, nchar(BIA$Date[i]) - 5)
    BIA$date[i] <- epoch(paste(Month, Day, Year, sep = "/"), format = "MDY")
  }
  #
  # Remove a few columns
  #
  BIA <- subset.data.frame(BIA, select = -Date)
  BIA <- subset.data.frame(BIA, select = -Heart.Rate)
  BIA <- subset.data.frame(BIA, select = -Cardiac.Index)
  BIA <- subset.data.frame(BIA, select = -Body.age)
  #
  # Remove units of measure. Note [^0-9.] and not [^0-9]: the latter deletes
  # the decimal point and turns 90.0kg into 900. Failed readings ("- -")
  # become NA.
  #
  for (j in 1:ncol(BIA)) {
    BIA[, j] <- gsub("[^0-9.]", "", as.character(BIA[, j]))
    BIA[, j][BIA[, j] == "" | BIA[, j] == "."] <- NA
    BIA[, j] <- as.numeric(BIA[, j])
  }
  #
  BIA$Skeletal.muscle.mass  <- BIA$Skeletal.muscle  * BIA$Weight / 100
  BIA$Subcutaneous.fat.mass <- BIA$Subcutaneous.fat * BIA$Weight / 100
  BIA$Body.fat.mass         <- BIA$Body.Fat         * BIA$Weight / 100
  BIA$Visceral.fat.mass     <- BIA$Body.fat.mass - BIA$Subcutaneous.fat.mass
  #
  # If more than one measure is present for the same day, keep only the first
  #
  BIA <- BIA[!duplicated(BIA$date), ]
  #
  for (j in 1:ncol(BIA)) {
    BIA[, j] <- round(BIA[, j], digits = 3)
  }
  rownames(BIA) <- NULL
  #
  return(BIA)
}
#
# Load and edit data from ShapeGanic app (macronitrients)
# 
ShapeGanic<-function() {
  #
  # Read the txt file with tabs 
  #
  current_dir<-getwd()
  sub_dir<-"ShapeGanic"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.name<-list.files(folder_path)
  file.name<-paste0(folder_path,"/",file.name)
  macro<-read.csv(file.name[1],header=T,sep=",")
  if (length(file.name)>1) {
    for (i in 2:length(file.name)) macro<-rbind(macro,read.csv(file.name[i],header=T,sep=","))
  }
  macro<-unique(macro)
  #
  # Add epoch
  #
  for (i in 1:nrow(macro)) {
    macro$date[i]<-epoch(macro$Date[i],format="YMD")
  }
  #
  # Calculate daily values
  #
  n.days<-length(unique(macro$date))
  macro.daily<-data.frame(kcal=rep(NA,n.days),carbo.g=rep(NA,n.days),
                          protein.g=rep(NA,n.days),fat.g=rep(NA,n.days),
                          date=rep(NA,n.days))
  for (i in 1:n.days) {
    df<-subset.data.frame(macro,date==unique(macro$date)[i])
    macro.daily$kcal[i]<-sum(df$Calories...Cal)
    macro.daily$carbo.g[i]<-sum(df$Carbohydrates...g)
    macro.daily$protein.g[i]<-sum(df$Protein...g)
    macro.daily$fat.g[i]<-sum(df$Fats...g)
    macro.daily$date[i]<-df$date[1]
  }
  #
  # Remove days with less than 1000 kcal (likely missing data)
  #
  for (i in 1:nrow(macro.daily)) {
    if (macro.daily$kcal[i]<1000) macro.daily[i,]<-rep(NA,ncol(macro.daily))
  }
  #
  return(macro.daily)
}
#
# Load and edit data from MUSE EEG
# 
EEG_func<-function() {
  #
  # Read
  #
  current_dir<-getwd()
  sub_dir<-"EEG"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file_names<-list.files(folder_path)
  days<-length(file_names)
  #
  # For each file, if uncompressed, compress it
  #
  for (i in 1:days) {
    if (!grepl(".gz",file_names[i])) {
      temp<-fread(file.path(folder_path,file_names[i]),header=T,sep=",",fill=T)
      if (nrow(temp)>1) {
        file.remove(file.path(folder_path,file_names[i]))
        file_names[i]<-paste0(file_names[i],".gz")
        fwrite(temp,file.path(folder_path,file_names[i]))
      }
    }
  }
  #
  # Store each recording in a list, only readings before 'hour_co' and only one for day
  # Remove fNIRS readings
  #
  Channels<-c("RAW_TP9","RAW_AF7","RAW_AF8","RAW_TP10")
  hour_co<-17
  eeg.data<-list()
  k<-1
  for (i in 1:days) {
    day<-substr(file_names[i],13,22)
    hour<-substr(file_names[i],25,26)
    if (as.numeric(hour)<hour_co) { 
      if (!day%in%names(eeg.data)) { # only one reading for day
        eeg.data[[k]]<-setDF(fread(file.path(folder_path,file_names[i]),
                                   select=c("TimeStamp",Channels,"Elements")))
        names(eeg.data)[k]<-substr(eeg.data[[k]][1,1],1,10) # name is the date as YYYY-MM-DD
        k<-k+1
      }
    } 
  }
  days<-length(eeg.data)
  #
  # Arrays for Frequency analysis
  #
  Bands<-list(delta=c(0.5,4),theta=c(4,8),alpha=c(8,13),beta=c(13,30),gamma=c(30,45)) 
  #
  # Prepare a data frame with the readings, one row per day
  #
  v<-rep(NA,days)
  mydata<-data.frame(date=v,
                     TP9_delta=v,TP9_theta=v,TP9_alpha=v,TP9_beta=v,TP9_gamma=v,
                     AF7_delta=v,AF7_theta=v,AF7_alpha=v,AF7_beta=v,AF7_gamma=v,
                     AF8_delta=v,AF8_theta=v,AF8_alpha=v,AF8_beta=v,AF8_gamma=v,
                     TP10_delta=v,TP10_theta=v,TP10_alpha=v,TP10_beta=v,TP10_gamma=v)
  #
  # For each day, edit time stamp
  #
  for (i in 1:days) {
    mydata$date[i]<-epoch(names(eeg.data)[i],"YMD")
    eeg.data[[i]]$TimeStamp<-as.POSIXct(eeg.data[[i]]$TimeStamp,format="%Y-%m-%d %H:%M:%OS",tz="UTC")
    eeg.data[[i]]$TimeStamp<-as.double(eeg.data[[i]]$TimeStamp)
  }
  #
  # Remove first and last seconds 
  # Estimate sampling frequency (Fs)
  # Mark window around event (blink or jaw clench) as NA
  # Remove event empty lines
  #
  Fs<-c()
  for (i in 1:days) {
    #
    # Remove first and last seconds
    #
    start.time<-eeg.data[[i]]$TimeStamp[1]+10
    end.time<-eeg.data[[i]]$TimeStamp[nrow(eeg.data[[i]])]-10
    eeg.data[[i]]<-subset.data.frame(eeg.data[[i]],TimeStamp>=start.time&TimeStamp<=end.time)    
    #
    # Estimate sampling frequency
    #
    df<-eeg.data[[i]] # make a copy
    event<-which(grepl("blink|jaw_clench",eeg.data[[i]]$Elements,ignore.case=T))
    if (length(event)>0)  df<-df[-event,]
    Fs[i]<-(nrow(df)-1)/(max(df$TimeStamp,na.rm=T)-min(df$TimeStamp,na.rm=T))
    #
    # Remove values around events and overwriting NAs
    #
    half_window<-round(Fs[i]*0.5)
    for (ev in event) {
      rows<-max(1,ev-half_window) : min(nrow(eeg.data[[i]]), ev+half_window)
      eeg.data[[i]][rows, Channels] <- NA
    }
    if (length(event)>0)  eeg.data[[i]]<-eeg.data[[i]][-event,]
    rownames(eeg.data[[i]])<-NULL # rename the rows from 1 to the end
    eeg.data[[i]]<-subset.data.frame(eeg.data[[i]],select=-Elements)
  }
  #
  # FFT analysis 
  #
  for (d in 1:days) {
    #
    ns<-10 # interval for FFT in seconds
    RawEEG<-subset.data.frame(eeg.data[[d]],select=Channels)
    num<-floor(nrow(RawEEG)/(ns*round(Fs[d]))) # number of ns-second recordings
    #
    # Work on ns-second intervals of day d
    #
    # Amplitude
    #
    Amp<-list()
    end.i<-0
    for (i in 1:num) {
      start.i<-end.i+1
      end.i<-start.i+ns*round(Fs[d])-1
      RawEEG2<-RawEEG[start.i:end.i,]
      Amp[[i]]<-eegfft(RawEEG2,Fs=Fs[d],lower=Bands[[1]][1],upper=Bands[[length(Bands)]][2]) 
      Amp[[i]]<-data.frame(Amp[[i]]$frequency,Amp[[i]]$strength[,1],Amp[[i]]$strength[,2],
                           Amp[[i]]$strength[,3],Amp[[i]]$strength[,4])
      colnames(Amp[[i]])<-c("frequency",Channels)
    }
    #
    # Relative Power 
    #
    v<-rep(NA,length(Amp))
    mydata2<-data.frame(date=v,
                        TP9_delta=v,TP9_theta=v,TP9_alpha=v,TP9_beta=v,TP9_gamma=v,
                        AF7_delta=v,AF7_theta=v,AF7_alpha=v,AF7_beta=v,AF7_gamma=v,
                        AF8_delta=v,AF8_theta=v,AF8_alpha=v,AF8_beta=v,AF8_gamma=v,
                        TP10_delta=v,TP10_theta=v,TP10_alpha=v,TP10_beta=v,TP10_gamma=v)
    NC<-1
    for (i in 1:length(Amp)) {
      count<-0
      for (j in 1:length(Channels)) {
        temp<-subset.data.frame(Amp[[i]],Amp[[i]]$frequency>=Bands[[1]][1]&Amp[[i]]$frequency<Bands[[length(Bands)]][2])
        TotPower<-trapz(temp$frequency,temp[,(j+1)]^2)
        for (k in 1:length(Bands)) {
          count<-count+1
          temp<-subset.data.frame(Amp[[i]],Amp[[i]]$frequency>=Bands[[k]][1]&Amp[[i]]$frequency<Bands[[k]][2])
          mydata2[i,(NC+count)]<-trapz(temp$frequency,temp[,(j+1)]^2)/TotPower # relative power
        }
      }
    }
    #
    # Calculate the mean values of relative powers
    #
    mydata2[is.nan(as.matrix(mydata2))]<-NA
    mycolumns<-setdiff(colnames(mydata2),"date")  
    for (c in mycolumns) {
      mydata[[c]][d]<-mean(mydata2[[c]],na.rm=T)  
    }
  }
  #
  # Save relative power for each combination band/channel 
  #
  for (i in 1:days) {
    mydata$Delta[i]<-mydata$AF7_delta[i]+mydata$AF8_delta[i]+mydata$TP9_delta[i]+mydata$TP10_delta[i]
    mydata$Delta_L[i]<-mydata$AF7_delta[i]+mydata$TP9_delta[i]
    mydata$Delta_R[i]<-mydata$AF8_delta[i]+mydata$TP10_delta[i]
    mydata$Theta[i]<-mydata$AF7_theta[i]+mydata$AF8_theta[i]+mydata$TP9_theta[i]+mydata$TP10_theta[i]
    mydata$Theta_L[i]<-mydata$AF7_theta[i]+mydata$TP9_theta[i]
    mydata$Theta_R[i]<-mydata$AF8_theta[i]+mydata$TP10_theta[i]
    mydata$Alpha[i]<-mydata$AF7_alpha[i]+mydata$AF8_alpha[i]+mydata$TP9_alpha[i]+mydata$TP10_alpha[i]
    mydata$Alpha_L[i]<-mydata$AF7_alpha[i]+mydata$TP9_alpha[i]
    mydata$Alpha_R[i]<-mydata$AF8_alpha[i]+mydata$TP10_alpha[i]
    mydata$Beta[i]<-mydata$AF7_beta[i]+mydata$AF8_beta[i]+mydata$TP9_beta[i]+mydata$TP10_beta[i]
    mydata$Beta_L[i]<-mydata$AF7_beta[i]+mydata$TP9_beta[i]
    mydata$Beta_R[i]<-mydata$AF8_beta[i]+mydata$TP10_beta[i]
    mydata$Gamma[i]<-mydata$AF7_gamma[i]+mydata$AF8_gamma[i]+mydata$TP9_gamma[i]+mydata$TP10_gamma[i]
    mydata$Gamma_L[i]<-mydata$AF7_gamma[i]+mydata$TP9_gamma[i]
    mydata$Gamma_R[i]<-mydata$AF8_gamma[i]+mydata$TP10_gamma[i]
  }
  #
  # Results
  #
  return(mydata)
}
#
# It plots readings of intensities for each receiver/wavelength for a session 
# provided as input (select = 1), for daily intensities (selection = 2), for 
# daily delta concentrations (selection = 3)
#
fNIRS_plot<-function(recording,select,day=NA) {
  #
  # create the folder for plots if it is not present
  #
  folder_path<-"Output/fNIRS_plotting"
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  # Plot all data
  #
  if (select==1) {
    figure.name<-paste0(day,"all_data.JPG")
    recording<-recording[sample(nrow(recording),min(1000,nrow(recording))),] # sample 
  } else if (select==2) {
    figure.name<-paste0("recordings_all_data.JPG")
  } else if (select==3) {
    figure.name<-paste0("concentrations_all_data.JPG")
  } 
  figure.name<-file.path(folder_path,figure.name)
  jpeg(figure.name,quality=100,res=200,width=2400,height=1200)
  plot(recording,pch=".")
  dev.off()
  #
  # Plot correlation table
  #
  if (select==1) {
    data.select<-day
  } else if (select==2) {
    data.select<-"recordings"
  } else if (select==3) {
    data.select<-"concentrations"
  } 
  Correlations(recording,1,data.select,folder_path,n.cex=10/ncol(recording),alpha=0.05)
} 
#
# It calculates relative concentration of Hb and Hb02 for each receiver. 
# The baseline is built using mean values over all the recordings
#
fNIRS_func<-function(year) {
  #
  # Read
  #
  current_dir<-getwd()
  sub_dir<-"EEG"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file_names<-list.files(folder_path)
  days<-length(file_names)
  #
  # Instructions for column rename, according to https://mind-monitor.com/FAQ.php#footnote-optics
  #
  optics_names<-c(
    Optics1  = "L.Outer.730", Optics2  = "R.Outer.730",
    Optics3  = "L.Outer.850", Optics4  = "R.Outer.850",
    Optics5  = "L.Inner.730", Optics6  = "R.Inner.730",
    Optics7  = "L.Inner.850", Optics8  = "R.Inner.850",
    Optics9  = "L.Outer.660", Optics10 = "R.Outer.660",
    Optics11 = "L.Outer.amb", Optics12 = "R.Outer.amb",
    Optics13 = "L.Inner.660", Optics14 = "R.Inner.660",
    Optics15 = "L.Inner.amb", Optics16 = "R.Inner.amb")
  #
  rename_optics<-function(df) {
    hit<-names(df) %in% names(optics_names)
    names(df)[hit]<-optics_names[names(df)[hit]]
    df
  }
  #
  # For each file, if uncompressed, compress it
  #
  for (i in 1:days) {
    if (!grepl(".gz",file_names[i])) {
      temp<-read.csv(file.path(folder_path,file_names[i]),header=T,sep=",")
      file.remove(file.path(folder_path,file_names[i]))
      file_names[i]<-paste0(file_names[i],".gz")
      fwrite(temp,file.path(folder_path,file_names[i]))
    }
  }
  #
  # Store each recording in a list, only readings before 'hour_co' and only one for day
  #
  hour_co<-17
  fNIRS.data<-list()
  k<-1
  for (i in 1:days) {
    day<-substr(file_names[i],13,22)
    hour<-substr(file_names[i],25,26)
    if (as.numeric(hour)<hour_co) { 
      if (!day%in%names(fNIRS.data)) { # only one reading for day
        temp<-fread(file.path(folder_path,file_names[i]),header=T,sep=",")
        test<-sum(temp$Optics9,na.rm=T)
        name<-substr(temp$TimeStamp[1],1,10) # name is the date as YYYY-MM-DD
        if (test>0) {
          temp<-temp[,c("TimeStamp",paste0("Optics",1:16))] # only optical data
          temp$TimeStamp<-as.numeric(temp$TimeStamp) # convert to number
          #
          # Remove first and last seconds
          #
          start.time<-temp$TimeStamp[1]+10
          end.time<-temp$TimeStamp[nrow(temp)]-10
          temp<-subset.data.frame(temp,TimeStamp>=start.time&TimeStamp<=end.time)    
          #
          fNIRS.data[[k]]<-temp
          names(fNIRS.data)[k]<-name
          k<-k+1  
        }
      }
    } 
  }
  days<-length(fNIRS.data)
  #
  # Plot regressions and correlations of intensities for one recording
  #
  d<-sample(seq(1:days),1) # select a day
  day<-names(fNIRS.data)[d]
  recording<-fNIRS.data[[d]]
  recording<-rename_optics(recording)
  recording$TimeStamp<-seq(1,nrow(recording))
  recording<-as.data.frame(recording)
  select<-1
  fNIRS_plot(recording,select,day)
  #
  # Save daily mean values for each channel
  #
  mymeans<-data.frame(matrix(NA, nrow=days, ncol=ncol(fNIRS.data[[1]])))
  colnames(mymeans)<-colnames(fNIRS.data[[1]])
  colnames(mymeans)[1]<-"date"
  #
  for (i in 1:days) {
    mymeans$date[i]<-epoch(names(fNIRS.data)[i], "YMD")
    for (j in 2:ncol(mymeans)) {
      mymeans[[j]][i]<-mean(fNIRS.data[[i]][[j]], na.rm=T)
    }
  }
  #
  # Rename optics
  #
  mymeans<-rename_optics(mymeans)
  #
  # add rows for missing measures (this is for plotting)
  #
  first.day<-min(mymeans$date)
  last.day<-max(mymeans$date)
  temp<-data.frame(date=seq(first.day,last.day,by=60*60*24))
  mymeans<-merge.data.frame(mymeans,temp,by="date",all.y=T)
  #
  # move date as last column
  #
  v<-mymeans$date
  mymeans$date<-NULL
  mymeans$date<-v
  #
  # Non-overlapping moving mean and remove NaN
  #
  nd<-5
  mydata<-mymeans[seq(nd,nrow(mymeans),by=nd),]
  for (j in 1:ncol(mymeans)) {
    mydata[,j]<-runner(mymeans[,j],k=nd,at=seq(nd,nrow(mymeans),by=nd),
                       f=function(x) mean(x,na.rm=TRUE)) 
  }
  mydata[is.nan(as.matrix(mydata))]<-NA
  #
  # Plot regression of daily mean intensities for all recordings
  #
  temp<-mydata 
  for (j in 1:ncol(mydata)) {
    if (colnames(mydata)[j]!="date")
    x<-mydata[,j]
    out<-boxplot.stats(x)$out
    temp[,j][x %in% out]<-NA
  }
  fNIRS_plot(temp,select=2)  
  #
  # Prepare reference values for the eight receiver/wavelength combinations
  # to use for concentrations calculation
  #
  myaverage<-mymeans[1,]
  for (j in 1:ncol(mymeans)) myaverage[,j]<-mean(mymeans[,j],na.rm=T)
  #
  # Derive differential path factor and paths for each receiver/wavelength
  #
  DPF_func<-function(lambda,age) {
    #
    alpha<-223.3
    beta<-0.05624
    gamma<-0.8493
    delta<--5.723*10^(-7)
    epsilon<-0.001245
    zeta<--0.9025
    #
    DPF<-alpha + beta*age^gamma + delta*lambda^3 + epsilon*lambda^2 + zeta*lambda
    return(DPF)
  }
  #
  # Store the absorption coefficients (1/(cm*mole/liter)) in a matrix
  #
  alpha<-function(C,lambda) {
    if (C=="Hb02"&lambda==730) a<-log(10)*390
    if (C=="Hb"&lambda==730) a<-log(10)*1102.2
    if (C=="Hb02"&lambda==850) a<-log(10)*1058
    if (C=="Hb"&lambda==850) a<-log(10)*691.32
    return(a)
  }
  #
  lambda_1<-730
  lambda_2<-850
  alpha_mat<-matrix(c(alpha("Hb",lambda_1),alpha("Hb",lambda_2),
                      alpha("Hb02",lambda_1),alpha("Hb02",lambda_2)),nrow=2,ncol=2)
  rownames(alpha_mat)<-c(lambda_1,lambda_2)
  colnames(alpha_mat)<-c("Hb","Hb02")
  #
  # Prepare a data frame for the concentrations of Hb and Hb02
  #
  days<-nrow(mymeans)
  v<-rep(NA_real_, days)
  myconc<-data.frame(L.Outer.Hb=v,R.Outer.Hb=v,L.Outer.Hb02=v,R.Outer.Hb02=v,
                     L.Inner.Hb=v,R.Inner.Hb=v,L.Inner.Hb02=v,R.Inner.Hb02=v,
                     date=v)
  #
  # Store concentrations
  #
  for (i in 1:days) {
    #
    # work only on days with data
    #
    if (!is.na(mymeans$L.Outer.730[i])) {
      #
      age<-as.numeric(substr(HRD(mymeans$date[i]),1,4))-year
      #
      for (rc in c(1,2,5,6)) {
        #
        # Store path for each receiver
        #
        mypath<-data.frame(L.Outer.730=3*DPF_func(730,age),
                           R.Outer.730=3*DPF_func(730,age),
                           L.Outer.850=3*DPF_func(850,age),
                           R.Outer.850=3*DPF_func(850,age),
                           L.Inner.730=1*DPF_func(730,age),
                           R.Inner.730=1*DPF_func(730,age),
                           L.Inner.850=1*DPF_func(850,age),
                           R.Inner.850=1*DPF_func(850,age))
        #
        # Calculate concentrations for channel rc
        #
        x.lambda.1<-mypath[,rc] # path for lambda 1 (730 nm)
        x.lambda.2<-mypath[,(rc+2)] # path for lambda 2 (850 nm)
        #
        delta.A.1<-log(myaverage[1,rc]/mymeans[i,rc]) # change in attenuation for lambda 1 (730 nm)
        delta.A.2<-log(myaverage[1,(rc+2)]/mymeans[i,(rc+2)]) # change in attenuation for lambda 2 (850 nm)
        #
        vector<-matrix(c(delta.A.1/x.lambda.1,delta.A.2/x.lambda.2),nrow=2,ncol=1)
        #
        delta.conc<-solve(alpha_mat,vector)
        rownames(delta.conc)<-c("delta.Hb","delta.Hb02")
        #
        # Store results
        #
        myconc$date[i]<-mymeans$date[i] # day
        myconc[i,rc]<-delta.conc[1,1]
        myconc[i,(rc+2)]<-delta.conc[2,1]
      }
    }
  }
  #
  # Add Ambient and Red light
  #
  myconc$L.Outer.amb<-mymeans$L.Outer.amb
  myconc$R.Outer.amb<-mymeans$R.Outer.amb
  myconc$L.Inner.amb<-mymeans$L.Inner.amb
  myconc$R.Inner.amb<-mymeans$R.Inner.amb
  #
  myconc$L.Outer.660<-mymeans$L.Outer.660
  myconc$R.Outer.660<-mymeans$R.Outer.660
  myconc$L.Inner.660<-mymeans$L.Inner.660
  myconc$R.Inner.660<-mymeans$R.Inner.660
  #
  # move date as last column
  #
  v<-myconc$date
  myconc$date<-NULL
  myconc$date<-v
  #
  # Non-overlapping moving mean and remove NaN
  #
  nd<-5
  mydata<-myconc[seq(nd,nrow(myconc),by=nd),]
  for (j in 1:ncol(myconc)) {
    mydata[,j]<-runner(myconc[,j],k=nd,at=seq(nd,nrow(myconc),by=nd),
                       f=function(x) mean(x,na.rm=TRUE)) 
  }
  mydata[is.nan(as.matrix(mydata))]<-NA
  #
  # Plot
  #
  temp<-mydata 
  for (j in 1:ncol(mydata)) {
    if (colnames(mydata)[j]!="date")
      x<-mydata[,j]
    out<-boxplot.stats(x)$out 
    temp[,j][x %in% out]<-NA
  }
  fNIRS_plot(temp,select=3)
  #
  # Add frontal Hb and HB02
  #
  myconc$Outer.Hb<-myconc$L.Outer.Hb+myconc$R.Outer.Hb
  myconc$Outer.Hb02<-myconc$L.Outer.Hb02+myconc$R.Outer.Hb02
  myconc$Inner.Hb<-myconc$L.Inner.Hb+myconc$R.Inner.Hb
  myconc$Inner.Hb02<-myconc$L.Inner.Hb02+myconc$R.Inner.Hb02
  #
  # Add [HbTot] for each receiver (perfusion index)
  #
  myconc$L.Outer.HbTot<-myconc$L.Outer.Hb + myconc$L.Outer.Hb02
  myconc$L.Inner.HbTot<-myconc$L.Inner.Hb + myconc$L.Inner.Hb02
  myconc$R.Outer.HbTot<-myconc$R.Outer.Hb + myconc$R.Outer.Hb02
  myconc$R.Inner.HbTot<-myconc$R.Inner.Hb + myconc$R.Inner.Hb02
  #
  # Add [Hb02] - [Hb] for each receiver (oxygenation index)
  #
  myconc$L.Outer.HbDiff<-myconc$L.Outer.Hb02 - myconc$L.Outer.Hb
  myconc$L.Inner.HbDiff<-myconc$L.Inner.Hb02 - myconc$L.Inner.Hb
  myconc$R.Outer.HbDiff<-myconc$R.Outer.Hb02 - myconc$R.Outer.Hb
  myconc$R.Inner.HbDiff<-myconc$R.Inner.Hb02 - myconc$R.Inner.Hb
  #
  # Add frontal perfusion and oxygenation indeces
  #
  myconc$Outer.HbTot<-myconc$L.Outer.HbTot+myconc$R.Outer.HbTot
  myconc$Outer.HbDiff<-myconc$L.Outer.HbDiff+myconc$R.Outer.HbDiff
  myconc$Inner.HbTot<-myconc$L.Inner.HbTot+myconc$R.Inner.HbTot
  myconc$Inner.HbDiff<-myconc$L.Inner.HbDiff+myconc$R.Inner.HbDiff
  #
  # Results
  #
  return(myconc)
}
#
#---- Activities and blood work ------------------------------------------------
#
# It calculates mechanical work in a day
#
WORK<-function(mydata) {
  #
  joule.step.1<-150 # joule for step when walking
  joule.step.2<-628 # joule for step when running
  watt.bike<-100 # watt for cycling
  joule.wl<-600000 # joule for smallest training session with weights
  #
  len<-nrow(mydata)
  #
  work.a<-rep(NA,len) # aerobic work (walking, running, cycling) in joule
  work.an<-rep(NA,len) # anaerobic work (walking, running, cycling) in joule
  #
  for (i in 1:len) {
    df<-subset.data.frame(mydata,date==mydata$date[i])
    if (!is.na(df$Steps.C[1])&!is.na(df$Bike[1])&!is.na(df$Run[1])) {
      work.a[i]<-0
      work.a[i]<-work.a[i]+df$Steps.C[1]*joule.step.1 # assuming joule.step.1 J/step
      work.a[i]<-work.a[i]+watt.bike*60*df$Bike[1]
      work.a[i]<-work.a[i]+(joule.step.2-joule.step.1)*150*df$Run[1] # assuming 150 steps for minute of running
    }
  }
  #
  for (i in 1:len) {
    df<-subset.data.frame(mydata,date==mydata$date[i])
    if (!is.na(df$WL[1])) {
      work.an[i]<-0
      if (df$WL[1]=="biceps") work.an[i]<-work.an[i]+1*joule.wl 
      if (df$WL[1]=="triceps") work.an[i]<-work.an[i]+1*joule.wl
      if (df$WL[1]=="delts") work.an[i]<-work.an[i]+1*joule.wl
      if (df$WL[1]=="pecs") work.an[i]<-work.an[i]+2*joule.wl
      if (df$WL[1]=="back") work.an[i]<-work.an[i]+3*joule.wl
      if (df$WL[1]=="legs") work.an[i]<-work.an[i]+3*joule.wl
      if (df$WL[1]=="total") work.an[i]<-work.an[i]+2*joule.wl
    }
  }
  #
  mydata$Work.a<-work.a
  mydata$Work.an<-work.an
  return(mydata)
}
#
# Read file from Activity watch. It measures total daily not-afk time
#
AWatch_func<-function() {
  #
  # Read the txt file with tabs 
  #
  current_dir<-getwd()
  sub_dir<-"ActivityWatch"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.name<-list.files(folder_path)
  file.name<-paste0(folder_path,"/",file.name)
  mydata<-read.csv(file.name[1],header=T,sep=",")
  if (length(file.name)>1) {
    for (i in 2:length(file.name)) mydata<-rbind(mydata,read.csv(file.name[i],header=T,sep=","))
  }
  mydata<-unique(mydata)
  #
  # Keep only not-afk durations 
  #
  mydata<-subset.data.frame(mydata,status=="not-afk")
  #
  # Add epoch
  #
  for (i in 1:nrow(mydata)) {
    mydata$date[i]<-epoch(mydata$timestamp[i],format="YMD")
  }
  #
  # For each timestamp keep only the longest measure
  #
  uts<-unique(mydata$timestamp)
  for (i in 1:length(uts)) {
    df<-subset.data.frame(mydata,timestamp==uts[i])
    max_not_afk<-max(df$duration)
    r<-which(mydata$timestamp%in%df$timestamp)
    mydata$duration[r]<-rep(max_not_afk,length(r))
  }
  mydata<-unique(mydata)
  #
  # Calculate daily values
  #
  n.days<-length(unique(mydata$date))
  mydata.daily<-data.frame(notafk=rep(NA,n.days),date=rep(NA,n.days))
  for (i in 1:n.days) {
    df<-subset.data.frame(mydata,date==unique(mydata$date)[i])
    mydata.daily$notafk[i]<-sum(df$duration)
    mydata.daily$date[i]<-df$date[1]
  }
  #
  # Results
  #
  return(mydata.daily)
  #
}
#
# Read and edit results of blood work
# 
Blood.work<-function() {
  #
  # read
  #
  current_dir<-getwd()
  sub_dir<-"BloodWork"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.name<-file.path(folder_path,"BloodWork.txt")
  mydata<-read.csv(file.name,header=T,sep=";")
  #
  # calculate Z scores
  #
  BloodWork<-data.frame(matrix(NA,nrow=nrow(mydata),ncol=1+(ncol(mydata)-1)/3))
  colnames(BloodWork)<-c("date",colnames(mydata)[seq(3,ncol(mydata),by=3)])
  for (i in 1:nrow(mydata)) {
    #
    BloodWork$date[i]<-epoch(mydata$MDY[i],format="MDY")
    #
    for (j in 1:(ncol(BloodWork)-1)) {
      mu<-(mydata[i,j*3+1]+mydata[i,(j*3+1-2)])/2
      range<-(mydata[i,j*3+1]-mydata[i,(j*3+1-2)])
      sigma<-range/6
      BloodWork[i,j+1]<-round((mydata[i,(j*3+1-1)]-mu)/sigma,3)
    }  
  }
  #
  # add T4/T3 
  #
  for (i in 1:nrow(BloodWork)) {
    BloodWork$fT4.fT3[i]<-round((1+BloodWork$fT4[i])/(1+BloodWork$fT3[i]),2)-1 # the normal reference ratio is zero
  }
  #
  # add binary variable for blood draw 
  #
  for (i in 1:nrow(BloodWork)) {
    BloodWork$blood.draw[i]<-1
  }
  #
  # Results
  #
  return(BloodWork)
}
#
#---- Interpolation and plotting -----------------------------------------------
#
# Interpolation with splines (it assumes equally distanced measures and no missing cells)
#
Spline.Me<-function(y) {
  #
  index<-which(is.na(y))
  if (length(index>0)) {
    #
    non_na_indices<-!is.na(y)
    #
    min.y<-min(y,na.rm=T)
    max.y<-max(y,na.rm=T)
    x<-seq(1,length(y),by=1)
    y[!non_na_indices]<-stats::spline(x[non_na_indices],
                                      y[non_na_indices],
                                      xout=x[!non_na_indices])$y
    for (i in 1:length(y)) {
      if (or(y[i]<min.y,y[i]>max.y)) y[i]<-NA
    }  
  }
  return(y)
}
#
# Plot histograms for time domain and frequency domain parameters
# 
Plot_Hist<-function(mydata,folder_path) {
  #
  # plotting time domain histograms
  #
  file.name<-"hist_time_domain.jpeg"
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=750)
  layout.mat<-matrix(c(1,2,3,4,5,6),3,2)
  layout(mat=layout.mat)
  names<-c("avgHR","HRVi","avgRR","pNN50","SDNN","rMSSD")
  bins<-rep(40,length(names))
  xlabels<-c("HR (beat/minute)","HRVi","RR (ms)",
             "pNN50","SDNN (ms)","rMSSD (ms)")
  #
  for (j in 1:length(names)) {
    for (i in 1:ncol(mydata)) {
      if (names[j]==colnames(mydata)[i]) {
        minimo<-min(mydata[,i],na.rm=T)
        massimo<-max(mydata[,i],na.rm=T)
        delta<-(massimo-minimo)/bins[j]
        vector<-seq(from=minimo,to=massimo,by=delta)
        hist(mydata[,i],vector,xlab=xlabels[j],ylab="count",main="")
        abline(v=mean(mydata[,i],na.rm=T),col="red",lwd=2)
      }
    }
  }
  dev.off()
  #
  # plotting histograms for power bands
  #
  file.name<-"hist_frequency_domain.jpeg"
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=750)
  layout.mat<-matrix(c(1,2,3,4,5,6),3,2)
  layout(mat=layout.mat)
  names<-c("LF","HF","LFonHF","TotPower","LFnu","HFnu")
  bins<-rep(40,length(names))
  xlabels<-c("LF (ms^2)","HF (ms^2)","LFonHF","TotPower (ms^2)","LFnu","HFnu")
  #
  for (j in 1:length(names)) {
    for (i in 1:ncol(mydata)) {
      if (names[j]==colnames(mydata)[i]) {
        minimo<-min(mydata[,i],na.rm=T)
        massimo<-max(mydata[,i],na.rm=T)
        delta<-(massimo-minimo)/bins[j]
        vector<-seq(from=minimo,to=massimo,by=delta)
        hist(mydata[,i],vector,xlab=xlabels[j],ylab="count",main="")
        abline(v=mean(mydata[,i],na.rm=T),col="red",lwd=2)
      }
    }
  }
  dev.off() 
  #
  # plotting cardiac histograms
  #
  file.name<-"hist_cardiac.jpeg"
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=750)
  layout.mat<-matrix(c(1,2,3,4),2,2)
  layout(mat=layout.mat)
  names<-c("Sys","Dia","HR","CPI")
  bins<-rep(40,length(names))
  xlabels<-c("Sys (mmHg)","Dia (mmHg)","HR (beat/min)","CPI (W/m^2)")
  #
  for (j in 1:length(names)) {
    for (i in 1:ncol(mydata)) {
      if (names[j]==colnames(mydata)[i]) {
        minimo<-min(mydata[,i],na.rm=T)
        massimo<-max(mydata[,i],na.rm=T)
        delta<-(massimo-minimo)/bins[j]
        vector<-seq(from=minimo,to=massimo,by=delta)
        hist(mydata[,i],vector,xlab=xlabels[j],ylab="count",main="")
        abline(v=mean(mydata[,i],na.rm=T),col="red",lwd=2)
      }
    }
  }
  dev.off() 
  #
  # plotting histograms for daily steps, morning temperature, glucose, ketones, score
  #
  file.name<-"hist_other.jpeg"
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1000,height=750)
  layout.mat<-matrix(c(1,2,3,4,5,4),2,3)
  layout(mat=layout.mat)
  names<-c("Steps.C","Temp","glucose","ketones","Score_m")
  bins<-c(40,20,10,5,8)
  xlabels<-c("steps/day","body temp (°C)","fasting glucose (mg/dL)",
             "fasting ketones (mg/dL)","mean daily score")
  #
  for (j in 1:length(names)) {
    for (i in 1:ncol(mydata)) {
      if (names[j]==colnames(mydata)[i]) {
        minimo<-min(mydata[,i],na.rm=T)
        massimo<-max(mydata[,i],na.rm=T)
        delta<-(massimo-minimo)/bins[j]
        vector<-seq(from=minimo,to=massimo,by=delta)
        hist(mydata[,i],vector,xlab=xlabels[j],ylab="count",main="")
        abline(v=mean(mydata[,i],na.rm=T),col="red",lwd=2)
      }
    }
  }
  dev.off()
}
#
# It plots normalized diagrams
# 
Data.plot<-function(mydata,nd,variables,start,end,folder_path,type,detrend=F) {
  #
  # select variables
  #
  mydata<-subset.data.frame(mydata,select=variables)
  #
  # select first and last day
  #
  if (start==0) start<-"2022-07-09"
  if (end==0) end<-HRD(mydata$date[nrow(mydata)])
  mydata<-subset.data.frame(mydata,mydata$date>=epoch(start,"YMD")&mydata$date<=epoch(end,"YMD"))
  NC<-ncol(mydata)
  L<-nrow(mydata)
  #
  # substitute NaN with NA
  #
  mydata[is.nan(as.matrix(mydata))]<-NA
  #
  # detrending
  #
  if (detrend) {
    t<-seq(1,nrow(mydata),1)
    for (j in 1:ncol(mydata)) {
      if (colnames(mydata)[j]!="date") {
        mydata[,j]<-residuals(lm(mydata[,j]~t,na.action=na.exclude))  
      }
    }
  }
  #
  # non-overlapping mean on blocks of nd days (at least half of the days available);
  # the date of each block is its central date
  #
  at<-seq(nd,nrow(mydata),by=nd)
  mydata<-as.data.frame(lapply(mydata,function(x) {
    runner(x,k=nd,at=at,f=function(v) if (sum(!is.na(v))>=1) mean(v,na.rm=TRUE) else NA)
  }))
  L<-nrow(mydata)
  #
  # scaling (date excluded)
  #
  mydata[,2:NC]<-scale(mydata[,2:NC])
  #
  # back to HRD
  #
  for (i in 1:nrow(mydata)) {
    mydata$date[i]<-HRD(mydata$date[i])
  }
  #
  # linear interpolation (for sparse data)
  #
  for (j in 2:ncol(mydata)) {
    mydata[,j]<-zoo::na.approx(mydata[,j],na.rm=FALSE)  
  }
  #
  # plotting
  #
  file.name<-paste0("normplot_",nd,"_",type,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1500,height=600)
  par(mar = c(5, 4, 4, 8))  # Increase the right margin
  minY<-min(mydata[,c(2:NC)],na.rm=T)
  maxY<-max(mydata[,c(2:NC)],na.rm=T)
  j<-2
  plot(as.Date(mydata$date),mydata[,j],type="l",main="",xlab="",ylab="",
       ylim=c(-3,3),col=j-1,lwd=3,xaxt="n")
  for (j in 3:NC) {
    par(new=T)
    plot(as.Date(mydata$date),mydata[,j],type="l",main="",xlab="",ylab="",
         ylim=c(-3,3),col=j-1,lwd=3,xaxt="n")
  }
  abline(h=seq(round(minY),round(maxY),by=1),lty=2)
  abline(v=seq(as.Date(mydata$date)[1],as.Date(mydata$date)[L],by="month"),lty=2)
  axis(1,at=seq(as.Date(mydata$date)[1],as.Date(mydata$date)[L],by="month"),
       labels=format(seq(as.Date(mydata$date)[1],as.Date(mydata$date)[L],by="month"),"%Y-%m-%d"))
  legend("topright",legend=colnames(mydata)[2:NC],col=seq(1,NC-1),
         lwd=rep(3,(NC-1)),lty=rep(1,(NC-1)),xpd = TRUE,
         inset = c(-0.1, 0))
  dev.off()
}
#
#---- Statistics ---------------------------------------------------------------
#
# Generate a table with means and standard deviations of the results
#
Table_results<-function(mydata,folder_path) {
  #
  names<-c("avgHR","avgRR","SDNN","pNN50","rMSSD",
                     "LF","HF","LFonHF","TotPower","LFnu","HFnu",
                     "Sys","Dia","PP","HR","CPI","Temp",
                     "glucose","ketones","Score_m","Steps","Steps.B","Steps.C")
  names.units<-c("beats/minute","ms","ms","-","ms",
                 "ms^2","ms^2"," ","ms^2","-","-",
                 "mmHg","mmHg","mmHg","beats/minute","W/m^2","C",
                 "mg/dL","mg/dL","-","n/day","n/day","n/day")
  names.range<-c("52-76","785-1160","32-93","-","19-75",
                 "193-1009","83-3630","1.1-11.6","-","-","-",
                 "90-120","60-80","30-50","52-76","0.5-0.7","35.9-36.7",
                 "70-100","<9","-","-","-","-")
  names.col<-c("variable","mean","sd","normal","unit","sample size") 
  #
  results<-matrix(nrow=length(names),ncol=length(names.col))
  colnames(results)<-names.col 
  #
  for (j in 1:length(names)) {
    for (i in 1:ncol(mydata)) {
      if (names[j]==colnames(mydata)[i]) {
        results[j,1]<-names[j]
        results[j,2]<-round(mean(mydata[,i],na.rm=T),digits=2)
        results[j,3]<-round(sd(mydata[,i],na.rm=T),digits=2)
        results[j,4]<-names.range[j]
        results[j,5]<-names.units[j]
        results[j,6]<-sum(!is.na(mydata[,i]))
      }
    }
  }
  #
  # Percentage of good days from 2021
  # 
  results<-rbind(results,rep(NA,ncol(results)))
  mydata$Year<-as.numeric(mydata$Year)
  mydata<-subset.data.frame(mydata,Year>2020)
  den<-sum(!is.na(mydata$Score_m))
  num<-0
  for (i in 1:nrow(mydata)) {
    if(!is.na(mydata$Score_1[i])&mydata$Score_1[i]==6) num<-num+1
    if(!is.na(mydata$Score_2[i])&mydata$Score_2[i]==6) num<-num+1
    if(!is.na(mydata$Score_3[i])&mydata$Score_3[i]==6) num<-num+1
  }
  num<-num/3
  good.days<-round(num/den,2)
  results[j+1,1]<-"good days"
  results[j+1,2]<-good.days*100
  results[j+1,3]<-"-"
  results[j+1,4]<-"-"
  results[j+1,5]<-"%"
  results[j+1,6]<-den
  #
  # Save
  #
  file.name<-"meta_data.csv"
  file.name<-file.path(folder_path,file.name)
  write.table(results,file=file.name,quote=F,row.names=F,col.names=T,sep=",")
  return(results)
}
#
# It calculates and plot the correlation table. Moving mean on nd days must be 
# computed before calling the function!
#
Correlations<-function(mydata,nd,data.select,folder_path,n.cex,alpha=0.05,mtd="spearman",p.cap=16) {
  #
  # Prepare figure name
  #
  file.name<-paste0(data.select,"_corr_",nd,".jpeg")
  file.name<-file.path(folder_path,file.name)
  #
  # Correlation coefficients
  #
  NC<-ncol(mydata)
  R<-cor(mydata[,1:NC],use="pairwise.complete.obs",method=mtd)
  #
  # p-values (only the upper triangle is computed, then mirrored)
  #
  p_mat<-matrix(1,nrow=NC,ncol=NC)
  for (i in 1:(NC-1)) {
    for (j in (i+1):NC) {
      temp<-try(cor.test(mydata[,i],mydata[,j],alternative="two.sided",method=mtd),
                silent=TRUE)
      p_mat[i,j]<-if (inherits(temp,"try-error")) NA else temp$p.value
      p_mat[j,i]<-p_mat[i,j]
    }
  }
  #
  # Correction for multiple comparisons (Benjamini-Hochberg, FDR)
  #
  ut<-upper.tri(p_mat)
  p_adj<-p_mat
  p_adj[ut]<-p.adjust(p_mat[ut],method="BH")
  p_adj[lower.tri(p_adj)]<-t(p_adj)[lower.tri(p_adj)]
  diag(p_adj)<-1
  p_mat<-p_adj
  #
  # Mirroring
  #
  rownames(p_mat)<-colnames(mydata)[1:NC]
  colnames(p_mat)<-colnames(mydata)[1:NC]
  #
  # -log10(p), capped so that one very small p does not flatten the scale
  #
  L<--log10(pmax(p_mat,1e-300))
  L<-pmin(L,p.cap)
  diag(L)<-0
  #
  # Plot: coefficients below the diagonal, significance above it
  #
  jpeg(file.name,quality=100,res=100,width=1000,height=1000)
  par(oma=c(0,2,2,0)) # outer margins (bottom, left, top, right)
  #
  corrplot(L,type="upper",method="circle",is.corr=FALSE,
           col.lim=c(0,p.cap),col=COL1("Reds"),
           addCoef.col="black",number.digits=1,number.cex=n.cex,
           tl.pos="n",cl.pos="r",mar=c(3,7,6,0))
  #
  p_mat[is.na(p_mat)]<-1
  R.plot<-R
  R.plot[is.na(R.plot)]<-0
  #
  corrplot(R.plot,type="lower",method="color",is.corr=TRUE,
           addCoef.col="black",number.digits=2,number.cex=n.cex,
           p.mat=p_mat,sig.level=alpha,insig="blank",cl.pos="b",
           tl.pos="lt",tl.col="black",add=TRUE)
  #
  mtext(paste("Mean on:",nd,"measures"),side=3,outer=TRUE,line=0.5,font=2,cex=1.2)
  #
  dev.off()
  #
  invisible(list(R=R,p=p_mat,logp=L))
}
#
# It plots polynomial regressions and gives a table with results
#
Regressions<-function(mydata,exit.v,nd,data.select,folder_path) {
  #
  # Settings 
  #
  if ("PEM"%in%colnames(mydata)) {
    mydata<-subset.data.frame(mydata,select=-PEM) # we consider only continuous variables  
  }
  NC<-ncol(mydata)
  score<-subset.data.frame(mydata,select=exit.v)
  colnames(score)<-"score"
  index<-which(colnames(mydata)==exit.v)
  mydata<-mydata[,-index]
  #
  for (j in 1:ncol(mydata)) {
    x<-mydata[,j]
    if (mean(x==median(x,na.rm=T),na.rm=T)<0.5) { # skip when more than 50% of values equal the median
      out<-boxplot.stats(x)$out
      mydata[,j][x %in% out]<-NA
    }
  }
  #
  mydata<-cbind(mydata,score) # we use this as dependent variable
  c.W<-4
  r.W<-ceiling((NC-1)/c.W)
  C<-12
  temp_v<-rep(NA,(NC-1)*C)
  R2adj<-matrix(temp_v,nrow=(NC-1),ncol=(C))
  rownames(R2adj)<-colnames(mydata)[1:(NC-1)]
  colnames(R2adj)<-c("1st order R2adj","1st order pval","corrected pval",
                     "2nd order R2adj","2nd order pval","corrected pval",
                     "3rd order R2adj","3rd order pval","corrected pval",
                     "best model","R2adj","cor. p val")
  #
  # Linear regressions 
  #
  file.name<-paste0(data.select,"_linear_",nd,"_",exit.v,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1500,height=2000)
  par(mfrow=c(r.W,c.W))
  for (i in 1:(NC-1)) {
    model<-lm(score~mydata[,i],data=mydata)
    r2adj<-summary(model)$adj.r.squared
    R2adj[i,1]<-r2adj
    TS<-summary(model)$fstatistic[1]
    df1<-summary(model)$fstatistic[2]
    df2<-summary(model)$fstatistic[3]
    pval<-pf(TS,df1=df1,df2=df2,lower.tail=F)
    R2adj[i,2]<-pval
    str_r2adj<-sprintf(paste0("%.",5,"f"),r2adj)
    str_pval<-sprintf("%.2e",pval)
    plot(mydata[,i],mydata$score,pch=21,col="black",bg="violet",
         xlim=c(min(mydata[,i],na.rm=T),max(mydata[,i],na.rm=T)),
         ylim=c(min(mydata$score,na.rm=T),max(mydata$score,na.rm=T)),
         xlab=colnames(mydata)[i],ylab=exit.v,main=paste("R2_adj=",str_r2adj,",p val=",str_pval))
    abline(model,lwd=2)  
  }
  title.str<-paste("linear regressions, average on",nd,"days")
  mtext(title.str,outer=TRUE,line=-1)
  dev.off()
  #
  # Second order regressions
  #
  file.name<-paste0(data.select,"_second_order_",nd,"_",exit.v,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1500,height=2000)
  par(mfrow=c(r.W,c.W))
  for (i in 1:(NC-1)) {
    model<-lm(mydata$score~mydata[,i]+I(mydata[,i]^2))
    r2adj<-summary(model)$adj.r.squared
    R2adj[i,4]<-r2adj
    TS<-summary(model)$fstatistic[1]
    df1<-summary(model)$fstatistic[2]
    df2<-summary(model)$fstatistic[3]
    pval<-pf(TS,df1=df1,df2=df2,lower.tail=F)
    R2adj[i,5]<-pval
    str_r2adj<-sprintf(paste0("%.",5,"f"),r2adj)
    str_pval<-sprintf("%.2e",pval)
    plot(mydata[,i],mydata$score,pch=21,col="black",bg="violet",
         xlim=c(min(mydata[,i],na.rm=T),max(mydata[,i],na.rm=T)),
         ylim=c(min(mydata$score,na.rm=T),max(mydata$score,na.rm=T)),
         xlab=colnames(mydata)[i],ylab=exit.v,main=paste("R2_adj=",str_r2adj,",p val=",str_pval))
    minimo<-min(mydata[,i],na.rm=T)
    massimo<-max(mydata[,i],na.rm=T)
    delta<-(massimo-minimo)/100
    abscissa<-seq(from=minimo,to=massimo,by=delta)
    ordinate<-c()
    B<-model$coefficients
    for (i in 1:length(abscissa)) ordinate[i]<-B[1]+B[2]*abscissa[i]+B[3]*(abscissa[i]^2)
    lines(abscissa,ordinate,lty=1,lwd=2)
  }
  title.str<-paste("second order polynomial regressions, average on",nd,"days")
  mtext(title.str,outer=TRUE,line=-1)
  dev.off()
  #
  # Third order regressions
  #
  file.name<-paste0(data.select,"_third_order_",nd,"_",exit.v,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=100,width=1500,height=2000)
  par(mfrow=c(r.W,c.W))
  for (i in 1:(NC-1)) {
    model<-lm(mydata$score~mydata[,i]+I(mydata[,i]^2)+I(mydata[,i]^3))
    r2adj<-summary(model)$adj.r.squared
    R2adj[i,7]<-r2adj
    TS<-summary(model)$fstatistic[1]
    df1<-summary(model)$fstatistic[2]
    df2<-summary(model)$fstatistic[3]
    pval<-pf(TS,df1=df1,df2=df2,lower.tail=F)
    R2adj[i,8]<-pval
    str_r2adj<-sprintf(paste0("%.",5,"f"),r2adj)
    str_pval<-sprintf("%.2e",pval)
    plot(mydata[,i],mydata$score,pch=21,col="black",bg="violet",
         xlim=c(min(mydata[,i],na.rm=T),max(mydata[,i],na.rm=T)),
         ylim=c(min(mydata$score,na.rm=T),max(mydata$score,na.rm=T)),
         xlab=colnames(mydata)[i],ylab=exit.v,main=paste("R2_adj=",str_r2adj,",p val=",str_pval))
    minimo<-min(mydata[,i],na.rm=T)
    massimo<-max(mydata[,i],na.rm=T)
    delta<-(massimo-minimo)/100
    abscissa<-seq(from=minimo,to=massimo,by=delta)
    ordinate<-c()
    B<-model$coefficients
    for (i in 1:length(abscissa)) ordinate[i]<-B[1]+B[2]*abscissa[i]+B[3]*(abscissa[i]^2)+B[4]*(abscissa[i]^3)
    lines(abscissa,ordinate,lty=1,lwd=2)
  }
  title.str<-paste("third order polynomial regressions, average on",nd,"days")
  mtext(title.str,outer=TRUE,line=-1)
  dev.off()
  #
  # calculate corrected p values
  #
  R2adj[,3]<-p.adjust(R2adj[,2], method = "BH")
  R2adj[,6]<-p.adjust(R2adj[,5], method = "BH")
  R2adj[,9]<-p.adjust(R2adj[,8], method = "BH")
  #
  # best model for each predictor
  #
  strs<-c("linear",NA,NA,"2nd order",NA,NA,"3rd order")
  for (i in 1:(NC-1)) {
    maxR2adj<-max(as.numeric(c(R2adj[i,1],R2adj[i,4],R2adj[i,7])),na.rm=T)
    for (j in c(1,4,7)) {
      if (as.numeric(R2adj[i,j])==maxR2adj) {
        R2adj[i,10]<-strs[j]
        R2adj[i,11]<-maxR2adj
        R2adj[i,12]<-R2adj[i,j+2]
      }
    }
  }
  #
  # plot regression and diagnostics for best model overall
  #
  i.max<-which.max(R2adj[,11])
  var.name<-colnames(mydata)[i.max]
  if (R2adj[i.max,10]=="linear") {
    str<-"linear model"
    formula<-as.formula(paste("score~",var.name))
    model<-lm(formula,data=mydata)
  }
  if (R2adj[i.max,10]=="2nd order") {
    str<-"2nd order polynomial model"
    formula<-as.formula(paste("score~",var.name,"+I(",var.name,"^2)"))
    model<-lm(formula,data=mydata)  
  }
  if (R2adj[i.max,10]=="3rd order") {
    str<-"3rd order polynomial model"
    formula<-as.formula(paste("score~",var.name,"+I(",var.name,"^2)+I(",var.name,"^3)"))
    model<-lm(formula,data=mydata)  
  }
  #
  file.name<-paste0(data.select,"_best_",nd,"_",exit.v,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=150,width=1800,height=900)
  layout.mat<-matrix(c(5,5,1,2,3,4),nrow=2,ncol=3) # column 1: best model; columns 2-3: diagnostics
  layout(mat=layout.mat,widths=c(1.3,1,1))
  par(oma=c(0,0,2,0))
  plot(model,which=1,sub.caption="")
  plot(model,which=2,sub.caption="")
  plot(model,which=3,sub.caption="")
  n<-length(mydata[,1])
  plot(model,which=5,cook.levels=c(4/(n+1),0.5,1),sub.caption="")
  plot(mydata[,i.max],mydata$score,pch=21,col="black",bg="black",cex=1,
       xlim=c(min(mydata[,i.max],na.rm=T),max(mydata[,i.max],na.rm=T)),
       ylim=c(min(mydata$score,na.rm=T),max(mydata$score,na.rm=T)),
       xlab=var.name,ylab=exit.v,main=str)
  minimo<-min(mydata[,i.max],na.rm=T)
  massimo<-max(mydata[,i.max],na.rm=T)
  delta<-(massimo-minimo)/100
  abscissa<-seq(from=minimo,to=massimo,by=delta)
  new.df<-data.frame(abscissa)
  colnames(new.df)<-var.name
  ordinate<-predict(model,newdata=new.df)
  lines(abscissa,ordinate,lty=1,lwd=3)
  title.str<-paste("Best model, average on",nd,"days")
  mtext(title.str,outer=TRUE,line=0.5,font=2)
  dev.off()
  #
  # save results
  #
  file.name<-paste0(data.select,"_R2adj_",nd,"_",exit.v,".csv")
  file.name<-file.path(folder_path,file.name)
  write.csv(R2adj,file=file.name,row.names=T)
}
#
# This function performs PCA analysis and plotting. It requires a variable containing
# a categorical variable with two or more than two levels (0,1,2...)
#
PCA_func<-function(mydata,exit.v,nd,data.select,folder_path) {
  #
  # store exit.v separately
  #
  exit.column<-mydata[[exit.v]]
  mydata<-mydata[,!(names(mydata) %in% exit.v)]
  #
  # scale and PCA
  #
  NC<-ncol(mydata)
  NR<-nrow(mydata)
  mydata.pca<-scale(mydata) # scale data 
  mydata.pca<-princomp(mydata.pca) # do not consider score 
  Var1<-mydata.pca$sdev[1]^2 
  Var2<-mydata.pca$sdev[2]^2 
  Var3<-mydata.pca$sdev[3]^2
  VarTot<-sum(mydata.pca$sdev^2) 
  Var1str<-as.character(round(Var1*100/VarTot,1)) 
  Var2str<-as.character(round(Var2*100/VarTot,1))
  Var3str<-as.character(round(Var3*100/VarTot,1))
  #
  # save loadings
  # 
  loadings<-round(mydata.pca$loadings,4)
  colnames(loadings)<-seq(1,ncol(mydata))
  rownames(loadings)<-colnames(mydata)
  file.name<-file.path(folder_path,paste0("loadings_",data.select,"_",nd,".csv"))
  write.csv(loadings,file=file.name,quote=F,row.names=T)
  #
  # define colors for each level of exit variable
  # 
  unique.exit.v<-sort(unique(as.vector(exit.column)),decreasing=TRUE)
  col.names<-c("green","red","blue","orange","yellow","sienna","purple","pink")
  mycolors<-data.frame(exit.v=unique.exit.v,color=col.names[1:length(unique.exit.v)])
  #
  # PCA plot
  #
  if (T) {
    library(ellipse) # for generating ellipse
    random<-data.frame(PC1=rep(NA,NR),PC2=rep(NA,NR),PC3=rep(NA,NR))
    col.vec<-c()
    pch.vec<-c()
    for (i in 1:NR) {
      df<-subset.data.frame(mycolors,exit.v==exit.column[i])
      col.vec[i]<-df$color
      pch.vec[i]<-21
      random$PC1[i]<-mydata.pca$scores[i,1]
      random$PC2[i]<-mydata.pca$scores[i,2]
      random$PC3[i]<-mydata.pca$scores[i,3]
    }
    file.name<-file.path(folder_path,paste0(data.select,"_","PCA","_",exit.v,"_",nd,".jpeg"))
    jpeg(file.name,quality=100,res=300,width=2000,height=2000)
    par(mfrow=c(2,2))
    #
    # plot plane PC1-PC2
    #
    multy<-5 # multiplier for axes
    minX1<-min(mydata.pca$scores[,1])
    maxX1<-max(mydata.pca$scores[,1])
    minX2<-min(mydata.pca$scores[,2])
    maxX2<-max(mydata.pca$scores[,2])
    plot(mydata.pca$scores[,1],mydata.pca$scores[,2],pch=pch.vec,
         xlim=c(minX1,maxX1),ylim=c(minX2,maxX2),col="black",bg=col.vec,cex=1.5,
         xlab=paste("PC1 - Var =",Var1str,"%"),ylab=paste("PC2 - Var =",Var2str,"%"),asp=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(1,2)]),centre=colMeans(random[index,c(1,2)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    for (j in 1:nrow(mycolors)) {
      mycolors$v[j]<-paste0(mycolors$exit.v[j],"-",mycolors$color[j])
    }
    title(main=paste0(mycolors$v,collapse=" "))
    #
    # plot plane PC1-PC3
    #
    minX3<-min(mydata.pca$scores[,3])
    maxX3<-max(mydata.pca$scores[,3])
    plot(mydata.pca$scores[,1],mydata.pca$scores[,3],pch=pch.vec,
         xlim=c(minX1,maxX1),ylim=c(minX3,maxX3),col="black",bg=col.vec,cex=1.5,
         xlab=paste("PC1 - Var =",Var1str,"%"),ylab=paste("PC3 - Var =",Var3str,"%"),asp=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(1,3)]),centre=colMeans(random[index,c(1,3)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    #
    # plot plane PC2-PC3
    #
    plot(mydata.pca$scores[,2],mydata.pca$scores[,3],pch=pch.vec,
         xlim=c(minX2,maxX2),ylim=c(minX3,maxX3),col="black",bg=col.vec,cex=1.5,
         xlab=paste("PC2 - Var =",Var2str,"%"),ylab=paste("PC3 - Var =",Var3str,"%"),asp=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(2,3)]),centre=colMeans(random[index,c(2,3)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    #
    # 3D plot
    #
    library(scatterplot3d)
    scatterplot3d(mydata.pca$scores[,c(1:3)],pch=pch.vec,color=col.vec,type="p",angle=30)
    dev.off()
  }
  #
  # alternative PCA plot
  #
  if (T) {
    file.name<-file.path(folder_path,paste0(data.select,"_","Biplots","_",exit.v,"_",nd,".jpeg"))
    jpeg(file.name,quality=100,res=300,width=2000,height=2000)
    par(mfrow=c(2,2))
    #
    biplot(mydata.pca,cex=1,col=c("blue"),choices=c(1,2),xlabs=rep("",NR),
           xlab=paste("PC1 - Var =",Var1str,"%"), 
           ylab=paste("PC2 - Var =",Var2str,"%")) 
    points(mydata.pca$scores[1:NR,1],
           mydata.pca$scores[1:NR,2],pch=21,col="black",cex=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(1,2)]),centre=colMeans(random[index,c(1,2)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    #
    biplot(mydata.pca,col=c("blue"),choices=c(1,3),xlabs=rep("",NR),
           xlab=paste("PC1 - Var =",Var1str,"%"), 
           ylab=paste("PC3 - Var =",Var3str,"%"))  
    points(mydata.pca$scores[1:NR,1],
           mydata.pca$scores[1:NR,3],pch=21,col="black",cex=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(1,3)]),centre=colMeans(random[index,c(1,3)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    #
    biplot(mydata.pca,col=c("blue"),choices=c(2,3),xlabs=rep("",NR),
           xlab=paste("PC2 - Var =",Var2str,"%"), 
           ylab=paste("PC3 - Var =",Var3str,"%"))
    points(mydata.pca$scores[1:NR,2],
           mydata.pca$scores[1:NR,3],pch=21,col="black",cex=1)
    for (i in 1:nrow(mycolors)) {
      index<-which(col.vec==mycolors$color[i])
      ellipse_random<-ellipse(cov(random[index,c(2,3)]),centre=colMeans(random[index,c(2,3)]),level=0.95)
      lines(ellipse_random[,1],ellipse_random[,2],col=mycolors$color[i],lwd=3)
    }
    #
    plot(cumsum(mydata.pca$sdev^2)/sum(mydata.pca$sdev^2),type="b",ylim=c(0,1),
         xlab=("number of components"),
         ylab=("cumulative variance/total variance")) 
    segments(1,0.8,NC,0.8,col="red")
    dev.off()
  }
}
#
# It generates a radar chart with binary exit variable (0, 1)
#
Radar<-function(mydata,exit.v,nd,data.select,folder_path) {
  #
  # Store exit.v separately
  #
  exit.column<-mydata[[exit.v]]
  mydata<-mydata[,!(names(mydata) %in% exit.v)]  
  NC<-ncol(mydata)
  NR<-nrow(mydata)
  #
  # Scale
  #
  mydata<-scale(mydata) # scale data 
  mydata<-cbind(mydata,exit.column)
  colnames(mydata)[NC+1]<-"class"
  mydata<-as.data.frame(mydata)
  #
  names.d<-c("Max","Min",paste("High",exit.v),paste("Low",exit.v)) # we need it later 
  #
  # two samples
  #
  mydata.1<-subset.data.frame(mydata,mydata$class==1)   
  mydata.1<-subset.data.frame(mydata.1,select=-class)
  mydata.0<-subset.data.frame(mydata,mydata$class==0) 
  mydata.0<-subset.data.frame(mydata.0,select=-class)
  #
  # p values
  #
  p.value<-matrix(nrow=ncol(mydata.0),ncol=2)
  for (j in 1:ncol(mydata.0)) {
    p.value[j,1]<-colnames(mydata.0)[j]
    p.value[j,2]<-t.test(mydata.0[,j],mydata.1[,j],alternative="two.sided",na.rm=T)$p.value
  }
  p.value[,2]<-as.numeric(p.value[,2])
  #
  # Benjamin-Hochberg correction
  #
  p.value[,2]<-p.adjust(p.value[,2],method="BH")
  #
  # Prepare the input for the chart
  #
  combined_data<-rbind(colMeans(mydata.1,na.rm=T),colMeans(mydata.0,na.rm=T))
  max.data<-max(combined_data)
  min.data<-min(combined_data)
  max_min<-rbind(rep(max.data,ncol(mydata.0)),rep(min.data,ncol(mydata.0)))  # Assuming scale from 1 to 5
  data<-as.data.frame(rbind(max_min,combined_data)) 
  rownames(data)<-names.d
  colors_border<-c("blue","red")
  colors_in<-c("#0000FF50","#FF000050")
  #
  # Create the chart
  #
  file.name<-paste0(data.select,"_radar_",nd,"_",exit.v,".jpeg")
  file.name<-file.path(folder_path,file.name)
  jpeg(file.name,quality=100,res=150,width=2000,height=1000)
  par(mar=rep(5,4),xpd=NA)
  radarchart(data,axistype=1,pcol=colors_border,pfcol=colors_in,plwd=4,
             cglcol="grey",cglty=1,axislabcol="grey",caxislabels=seq(0,5,1),
             vlcex=0.8)
  #
  # Add a legend
  #
  legend(x="topright",legend=names.d[3:4],bty="n",pch=20,col=colors_border,text.col="black",cex=1.2,pt.cex=3)
  #
  # Add p values
  #
  angle<-seq(0,2*pi,by=2*pi/ncol(mydata.0))+pi/2
  label_pos<-apply(data[1:2, ],2,mean)
  alpha<-0.05
  for (i in 1:ncol(data)) {
    text(1.5*cos(angle[i]),1.5*sin(angle[i]),labels=sprintf("%.2e",as.numeric(p.value[i,2])),
         cex=0.8)
  }
  dev.off()
}
#
# It performs Correlations, regressions, radar, and PCA. For time series 
# use data detrended and prewhitened!
#
CoReRaPCA<-function(mydata,exit.v,data.select,nd) {
  #
  # Create folders, if absent
  #
  folder<-"Output/PCA"
  folder_path<-file.path(current_dir,folder)  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  } 
  folder<-"Output/Radar"
  folder_path<-file.path(current_dir,folder)  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  } 
  folder<-"Output/Correlations"
  folder_path<-file.path(current_dir,folder)  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  } 
  folder<-"Output/Regressions"
  folder_path<-file.path(current_dir,folder)  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  } 
  #
  # Remove NaN, interpolate
  #
  for (j in 1:ncol(mydata)) {
    mydata[,j][is.nan(mydata[,j])]<-NA
    mydata[,j]<-zoo::na.approx(mydata[,j],maxgap=5,na.rm=FALSE)
  }
  #
  # Non-overlapping moving mean and remove NaN
  #
  if (nd>1) {
    mydata2<-mydata[seq(nd,nrow(mydata),by=nd),]
    for (j in 1:ncol(mydata)) {
      mydata2[,j]<-runner(mydata[,j],k=nd,at=seq(nd,nrow(mydata),by=nd),
                          f=function(x) mean(x,na.rm=TRUE)) 
    }
    mydata<-mydata2  
  }
  #
  # Detrend and pre-whitening
  #
  t<-seq_len(nrow(mydata))
  for (j in 1:ncol(mydata)) {
    lag1<-dplyr::lag(mydata[,j],1)
    lag2<-dplyr::lag(mydata[,j],2)
    lag3<-dplyr::lag(mydata[,j],3)
    mydata[,j]<-residuals(lm(mydata[,j]~t+lag1+lag2+lag3,na.action=na.exclude))
  }
  #
  # Statistics
  #
  folder<-"Output/Correlations"
  folder_path<-file.path(current_dir,folder)  
  try(Correlations(mydata,nd,data.select,folder_path,15/ncol(mydata),
               alpha=0.05,mtd="spearman",p.cap=16))
  #
  folder<-"Output/Regressions"
  folder_path<-file.path(current_dir,folder)  
  try(Regressions(mydata,exit.v,nd,data.select,folder_path))
  #
  mydata[[exit.v]]<-as.integer(mydata[[exit.v]]>=quantile(mydata[[exit.v]],0.75,na.rm=T))
  #
  folder<-"Output/Radar"
  folder_path<-file.path(current_dir,folder)  
  try(Radar(mydata,exit.v,nd,data.select,folder_path))
  #
  folder<-"Output/PCA"
  folder_path<-file.path(current_dir,folder)  
  mydata<-na.omit(mydata)
  try(PCA_func(mydata,exit.v,nd,data.select,folder_path),silent=T)
}
#
#---- Dynamic regressions ------------------------------------------------------
#
DynaReg<-function(mydata0,exit.vs,covariates,type="marginal",alpha=0.05,
                  data.select=NA,days=days<-c(0,1,2,3,4,8,16,32)) {
  #
  # Create output folders, if not present
  #
  folder<-"Output/Dynamic Regression"
  folder_path<-file.path(current_dir,folder)  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  if (type=="marginal") {
    folder<-"Output/Dynamic Regression/Marginal"
    folder_path<-file.path(current_dir,folder)  
    if(!dir.exists(folder_path)) {
      dir.create(folder_path) 
    }  
  } else if (type=="multinomial") {
    folder<-"Output/Dynamic Regression/Multinomial"
    folder_path<-file.path(current_dir,folder)  
    if(!dir.exists(folder_path)) {
      dir.create(folder_path) 
    } 
    folder_path<-file.path(folder_path,data.select)  
    if(!dir.exists(folder_path)) {
      dir.create(folder_path) 
    } 
  }
  #
  # Work on covariates
  #
  mydata.cov<-subset.data.frame(mydata0,select=covariates)
  #
  # Create the matrix that will contain regression coefficients
  #
  NC<-length(covariates)
  results<-matrix(NA,nrow=NC,ncol=(length(days)+5))
  colnames(results)<-c("variable",days,"days on","days off","size","sd/mean")
  results[,1]<-covariates
  #
  # Add days on, days off, and coefficient of variance 
  #
  for (j in 1:NC) { 
    #
    if (sum(mydata.cov[,j]<0,na.rm=T)==0) {
      results[j,(length(days)+2)]<-sum(mydata.cov[,j]>0,na.rm=T)
      results[j,(length(days)+3)]<-sum(mydata.cov[,j]==0,na.rm=T)  
    }
    results[j,(length(days)+4)]<-sum(!is.na(mydata.cov[,j]))
    #
    if (mean(mydata.cov[,j],na.rm=T)!=0) {
      results[j,(length(days)+5)]<-round(sd(mydata.cov[,j],na.rm=T)/mean(mydata.cov[,j],na.rm=T),2)  
    } else {
      results[j,(length(days)+5)]<-NA
    }
  }
  #
  # Remove NaN, impute missing values, and detrend
  #
  t<-seq(1,nrow(mydata.cov),by=1)
  for (j in 1:ncol(mydata.cov)) {
    mydata.cov[,j][is.nan(mydata.cov[,j])]<-NA
    mydata.cov[,j]<-zoo::na.approx(mydata.cov[,j],maxgap=5,na.rm=FALSE)
    mydata.cov[,j]<-residuals(lm(mydata.cov[,j]~t,na.action=na.exclude))
  }
  #
  # Scale 
  #
  for (j in 1:ncol(mydata.cov)) {
    mydata.cov[,j]<-as.numeric(scale(mydata.cov[,j]))
  }
  #
  # Correction for multiple comparisons (for marginal regression only)
  #
  if (type=="marginal") {
#    alpha<-alpha/ncol(mydata.cov)  
  }
  #
  # Select exit variable
  #
  for (exit.v in exit.vs) {
    #
    # Clear folder from previous run
    #
    exit_path<-file.path(folder_path,exit.v)
    if (dir.exists(exit_path)) unlink(exit_path,recursive=T)
    dir.create(exit_path)
    #
    # Add exit variable to the covariates
    #
    mydata<-mydata.cov
    mydata[[exit.v]]<-mydata0[[exit.v]]
    temp<-setdiff(covariates,exit.v)
    mydata<-subset.data.frame(mydata,select=c(temp,exit.v))
    #
    # Edit matrices for results
    #
    index<-which(results[,1]%in%temp)
    results2<-results[index,]
    #
    # Impute missing values, detrend, pre-white, and scale exit.v
    #
    t<-seq(1,nrow(mydata),by=1)
    mydata[[exit.v]]<-zoo::na.approx(mydata[[exit.v]],maxgap=4,na.rm=FALSE)
    mydata[[exit.v]]<-residuals(lm(mydata[[exit.v]]~t,na.action=na.exclude))
    fit<-ar(mydata[[exit.v]],order.max=30,na.action=na.pass) # autoregression
    mydata[[exit.v]]<-fit$resid # we take residuals of autoregression
    mydata[[exit.v]]<-as.numeric(scale(mydata[[exit.v]]))
    #
    # Convert exit.v into a binary variable and add it
    #
    mydata$exit.v.b<-as.integer(mydata[[exit.v]]>=mean(mydata[[exit.v]],na.rm=T))
    #
    # Regressions 
    #
    for (fam in c("binomial","gaussian")) {
      #
      # Regressions
      #
      myresults<-results2
      myresults.all<-myresults
      p.mat<-myresults
      p.mat.all<-p.mat
      index<-1 # column for the number of days
      for (nd in days) { 
        k<-max(nd,1) # window length: 1 day for nd = 0
        lag<-if (nd==0) 0 else 1 # nd = 0: today; nd >= 1: the nd days before today
        weights<-seq(from=0,to=1,length.out=max(days)) # for the weighted mean
        weights<-weights[(max(days)-k+1):max(days)] # take last k values
        mydata2<-mydata
        #
        # weighted mean and normalization
        #
        for (j in seq_along(temp)) {
          mydata2[,j]<-runner(x=mydata[,j],k=k,lag=lag,f=function(x) {
            if (length(x)<k || all(is.na(x))) return(NA_real_)
            mask<-!is.na(x)
            wx<-weights[mask]
            weighted.mean(x[mask],wx)
          })
          mydata2[,j]<-as.numeric(scale(mydata2[,j])) # normalization
        }
        #
        # Regressions
        #
        for (j in seq_along(temp)) {
          #
          if (type=="marginal") {
            #
            if (fam=="binomial") {
              form<-as.formula(paste("exit.v.b~",colnames(mydata2)[j]))
              model<-try(glm(form,data=mydata2,family=binomial))
              if (inherits(model,"try-error")) {
                message("Skipped: ",exit.v," ~ ",colnames(mydata2)[j],", nd = ",nd," (",trimws(model),")")
                p.mat[j,(index+1)]<-"NE" 
                p.mat.all[j,(index+1)]<-"NE" 
                myresults[j,(index+1)]<-"NE"
                myresults.all[j,(index+1)]<-"NE"
              } 
            } else if (fam=="gaussian") {
              form<-as.formula(paste(exit.v,"~",colnames(mydata2)[j]))
              model<-try(lm(form,data=mydata2))
              if (inherits(model,"try-error")) {
                message("Skipped: ",exit.v," ~ ",colnames(mydata2)[j],", nd = ",nd," (",trimws(model),")")
                p.mat[j,(index+1)]<-"NE" 
                p.mat.all[j,(index+1)]<-"NE" 
                myresults[j,(index+1)]<-"NE"
                myresults.all[j,(index+1)]<-"NE"
              }
            }
            if (is.na(p.mat[j,(index+1)])) {
              cf<-summary(model)$coefficients
              if (nrow(cf)>1) {
                p.mat.all[j,(index+1)]<-cf[2,4] 
                myresults.all[j,(index+1)]<-cf[2,1] 
                if (cf[2,4]<alpha) {
                  p.mat[j,(index+1)]<-cf[2,4] 
                  myresults[j,(index+1)]<-round(cf[2,1],2) 
                } else {
                  p.mat[j,(index+1)]<-"NS"
                  myresults[j,(index+1)]<-"NS" 
                }
              }
            }
            #
          } else if (type=="multinomial") {
            #
            string<-paste(temp,collapse="+")
            #
            if (fam=="binomial") {
              form<-as.formula(paste("exit.v.b~",string))
              model<-try(glm(form,data=mydata2,family=binomial))
              if (inherits(model,"try-error")) {
                message("Skipped: ",exit.v,"multinomial , nd = ",nd," (",trimws(model),")")
                p.mat[j,(index+1)]<-"NE"
                p.mat.all[j,(index+1)]<-"NE" 
                myresults[j,(index+1)]<-"NE"
                myresults.all[j,(index+1)]<-"NE" 
              } 
            } else if (fam=="gaussian") {
              form<-as.formula(paste(exit.v,"~",string))
              model<-try(lm(form,data=mydata2))
              if (inherits(model,"try-error")) {
                message("Skipped: ",exit.v,"multinomial , nd = ",nd," (",trimws(model),")")
                p.mat[j,(index+1)]<-"NE"
                p.mat.all[j,(index+1)]<-"NE" 
                myresults[j,(index+1)]<-"NE"
                myresults.all[j,(index+1)]<-"NE" 
              } 
            }
            if (is.na(p.mat[j,(index+1)])) {
              cf<-summary(model)$coefficients
              a<-which(rownames(cf)==colnames(mydata2)[j])
              if (length(a)!=0) {
                p.mat.all[j,(index+1)]<-cf[a,4] 
                myresults.all[j,(index+1)]<-cf[a,1] 
                if (cf[a,4]<alpha) {
                  p.mat[j,(index+1)]<-cf[a,4] 
                  myresults[j,(index+1)]<-round(cf[a,1],2) 
                } else {
                  p.mat[j,(index+1)]<-"NS" 
                  myresults[j,(index+1)]<-"NS" 
                }
              }
            }
          }
          #
          # Plot significant regressions
          #
          if(!dir.exists(file.path(folder_path,exit.v,"plots"))) {
            dir.create(file.path(folder_path,exit.v,"plots")) 
          }
          v<-p.mat[j,(index+1)]
          if (!is.na(v) && !(v %in% c("NS","NE"))) {
            p1<-as.numeric(p.mat[j,(index+1)])
            if (p1<alpha) {
              file.name<-paste0(exit.v,"_",colnames(mydata2)[j],"_",nd,"_",fam,".png")
              file.name<-file.path(folder_path,exit.v,"plots",file.name)
              png(file.name,width=800,height=400)
              a<-try(crPlots(model,terms=colnames(mydata2)[j],ylab=exit.v),silent=T)
              if (!inherits(a,"try-error")) {
                title(main=paste0("Dynamic regression for ",nd," days (",fam,")"))
              }
              dev.off()
              if (inherits(a,"try-error")) unlink(file.name)
            }  
          }
        }
        index<-index+1
      }
      #
      # For multinomial regression, plot signed -log10 p values
      # (above zero: positive coefficient, below zero: negative coefficient)
      #
      if (type=="multinomial") {
        cols<-2:(length(days)+1)
        pv<-suppressWarnings(matrix(as.numeric(p.mat.all[,cols]),nrow=nrow(p.mat.all)))
        cf.all<-suppressWarnings(matrix(as.numeric(myresults.all[,cols]),nrow=nrow(myresults.all)))
        p.cap<-16 # maximum -log10(p) shown
        lp<- -log10(pmax(pv,10^-p.cap))*sign(cf.all)
        nr<-nrow(lp)
        n.col<-ceiling(sqrt(nr))
        n.row<-ceiling(nr/n.col)
        file.name<-file.path(folder_path,exit.v,paste0(type,"_log10p_",exit.v,"_",fam,".png"))
        png(file.name,width=4*n.col,height=3.3*n.row,units="in",res=300)
        par(mfrow=c(n.row,n.col),mar=c(4,4,3,1),oma=c(0,0,3,0))
        for (i in 1:nr) {
          h<-lp[i,]
          h[is.na(h)]<-0 # models not estimated: no bar
          ym<-max(c(abs(h),-log10(alpha)),na.rm=TRUE)*1.1 # own scale for each panel, 10% margin
          barplot(h,names.arg=days,ylim=c(-ym,ym),main=p.mat.all[i,1],
                  xlab="days",ylab="signed -log10(p)",col=ifelse(h>=0,"steelblue","firebrick"))
          abline(h=0)
          abline(h=c(-1,1)*-log10(alpha),col="red",lty=2)
        }
        mtext(paste0(exit.v,"~",data.select," (",fam,")"),side=3,outer=TRUE,line=1,cex=1.3,font=2)
        dev.off()
      }
      #
      # Save
      #
      file.name<-paste0(type,"_reg_coeff_",exit.v,"_",fam,".csv")
      file.name<-file.path(folder_path,exit.v,file.name)
      write.table(myresults,file=file.name,quote=F,row.names=F,col.names=T,sep=",")
      #
      file.name<-paste0(type,"_pval_",exit.v,"_",fam,".csv")
      file.name<-file.path(folder_path,exit.v,file.name)
      write.table(p.mat,file=file.name,quote=F,row.names=F,col.names=T,sep=",")
    }
  }
}