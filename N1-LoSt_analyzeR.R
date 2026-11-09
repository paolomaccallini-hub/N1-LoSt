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
#---- Study of exit variables --------------------------------------------------
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
# Run decomposition, autocorrelation, Holt-Winters forecast and FFT
#
Outcome_analysis(mydata,exit.v="Score_m",start="2021-03-01",maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="Score_T",start="2021-03-01",maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="notafk",start=span(mydata,"notafk")[1],maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="Steps.C",start=span(mydata,"Steps.C")[1],maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="Alpha",start=span(mydata,"Alpha")[1],maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="Sys",start="2023-06-29",maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="Dia",start="2023-06-29",maxgap=10,
                 window=3*365,n.harm=10)
Outcome_analysis(mydata,exit.v="HR",start="2023-06-29",maxgap=10,
                 window=3*365,n.harm=10)
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
##--- Good days by year -------------------------------------------------------
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
##--- General plotting ---------------------------------------------------------
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
deg<-4
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
variables<-c("Sys","Sys.1","Sys.2","Sys.3","Sys.4","Sys.5","Sys.6","Sys.7",
             "Sys.8","Sys.9","Sys.10")
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
variables<-c("Dia","Dia.1","Dia.2","Dia.3","Dia.4","Dia.5","Dia.6","Dia.7",
             "Dia.8","Dia.9","Dia.10")
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
# Plotting: PP
#
variables<-c("PP","PP.1","PP.2","PP.3","PP.4","PP.5","PP.6","PP.7","PP.8","PP.9","PP.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Pulse Pressure (Sys - Dia)"
deg<-2
LeanTest_plot(mydata,data.select,deg)  
#
# Plotting: MAP
#
variables<-c("MAP","MAP.1","MAP.2","MAP.3","MAP.4","MAP.5","MAP.6","MAP.7",
             "MAP.8","MAP.9","MAP.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Mean Arterial Pressure (mmHg)"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
# Plotting: SV
#
variables<-c("SV","SV.1","SV.2","SV.3","SV.4","SV.5","SV.6","SV.7","SV.8",
             "SV.9","SV.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Stroke Volume"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
# Plotting: CO
#
variables<-c("COut","COut.1","COut.2","COut.3","COut.4","COut.5","COut.6","COut.7","COut.8","COut.9","COut.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Cardiac Output"
deg<-2
LeanTest_plot(mydata,data.select,deg)   
#
# Plotting: CP
#
variables<-c("CP","CP.1","CP.2","CP.3","CP.4","CP.5","CP.6","CP.7","CP.8","CP.9","CP.10")
mydata<-subset.data.frame(mydata0,select=variables)
data.select="Cardiac Power"
deg<-2
LeanTest_plot(mydata,data.select,deg)     
#
##--- Study of positivity to POTS and/or OI ------------------------------------
#
# Function for plotting studies of Lean Tests by minute standing
#
pots.bar<-function(y,lab,ylab,main,file.name) {
  jpeg(file.path("Output","LeanTest",file.name),quality=100,res=100,width=800,height=800)
  ymax<-max(y,na.rm=TRUE)+10
  bp<-barplot(y,names.arg=1:10,xlab="minute",ylab=ylab,ylim=c(0,ymax),
              col=NA,border=NA) # empty plot: axes and positions only
  abline(h=seq(0,round(ymax),by=5),col="grey85") # horizontal grid
  barplot(y,ylim=c(0,ymax),add=TRUE,axes=FALSE,axisnames=FALSE) # bars over the grid
  text(bp,y+2,labels=lab,col="red")
  title(main=main,cex.main=0.9)
  dev.off()
}
#
# Number of POTS positive lean tests by minute of standing
#
variables<-c("Delta.HR.1","Delta.HR.2","Delta.HR.3","Delta.HR.4","Delta.HR.5",
             "Delta.HR.6","Delta.HR.7","Delta.HR.8","Delta.HR.9","Delta.HR.10")
mydata<-subset.data.frame(mydata0,select=variables)
x<-as.matrix(mydata)
last<-apply(!is.na(x),1,function(v) if (any(v)) max(which(v)) else 0) # test duration (minutes)
run.max<-t(apply(x,1,function(v) cummax(ifelse(is.na(v),-Inf,v)))) # max Delta HR up to each minute
count.TOT<-sapply(1:10,function(m) sum(last>=m)) # tests lasting at least m minutes
count.POS<-sapply(1:10,function(m) sum(last>=m & run.max[,m]>=30)) # ...and positive by minute m
count.PER<-round(100*count.POS/count.TOT)
Lab<-paste0(count.POS,"/",count.TOT)
#
pots.bar(count.PER,Lab,"Max Delta HR >= 30 bpm (%)",
         "N of positive lean tests at minute m / total available lean tests at minute m",
         "POTS_positive_by_duration.jpeg")
#
# Sensitivity of measure at each minute for positive lean tests of 10 minutes
#
mydata<-subset.data.frame(mydata0,!is.na(Delta.HR.10)) # only 10-minute lean tests
mydata<-subset.data.frame(mydata,max.Delta.HR>=30) # only positive 10-minute lean tests
mydata<-subset.data.frame(mydata,select=variables)
x<-as.matrix(mydata)
last<-apply(!is.na(x),1,function(v) if (any(v)) max(which(v)) else 0) # test duration (minutes)
run.max<-t(apply(x,1,function(v) cummax(ifelse(is.na(v),-Inf,v)))) # max Delta HR up to each minute
count.POS<-sapply(1:10,function(m) sum(last>=m & run.max[,m]>=30)) # ...and positive by minute m
count.TOT<-rep(nrow(mydata),10)
count.PER<-round(100*count.POS/count.TOT)
Lab<-paste0(count.POS,"/",count.TOT)
#
pots.bar(count.PER,Lab,"Max Delta HR >= 30 bpm (%)",
         "N of positive lean tests at minute m / total positive 10-minute lean tests",
         "POTS_sensitivity_by_minute.jpeg")
#
# Number of OH positive lean tests by minute of standing
#
xs<-as.matrix(mydata0[paste0("Delta.Sys.",1:10)]) # change in systolic pressure
xd<-as.matrix(mydata0[paste0("Delta.Dia.",1:10)]) # change in diastolic pressure
has<-!is.na(xs) | !is.na(xd) # minutes with at least one reading
last<-apply(has,1,function(v) if (any(v)) max(which(v)) else 0) # test duration (minutes)
run.min<-function(x) t(apply(x,1,function(v) cummin(ifelse(is.na(v),Inf,v)))) # min up to each minute
min.s<-run.min(xs)
min.d<-run.min(xd)
pos<-(min.s<=-20) | (min.d<=-10) # OH criterion reached by each minute
count.TOT<-sapply(1:10,function(m) sum(last>=m)) # tests lasting at least m minutes
count.POS<-sapply(1:10,function(m) sum(last>=m & pos[,m])) # ...and positive by minute m
count.PER<-round(100*count.POS/count.TOT)
Lab<-paste0(count.POS,"/",count.TOT)
#
pots.bar(count.PER,Lab,"Min Delta Sys <= -20 or Min Delta Dia <= -10 mmHg or (%)",
         "N of positive lean tests at minute m / total available lean tests at minute m",
         "OH_sensitivity_by_minute.jpeg")
#
# Number of OI positive lean tests by minute of standing (OI = OH or POTS)
#
xs<-as.matrix(mydata0[paste0("Delta.Sys.",1:10)]) # change in systolic pressure
xd<-as.matrix(mydata0[paste0("Delta.Dia.",1:10)]) # change in diastolic pressure
xr<-as.matrix(mydata0[paste0("Delta.HR.",1:10)]) # change in heart rate
has<-!is.na(xs) | !is.na(xd) | !is.na(xr) # minutes with at least one reading
last<-apply(has,1,function(v) if (any(v)) max(which(v)) else 0) # test duration (minutes)
run.min<-function(x) t(apply(x,1,function(v) cummin(ifelse(is.na(v),Inf,v)))) # min up to each minute
run.max<-function(x) t(apply(x,1,function(v) cummax(ifelse(is.na(v),-Inf,v)))) # max up to each minute
min.s<-run.min(xs)
min.d<-run.min(xd)
max.r<-run.max(xr)
pos<-(min.s<=-20) | (min.d<=-10) | (max.r>=30) # OI criterion reached by each minute
count.TOT<-sapply(1:10,function(m) sum(last>=m)) # tests lasting at least m minutes
count.POS<-sapply(1:10,function(m) sum(last>=m & pos[,m])) # ...and positive by minute m
count.PER<-round(100*count.POS/count.TOT)
Lab<-paste0(count.POS,"/",count.TOT)
#
pots.bar(count.PER,Lab,"OH or POTS criterion met (%)",
         "N of positive lean tests by minute m / total lean tests lasting at least m minutes",
         "OI_positive_by_duration.jpeg")
#
# POTS and OH positivity of each lean test, and their association
#
mydata<-subset.data.frame(mydata0,!is.na(HR.10)) # only 10-minute lean tests
m.max<-10 # minutes considered (e.g. 6 to compare all tests on the same duration)
xh<-as.matrix(mydata[paste0("Delta.HR.",1:m.max)]) # change in heart rate
xs<-as.matrix(mydata[paste0("Delta.Sys.",1:m.max)]) # change in systolic pressure
xd<-as.matrix(mydata[paste0("Delta.Dia.",1:m.max)]) # change in diastolic pressure
#
has<-!is.na(xh) | !is.na(xs) | !is.na(xd) # minutes with at least one reading
done<-rowSums(has)>0 # days with a lean test
#
pots.any<-apply(xh,1,function(v) any(v>=30,na.rm=TRUE)) # HR rise >= 30 bpm at any minute
oh.any<-apply(xs,1,function(v) any(v<=-20,na.rm=TRUE)) | # systolic drop >= 20 mmHg
  apply(xd,1,function(v) any(v<=-10,na.rm=TRUE))   # or diastolic drop >= 10 mmHg
#
tab<-table(POTS=pots.any[done],OH=oh.any[done])
print(tab) # counts
print(round(prop.table(tab,1),2)) # share of OH-positive tests by POTS status
print(fisher.test(tab))
#
##--- Search for a predictor of LTs of 10 minutes from shorter LTs -------------
#
# type: slope, ratio.max, max, ratio.mean, mean
#
Find_Error_func<-function(mydata0,type) {
  #
  mydata<-subset.data.frame(mydata0,!is.na(Delta.HR.10)) # only 10-minute lean tests
  mydata<-mydata[paste0("Delta.HR.",1:10)]
  #
  # Build the regression coefficients for duration from 3 minutes to ten
  #
  for (m in 3:10) {
    mydata[[paste0("Delta.HR.",m,".",type)]]<-NA
    for (i in 1:nrow(mydata)) {
      v<-as.vector(transpose(mydata[i,1:m]))
      #
      if (type=="slope") {
        df<-data.frame(time=seq(1,m,1),Delta.HR=v$V1)
        model<-try(lm(Delta.HR~time,data=df),silent=T)
        if (is.list(model)) {
          beta<-model$coefficients[2]
          if (!is.na(beta)) {
            mydata[[paste0("Delta.HR.",m,".",type)]][i]<-beta
          }
        }  
      } else if (type=="max.ratio") {
        if (!is.infinite(max(v$V1,na.rm=T))) {
          mydata[[paste0("Delta.HR.",m,".",type)]][i]<-max(v$V1,na.rm=T)/m
        } 
      } else if (type=="max") {
        if (!is.infinite(max(v$V1,na.rm=T))) {
          mydata[[paste0("Delta.HR.",m,".",type)]][i]<-max(v$V1,na.rm=T)
        } 
      } else if (type=="mean.ratio") {
        if (!is.infinite(max(v$V1,na.rm=T))) {
          mydata[[paste0("Delta.HR.",m,".",type)]][i]<-mean(v$V1,na.rm=T)/m
        } 
      } else if (type=="mean") {
        if (!is.infinite(max(v$V1,na.rm=T))) {
          mydata[[paste0("Delta.HR.",m,".",type)]][i]<-mean(v$V1,na.rm=T)
        } 
      }
    }
  }
  #
  # Calculate absolute relative error
  #
  mydata<-mydata[paste0("Delta.HR.",3:10,".",type)] 
  for (m in 3:9) {
    cname<-paste0("Delta.HR.",m,".",type)
    mydata[[cname]]<-round(100*abs(mydata[[cname]]-mydata[[paste0("Delta.HR.10.",type)]])/abs(mydata[[paste0("Delta.HR.10",".",type)]]))
  }
  for (j in 1:ncol(mydata)) {
    x<-mydata[,j]
    out<-boxplot.stats(x)$out # outliers
    mydata[,j][x %in% out]<-NA  # remove outliers
  }
  #
  # Plot the errors
  #
  file.name<-paste0("Delta_HR_",type,"_RAE_by_minute_10_min_LT.jpeg")
  jpeg(file.path("Output","LeanTest",file.name),quality=100,res=100,width=800,height=800)
  max.err<-quantile(mydata[paste0("Delta.HR.",3:9,".",type)],na.rm=T)
  bp<-boxplot(mydata[paste0("Delta.HR.",3:9,".",type)],names=seq(3,9),xlab="Minutes standing",
              ylab=paste0("Delta.HR ",type," RAE (%)"),ylim=c(0,max.err[4]))
  abline(h=c(10,20),col="red",lty=2,lwd=2)
  dev.off()
  #
}
#
# Is the regression coefficient of Delta.HR of 10-minute Lean test predicted 
# by a test of shorter duration?
#
Find_Error_func(mydata0,type="slope")
#
# Is max.Delta.HR/duration of 10-minute Lean test predicted 
# by a test of shorter duration?
#
Find_Error_func(mydata0,type="max.ratio")
#
# Is max.Delta.HR of 10-minute Lean test predicted 
# by a test of shorter duration?
#
Find_Error_func(mydata0,type="max")
#
# Is mean.Delta.HR/duration of 10-minute Lean test predicted 
# by a test of shorter duration?
#
Find_Error_func(mydata0,type="mean.ratio")
#
# Is mean.Delta.HR of 10-minute Lean test predicted 
# by a test of shorter duration?
#
Find_Error_func(mydata0,type="mean")
#
##---- Time series -------------------------------------------------------------
#
variables<-c("date","Score_m","Max.Delta.HR","Min.Delta.Sys","Steps.C")
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
nd<-7
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Climate_BG"
exit.v<-"Score_m"
variables<-c(exit.v,c("Temp","MRT","TG","TA","AH","RH","TAmean","Pressure","Light",
                      "Clouds","Precipitation","PM10","PM2.5","SO2","NO","AQI","CO"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-7
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain_L"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF","HRV.hour",
                      "Delta_L","Theta_L","Alpha_L","Beta_L","Gamma_L","EEG.hour",
                      "L.Outer.HbDiff","L.Outer.HbTot","fNIRS.hour"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-5
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain_R"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF",
                      "Delta_R","Theta_R","Alpha_R","Beta_R","Gamma_R","EEG.hour",
                      "R.Outer.HbDiff","R.Outer.HbTot","fNIRS.hour"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-5
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Brain"
exit.v<-"Score_m"
variables<-c(exit.v,c("Sys","Dia","SpO2","Temp",
                      "SDNN","pNN50","rMSSD","LFnu","HFnu","LFonHF","HRV.hour",
                      "Delta","Theta","Alpha","Beta","Gamma","EEG.hour",
                      "Outer.HbDiff","Outer.HbTot","fNIRS.hour"))
mydata<-subset.data.frame(mydata0,select=variables)
nd<-7
CoReRaPCA(mydata,exit.v,data.select,nd)
#
# Select covariates
#
data.select<-"Lean Test"
exit.v<-"notafk"
variables<-c(exit.v,c("HR","Delta.HR.1","Delta.HR.2","Delta.HR.3","Delta.HR.4","Delta.HR.5",
                      "Delta.HR.6","Delta.HR.7","Delta.HR.8","Delta.HR.9","Delta.HR.10","mean.Delta.HR","LT.hour",
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
# Select model 
#
evs<-c("Score_m","notafk","Score_T","Alpha")
dsel<-"Climate_BG"
cv<-c("MRT","AH","TAmean","Pressure","Light","Clouds","Precipitation","AQI","Steps.C")
tp<-"multinomial"
al<-0.05
ds<-c(0)
fm<-c("binomial","gaussian")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm,unit.t=4)
#
# Select model
#
evs<-c("Score_m","Score_T","max.Delta.HR","notafk","Alpha")
dsel<-"fNIRS"
cv<-c("Outer.Hb02","Outer.HbTot",
      "Inner.Hb02","Inner.HbTot",
      "Temp","SpO2","TA.in")
tp<-"multinomial"
al<-0.05
ds<-c(0)
fm<-c("binomial","gaussian")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm,unit.t=3)
#
# Select model
#
evs<-c("Score_m")
dsel<-"fNIRS_L"
cv<-c("L.Outer.Hb02","L.Outer.HbTot",
      "L.Inner.Hb02","L.Inner.HbTot",
      "Temp","SpO2")
tp<-"multinomial"
al<-0.05
ds<-c(0,1)
fm<-c("binomial","gaussian")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select model
#
evs<-c("Score_m")
dsel<-"fNIRS_R"
cv<-c("R.Outer.Hb02","R.Outer.HbTot",
      "R.Inner.Hb02","R.Inner.HbTot",
      "Temp","SpO2")
tp<-"multinomial"
al<-0.05
ds<-c(0,1)
fm<-c("binomial","gaussian")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select model
#
evs<-c("Score_m","Score_T","max.Delta.HR","notafk")
dsel<-"EEG"
cv<-c("Theta","Alpha","Beta","Gamma")
tp<-"multinomial"
al<-0.05
ds<-c(0)
fm<-c("binomial","gaussian")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm,unit.t=1)
#
# Select model
#
evs<-c("Score_m","Score_T","Steps.C")
dsel<-"Drugs"
cv<-c("Provigil","Bupropion",
      "Aniracetam","D_serine","Sarcosine","Destrometorfano",
      "Beta.Lactam",
      "Tachipirine","Ibuprofene",
      "Corticosteroid","Pregnenolone","DHEA","Melatonin","T4","T3",
      "Ketotifen","Desloratadina",
      "Pioglitazone","Deursil",
      "Mestinon","Etilefrine","Acetazolamide",
      "Nimodipine","Compression.suit","TAmean","Steps.C")
tp<-"multinomial"
al<-0.05
ds<-c(0,1,3,7,14,21,28)
fm<-c("gaussian","binomial")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select model
#
evs<-c("Score_m","Score_T","Steps.C")
dsel<-"Supplements"
cv<-c("Creatine","ALCar",
      "VitE","VitD3","B12.Tot",
      "Thiamin","Riboflavin","Niacin","VitB6","Folic.acid",          
      "Folinic.acid","Biotin",   
      "Pantothenic.acid","Choline.bitartrate","Myo.Inositol",  
      "Omega3","Steps.C")
tp<-"multinomial"
al<-0.05
ds<-c(0,1,3,7,14,21,28)
fm<-c("gaussian","binomial")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select model
#
evs<-c("Score_T","Score_m","notafk","SDNN","max.Delta.HR","Outer.Hb02","Alpha")
dsel<-"Exercise"
cv<-c("Steps.C*Bike","TAmean","TA","Provigil")
tp<-"multinomial"
al<-0.05
ds<-c(1)
fm<-c("gaussian","binomial")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select covariates 
#
evs<-c("Score_m","Score_T","max.Delta.HR","notafk")
dsel<-"HRV"
cv<-c("SDNN","pNN50","rMSSD","LFonHF","Steps.C")
tp<-"multinomial"
al<-0.05
ds<-c(0,1)
fm<-c("gaussian","binomial")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
#
# Select covariates 
#
evs<-c("Score_m")
dsel<-"Lean test"
cv<-c("max.Delta.HR","ST","LT.hour","Steps.C","TA.in")
tp<-"multinomial"
al<-0.05
ds<-c(0)
fm<-c("gaussian","binomial")
DynaReg(mydata0,exit.vs=evs,covariates=cv,type=tp,alpha=al,data.select=dsel,days=ds,family=fm)
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
