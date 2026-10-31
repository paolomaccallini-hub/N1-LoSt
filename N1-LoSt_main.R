# file name: N1-LoSt_main
#
# Rome, 19th June 2024
#
# Main unit of the Monster Data Study: it builds the file with all the daily 
# variables: "all_daily_data.csv"
#
# Sections can be collapsed in RStudio (Alt+O collapses all, Shift+Alt+O expands all)
#
#---- Directory and subroutines ------------------------------------------------
#
current_dir<-getwd()
source("N1-LoSt_func.R") # Compile subroutines
set.seed(12)
#
# create output folder for plots if it is not present
#
folder_path<-"Output"
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
#---- Read score, drugs, etc and edit ---------------------------------------------------
#
##--- Read ---------------------------------------------------------------------
#
# Read data as strings to avoid a few problems later
#
Score<-read_excel("Data/N1-LoSt_score.xlsx",col_names=T,col_types="guess")
number<-ncol(Score)
Score<-read_excel("Data/N1-LoSt_score.xlsx",col_names=T,col_types=rep("text",number))
Score<-Score[1:(nrow(Score)-3),]
#
##--- Add a column with date as epoch ------------------------------------------
#
Score$date<-as.Date(paste(Score$Year,Score$Month,Score$Day,sep="-"))
Score$date<-as.character(Score$date)
for (i in 1:nrow(Score)) {
  Score$date[i]<-epoch(Score$date[i],"YMD")
}
#
##--- Conversion into numeric values -------------------------------------------
#
Score$Year<-as.numeric(Score$Year)
Score$Month<-as.numeric(Score$Month)
Score$Day<-as.numeric(Score$Day)
Score$date<-as.numeric(Score$date)
#
Score$Steps<-as.numeric(Score$Steps)
Score$Steps.B<-as.numeric(Score$Steps.B)
Score$reps<-as.numeric(Score$reps)
Score$Run<-as.numeric(Score$Run)
Score$Bike<-as.numeric(Score$Bike)
#
Score$Sys<-as.numeric(Score$Sys)
Score$Dia<-as.numeric(Score$Dia)
Score$HR<-as.numeric(Score$HR)
Score$Temp<-as.numeric(Score$Temp)
Score$SpO2<-as.numeric(Score$SpO2)
Score$glucose<-as.numeric(Score$glucose)
Score$ketones<-as.numeric(Score$ketones)
Score$ST<-as.numeric(Score$ST)
Score$ST.hour<-as.numeric(Score$ST.hour)
Score$Sys.1<-as.numeric(Score$Sys.1)
Score$Dia.1<-as.numeric(Score$Dia.1)
Score$HR.1<-as.numeric(Score$HR.1)
Score$Sys.2<-as.numeric(Score$Sys.2)
Score$Dia.2<-as.numeric(Score$Dia.2)
Score$HR.2<-as.numeric(Score$HR.2)
Score$Sys.3<-as.numeric(Score$Sys.3)
Score$Dia.3<-as.numeric(Score$Dia.3)
Score$HR.3<-as.numeric(Score$HR.3)
Score$Sys.4<-as.numeric(Score$Sys.4)
Score$Dia.4<-as.numeric(Score$Dia.4)
Score$HR.4<-as.numeric(Score$HR.4)
Score$Sys.5<-as.numeric(Score$Sys.5)
Score$Dia.5<-as.numeric(Score$Dia.5)
Score$HR.5<-as.numeric(Score$HR.5)
Score$Sys.6<-as.numeric(Score$Sys.6)
Score$Dia.6<-as.numeric(Score$Dia.6)
Score$HR.6<-as.numeric(Score$HR.6)
Score$Sys.7<-as.numeric(Score$Sys.7)
Score$Dia.7<-as.numeric(Score$Dia.7)
Score$HR.7<-as.numeric(Score$HR.7)
Score$Sys.8<-as.numeric(Score$Sys.8)
Score$Dia.8<-as.numeric(Score$Dia.8)
Score$HR.8<-as.numeric(Score$HR.8)
Score$Sys.9<-as.numeric(Score$Sys.9)
Score$Dia.9<-as.numeric(Score$Dia.9)
Score$HR.9<-as.numeric(Score$HR.9)
Score$Sys.10<-as.numeric(Score$Sys.10)
Score$Dia.10<-as.numeric(Score$Dia.10)
Score$HR.10<-as.numeric(Score$HR.10)
#
# Hormones
#
Score$Finasteride.Topic<-as.numeric(Score$Finasteride.Topic)
Score$Pregnenolone<-as.numeric(Score$Pregnenolone)
Score$Nasonex<-as.numeric(Score$Nasonex)
Score$DHEA<-as.numeric(Score$DHEA)
Score$Prednisone<-as.numeric(Score$Prednisone)
Score$Deflazacort<-as.numeric(Score$Deflazacort)
Score$Melatonin<-as.numeric(Score$Melatonin)
Score$T4<-as.numeric(Score$T4)
Score$T3<-as.numeric(Score$T3)
#
# Psicostimulants
#
Score$Bupropion<-as.numeric(Score$Bupropion)
Score$Provigil<-as.numeric(Score$Provigil)
Score$Sulbutiamine<-as.numeric(Score$Sulbutiamine)
Score$Aniracetam<-as.numeric(Score$Aniracetam)
Score$D_serine<-as.numeric(Score$D_serine)
Score$Sarcosine<-as.numeric(Score$Sarcosine)
#
# Antibiotics
#
Score$Amoxicillin<-as.numeric(Score$Amoxicillin)
Score$Ospen<-as.numeric(Score$Ospen)
Score$Penilevel<-as.numeric(Score$Penilevel)
Score$Cefuroxime<-as.numeric(Score$Cefuroxime)
#
# Liver
#
Score$Pioglitazone<-as.numeric(Score$Pioglitazone)
Score$Deursil<-as.numeric(Score$Deursil)
#
# Orthostatic intolleance/vascular
#
Score$Mestinon<-as.numeric(Score$Mestinon)
Score$Etilefrine<-as.numeric(Score$Etilefrine)
Score$Midodrine<-as.numeric(Score$Midodrine)
Score$Nimodipine<-as.numeric(Score$Nimodipine)
Score$Compression.suit<-as.numeric(Score$Compression.suit)
Score$Compression.socks<-as.numeric(Score$Compression.socks)
#
# PEM treatments (putative)
#
Score$Acetazolamide<-as.numeric(Score$Acetazolamide)
Score$Bismuth.subsalicylate<-as.numeric(Score$Bismuth.subsalicylate)
Score$Destrometorfano<-as.numeric(Score$Destrometorfano)
#
# Vitamines and minerals
#
Score$VitC<-as.numeric(Score$VitC)
Score$VitE<-as.numeric(Score$VitE)
Score$VitD3<-as.numeric(Score$VitD3)
Score$Se<-as.numeric(Score$Se)
Score$Omega3<-as.numeric(Score$Omega3)
Score$Methylene.blue<-as.numeric(Score$Methylene.blue)
Score$EPA<-as.numeric(Score$EPA)
Score$DHA<-as.numeric(Score$DHA)
#
# B vitamines
#
Score$B.complex<-as.numeric(Score$B.complex)
Score$Thiamin<-as.numeric(Score$Thiamin)
Score$Riboflavin<-as.numeric(Score$Riboflavin)
Score$Niacin<-as.numeric(Score$Niacin)
Score$VitB6<-as.numeric(Score$VitB6)
Score$Folic.acid<-as.numeric(Score$Folic.acid)
Score$Folinic.acid<-as.numeric(Score$Folinic.acid)
Score$Cyanocobalamin<-as.numeric(Score$Cyanocobalamin)
Score$HydroxycobalaminIM<-as.numeric(Score$HydroxycobalaminIM)
Score$HydroxycobalaminPO<-as.numeric(Score$HydroxycobalaminPO)
Score$Biotin<-as.numeric(Score$Biotin)
Score$Pantothenic.acid<-as.numeric(Score$Pantothenic.acid)
Score$Choline.bitartrate<-as.numeric(Score$Choline.bitartrate)
Score$Myo.Inositol<-as.numeric(Score$Myo.Inositol)
#
# Energy
#
Score$Creatine<-as.numeric(Score$Creatine)
Score$ALCar<-as.numeric(Score$ALCar)
#
# NSAIs
#
Score$Tachipirine<-as.numeric(Score$Tachipirine)
Score$Ibuprofene<-as.numeric(Score$Ibuprofene)
#
# Anti-H1 and anti-H2
#
Score$Ketotifen<-as.numeric(Score$Ketotifen)
Score$Desloratadina<-as.numeric(Score$Desloratadina)
Score$Famotidine<-as.numeric(Score$Famotidine)
#
# Stressors
#
Score$Assistant<-as.numeric(Score$Assistant)
Score$Cleaning<-as.numeric(Score$Cleaning)
Score$Relatives<-as.numeric(Score$Relatives)
Score$Travel<-as.numeric(Score$Travel)
Score$FLU<-as.numeric(Score$FLU)
#
# Exit variables
#
Score$Score_1<-as.numeric(Score$Score_1)
Score$Score_2<-as.numeric(Score$Score_2)
Score$Score_3<-as.numeric(Score$Score_3)
Score$Score_m<-as.numeric(Score$Score_m)
Score$PEM<-as.numeric(Score$PEM)
#
##--- Change of units ----------------------------------------------------------
#
Score$glucose<-Score$glucose*18 # from mmol/L to mg/dL
Score$ketones<-Score$ketones*18 # from mmol/L to mg/dL
#
##--- Calculate again Score_m and check for values out of scale ----------------
#
for (i in 1:nrow(Score)) {
  if (!is.na(Score$Score_1[i])&Score$Score_1[i]<4) Score$Score_1[i]<-4
  if (!is.na(Score$Score_2[i])&Score$Score_2[i]<4) Score$Score_2[i]<-4
  if (!is.na(Score$Score_3[i])&Score$Score_3[i]<4) Score$Score_3[i]<-4
  #
  if (!is.na(Score$Score_1[i])&Score$Score_1[i]>6) Score$Score_1[i]<-6
  if (!is.na(Score$Score_2[i])&Score$Score_2[i]>6) Score$Score_2[i]<-6
  if (!is.na(Score$Score_3[i])&Score$Score_3[i]>6) Score$Score_3[i]<-6
  #
  Score$Score_m[i]<-mean(c(Score$Score_1[i],Score$Score_2[i],Score$Score_3[i]),na.rm=T)
}
#
##--- Remove earlier measures, not used in this study --------------------------
#
# Score<-subset.data.frame(Score,Score$date>=as.integer(epoch("2021-03-01","YMD")))
#
##--- Check for missing days or inconsistencies in dates -----------------------
#
h<-1
lines<-c()
for (i in 2:nrow(Score)) {
  delta<-Score$date[i]-Score$date[i-1]
  if (delta!=60*60*24) {
    lines[h]<-i
    h<-h+1
  }
}
#
if (length(lines)>0) {
  stop("There is at least one date missing or duplicated! Check array lines")
}
#
##--- Create binary variable for WL --------------------------------------------
#
Score$WL.bin<-Score$WL
for (i in 1:nrow(Score)) {
  if (!is.na(Score$WL[i])) {
    if (Score$WL[i]!=0) Score$WL.bin[i]<-1
  }  
}
#
##--- Add equivalent drugs for class of drugs ----------------------------------
#
Score<-Beta.Lactam(Score) # an equivalent oral betalactam including (penicilline, amoxicilline, cftriaxone...)
Score<-Corticosteroid(Score) # includes prednison and deflaxacort
Score<-Hydroxycobalamin.Tot(Score) # includes oral and IM hydroxycobalamin
Score<-B12.Tot(Score) # includes oral and IM hydroxicobalamin and cyanocobalamin
#
##--- Create a score with the daily number of top scores (0, 1, 2, 3) ----------
#
Score$Score_T<-Score$Score_m
for (i in 1:nrow(Score)) {
  if (!is.na(Score$Score_1[i])&!is.na(Score$Score_2[i])&!is.na(Score$Score_3[i])) {
    count<-0
    if (Score$Score_1[i]==6) {
      count<-count+1
    }
    if (Score$Score_2[i]==6) {
      count<-count+1
    }
    if (Score$Score_3[i]==6) {
      count<-count+1
    }
    Score$Score_T[i]<-count  
  } else {
    Score$Score_T[i]<-NA
  }
}
#
#
#---- Estimate watch steps (Steps.B) from phone steps (Steps): Steps.C ---------
#
# Estimate missing values of watch steps (Steps.B) from phone steps (Steps) and 
# use the prediction to build a complete time series Steps.C
#
###-- Missing values and interpolation -----------------------------------------
#
# Count missing values to Steps
#
first.step<-which(!is.na(Score$Steps))[1]
missing.steps<-length(which(is.na(Score$Steps[first.step:nrow(Score)])))
missing.steps<-round(100*missing.steps/(nrow(Score)-first.step),1)
print(paste("Days with missing step count by cell phone:",missing.steps,"%"))
#
# Impute missing values to Steps if less than 10%
#
if (missing.steps<10) {
  Score$Steps<-zoo::na.approx(Score$Steps,maxgap=5,na.rm=FALSE)
}
#
# Steps.B should more than Steps 
#
index<-which(Score$Steps.B<(Score$Steps+50))
Score$Steps.B[index]<-NA
#
# Count missing values to Steps.B
#
first.step<-which(!is.na(Score$Steps.B))[1]
missing.steps<-length(which(is.na(Score$Steps.B[first.step:nrow(Score)])))
missing.steps<-round(100*missing.steps/(nrow(Score)-first.step),1)
print(paste("Days with missing step count by watch:",missing.steps,"%"))
#
# Impute missing values to Steps.B if less than 10%
#
if (missing.steps<10) {
  Score$Steps.B<-zoo::na.approx(Score$Steps.B,maxgap=5,na.rm=FALSE)
}
#
# Create Steps.C: watch steps imputed from Steps for days without Steps.B
#
df<-Steps_impute(Score)
Score$Steps.C<-df$Steps.C
#
# Plot time series
#
variables<-c("date","Steps","Steps.B","Steps.C")
start<-"2023-01-01"
end<-0
folder_path<-file.path(current_dir,"Output","Steps") 
type<-"Steps"
nd<-10
Data.plot(Score,nd,variables,start,end,folder_path,type,detrend=F,scaling=F) 
#
#---- Total daily mechanical work (joule) -------------------------------------
#
Score<-WORK(Score)
Score$Work<-Score$Work.a+Score$Work.an 
#
#
#---- Derive cardiac power index (CPI) -----------------------------------------
#
Score<-Cardiac.Derived(Score)
#
#---- Derive deltas from lean test, max and mean values ------------------------
#
LeanTest<-LeanTest_func(Score)
#
#---- RR series analysis -------------------------------------------------------
#
##--- Read HRV recordings ------------------------------------------------------
#
# acquire the names of the available HRV recordings
#
current_dir<-getwd()
sub_dir<-"Data/HRV"
folder_path<-file.path(current_dir,sub_dir)
file_names<-list.files(folder_path)
days<-length(file_names)
#
# store each recording in its RHRV list
#
hrv.data<-list()
for (i in 1:days) {
  hrv.data[[i]]<-CreateHRVData()
  hrv.data[[i]]<-SetVerbose(hrv.data[[i]],F)
  date.str<-gsub(".txt","",file_names[[i]]) # remove ".txt" from file name
  date.str<-strptime(date.str,format="%Y-%m-%d %H-%M-%S")
  date.str<-format(date.str,format="%d/%m/%Y %H:%M:%S")
  hrv.data[[i]]<-LoadBeatRR(hrv.data[[i]],file_names[[i]],
                            RecordPath=folder_path,scale=0.001,datetime=date.str)
}
#
##--- Time domain and frequency domain analysis --------------------------------
#
hrv.data<-HRV(hrv.data,days)
#
# remove recordings of less than 4 minutes and more than 6 minutes
#
index<-c(rep(NA,days))
for (i in 1:days) {
  length<-max(hrv.data[[i]]$Beat$Time)
  if (length<240) index[i]<-i
  if (length>360) index[i]<-i
}
index<-index[is.na(index)==F]
if (length(index)>0) hrv.data<-hrv.data[-index]
days<-length(hrv.data)
#
##--- Store relevant data in a data frame --------------------------------------
#
myHRVdata<-data.frame(date=c(rep(NA,days)),time=c(rep(NA,days)),avgHR=c(rep(NA,days)),
                      avgRR=c(rep(NA,days)),SDNN=c(rep(NA,days)),HRVi=c(rep(NA,days)),
                      pNN50=c(rep(NA,days)),rMSSD=c(rep(NA,days)),LF=c(rep(NA,days)),
                      HF=c(rep(NA,days)),LFonHF=c(rep(NA,days)),TotPower=c(rep(NA,days)),
                      LFnu=c(rep(NA,days)),HFnu=c(rep(NA,days)))
#
for (i in 1:days) {
  myHRVdata$date[i]<-as.integer(epoch(substr(hrv.data[[i]]$datetime,1,10),"YMD"))
  myHRVdata$time[i]<-substr(hrv.data[[i]]$datetime,12,19)
  myHRVdata$avgHR[i]<-mean(hrv.data[[i]]$Beat$niHR,na.rm=T)
  myHRVdata$avgRR[i]<-mean(hrv.data[[i]]$Beat$RR,na.rm=T)
  myHRVdata$SDNN[i]<-hrv.data[[i]]$TimeAnalysis[[1]]$SDNN
  myHRVdata$HRVi[i]<-hrv.data[[i]]$TimeAnalysis[[1]]$HRVi
  myHRVdata$pNN50[i]<-hrv.data[[i]]$TimeAnalysis[[1]]$pNN50
  myHRVdata$rMSSD[i]<-hrv.data[[i]]$TimeAnalysis[[1]]$rMSSD
  myHRVdata$LF[i]<-mean(hrv.data[[i]]$FreqAnalysis[[1]]$LF,na.rm=T)
  myHRVdata$HF[i]<-mean(hrv.data[[i]]$FreqAnalysis[[1]]$HF,na.rm=T)
  myHRVdata$LFonHF[i]<-myHRVdata$LF[i]/myHRVdata$HF[i]
  myHRVdata$TotPower[i]<-myHRVdata$LF[i]+myHRVdata$HF[i]
  myHRVdata$LFnu[i]<-myHRVdata$LF[i]/myHRVdata$TotPower[i]
  myHRVdata$HFnu[i]<-myHRVdata$HF[i]/myHRVdata$TotPower[i]
}
#
# remove eventual second recording for same day
#
myHRVdata<-myHRVdata[!duplicated(myHRVdata$date),]
#
#---- Data from BIA scale ------------------------------------------------------
#
BIA<-Fitdays()
#
# Imputation of data from BIA scale
#
first.day<-BIA$date[1]
last.day<-epoch(Sys.Date(),"YMD")
temp<-data.frame(date=seq(first.day,last.day,by=60*60*24))
BIA<-merge.data.frame(temp,BIA,by.y="date",all.x=T)
#
# add missing days with splines and means
#
for (j in 2:ncol(BIA)) {
  BIA[,j]<-approx(seq_along(BIA[,j]),BIA[,j],xout=seq_along(BIA[,j]))$y
  BIA[,j]<-Spline.Me(BIA[,j])
}
#
#---- Data from ShapeGanic App -------------------------------------------------
#
macro<-ShapeGanic()
#
#---- Data from Black Globe thermometer (Tenmars, TM-188D) ---------------------
#
BG<-Tenmars()
#
#---- Data from Bresser thermometer --------------------------------------------
#
Bresser.data<-Bresser()
#
#---- Data from Open Weather API -------------------
#
##--- Read saved data and retrieve missing days --------------------------------
#
for (type in c("weather","pollution")) {
  print(type)
  current_dir<-getwd()
  sub_dir<-"OpenWeather"
  folder_path<-file.path(current_dir,"Data",sub_dir)
  file.names<-list.files(folder_path)
  file.names<-file.path(folder_path,file.names)
  #
  file.name<-file.path(folder_path,paste0("OpenWeather_",type,".csv"))
  if (file.name %in% file.names) {
    df<-read.csv(file.name,sep=";",header=T)
    Score.todo<-anti_join(Score,df,by="date") # days we still have to retrieve environmental data for
  } else {
    df<-data.frame()
    Score.todo<-Score
  }
  #
  # retrieve
  #
  k<-1
  for (i in 1:nrow(Score.todo)) {
    Sys.sleep(0.1)
    timestamp<-as.integer(Score.todo$date[i])
    print(HRD(timestamp))
    city.name<-Score.todo$Place[i]
    print(city.name)
    if (is.na(city.name)) {
      next
    } else if (city.name=="Pordenone"|city.name=="Roma") {
      country.name<-"IT"
      time.zone<-"+01:00"
    } else if (city.name=="Rosario"|city.name=="Buenos Aires"|city.name=="Corrientes") {
      country.name<-"AR"
      time.zone<-"-03:00"
    } else if (city.name=="Asuncion") {
      country.name<-"PY"
      time.zone<-"-04:00"
    } else if (city.name=="Arona") {
      country.name<-"ES"
      time.zone<-"+00:00"
    } else {
      print(paste("You must indicate country and time zone of",city.name))
    }
    #
    k<-k+1 # select API key
    if(k>length(OWM_key)) k<-1
    api_key<-OWM_key[k]
    #
    limit<-5
    coordinates<-Get.coord(city.name,country.name,api_key,limit)
    lat<-coordinates$lat
    lon<-coordinates$lon
    place<-coordinates$name
    country<-coordinates$country
    #
    if (type=="weather") {
      a<-try(Get.weather(timestamp,time.zone,lat,lon,place,country,api_key),silent=T)
      if (is.data.frame(a)) df<-rbind(df,a)  
    } else if (type=="pollution") {
      a<-try(Get.pollution(timestamp,time.zone,lat,lon,place,country,api_key),silent=T)
      if (is.data.frame(a)) {
        index<-which(a[1,]<0)
        if (length(index)>0) a[1,index]<-NA 
        df<-rbind(df,a)
      }
    }
  }
  #
  if (type=="weather") {
    #
    OW.weather<-distinct(df)
    OW.weather<-OW.weather[order(OW.weather$date), ] # order by date  
    file.name<-file.path(folder_path,"OpenWeather_weather.csv")
    write.table(OW.weather,file=file.name,sep=";",col.names=T,row.names=F,quote=F)
    #
  } else if (type=="pollution") {
    #
    OW.pollution<-distinct(df)
    OW.pollution<-OW.pollution[order(OW.pollution$date), ] # order by date
    file.name<-file.path(folder_path,"OpenWeather_pollution.csv")
    write.table(OW.pollution,file=file.name,sep=";",col.names=T,row.names=F,quote=F)
    #
  }
}
#
##--- TAmean, absolute humidity and water vapor pressure -----------------------
#
# Calculate TAmean using TAmin and TAmax (there are mistakes in TAmean from API, because of missing data)
#
for (i in 1:nrow(OW.weather)) {
  if(!is.na(OW.weather$TAmin[i])&!is.na(OW.weather$TAmax[i])) {
    OW.weather$TAmean[i]<-mean(c(OW.weather$TAmin[i],OW.weather$TAmax[i]))  
  } else {
    OW.weather$TAmean[i]<-NA
  }
}
#
# Calculate absolute humidity (g/m^3) and water vapor pressure (hPa)
#
L<-nrow(OW.weather)
Abs.Humidity<-c()
ep<-c()
for (i in 1:L) {
  ep[i]<-(OW.weather$Humidity[i]/100)*6.112*exp(17.62*OW.weather$TAmean[i]/(243.12+OW.weather$TAmean[i]))
  value<-mu_v/(R*(273.15+OW.weather$TAmean[i]))
  Abs.Humidity[i]<-value*100*ep[i]
}
OW.weather$Abs.Humidity<-round(Abs.Humidity,2)
OW.weather$ep<-round(ep,2)
#
OW.weather<-subset.data.frame(OW.weather,select=-Place)
OW.weather<-subset.data.frame(OW.weather,select=-Country)
#
OW.pollution<-subset.data.frame(OW.pollution,select=-Place)
OW.pollution<-subset.data.frame(OW.pollution,select=-Country)
#
#---- Data from blood work -------------------------------------------------
#
BloodWork<-Blood.work() 
#
#---- Data from ActivityWatch ----------------------------------------------
#
ActivityWatch<-AWatch_func() 
#
#---- Data from Turbo Lamp + Shelly and calculation of wavelength at peak power ----
#
ShellyTurbo<-ShellyTurbo_func()
#
#---- Data from EEG readings ---------------------------------------------------
#
EEG<-EEG_func()
#
#---- Data from fNIRS readings -------------------------------------------------
#
year<-1979 # year of birth
fNIRS<-fNIRS_func(year)
#
#---- Merge all data -----------------------------------------------------------
#
mydata<-merge.data.frame(Score,myHRVdata,by.y="date",all.x=T)
mydata<-merge.data.frame(mydata,BIA,by.y="date",all.x=T)
mydata<-merge.data.frame(mydata,macro,by.y="date",all.x=T)
mydata<-merge.data.frame(mydata,BG,by.y="date",all.x=T)
mydata<-merge.data.frame(mydata,Bresser.data,by.y="date",all.x=T)
mydata<-merge.data.frame(mydata,OW.weather,by="date",all.x=T)
mydata<-merge.data.frame(mydata,OW.pollution,by="date",all.x=T)
mydata<-merge.data.frame(mydata,BloodWork,by="date",all.x=T)
mydata<-merge.data.frame(mydata,ActivityWatch,by="date",all.x=T)
mydata<-merge.data.frame(mydata,ShellyTurbo,by="date",all.x=T)
mydata<-merge.data.frame(mydata,EEG,by="date",all.x=T)
mydata<-merge.data.frame(mydata,fNIRS,by="date",all.x=T)
mydata<-merge.data.frame(mydata,LeanTest,by="date",all.x=T)
mydata<-mydata[order(mydata$date), ] # order by date
#
# Create a variable with total daily seconds of calls (cellular network + Whatsapp)
#
for (i in 1:nrow(mydata)) {
  mydata$Call.length.tot[i]<-mydata$WA.call[i]+mydata$Call.length[i]
}
#
# Check for missing days or inconsistencies in dates
#
h<-1
lines<-c()
for (i in 2:nrow(mydata)) {
  delta<-mydata$date[i]-mydata$date[i-1]
  if (delta!=60*60*24) {
    lines[h]<-i
    h<-h+1
  }
}
#
if (length(lines)>0) {
  stop("There is at leat one date missing or duplicated! Check array lines")
}
#
#---- Add binary variable for blood draw (zeros only) --------------------------
#
mydata$blood.draw[is.na(mydata$blood.draw)]<-0
#
#---- Save ---------------------------------------------------------------------
#
file.name<-"Output/all_daily_data.csv"
write.table(mydata,file=file.name,quote=F,row.names=F,col.names=T,sep=",")

