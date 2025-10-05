######################################
## EFI script to run Multi-model v1 ##
######################################
#--------------------------------------------
#Developed by Nicola Bozzolan (nicola.bozzolan@efi.int)
#             Sergey Zudin (sergey.zudin@efi.int)
#June 2025
#--------------------------------------------

#Upload required packages
#library(shiny)
#library(shinydashboard)
library(SoilR)
#library(shinyWidgets)
library(readxl)
library(dplyr)
library(scales)
#library(shinyvalidate)
library(ncdf4)
library(ncdf4.helpers)
library(PCICt)
library(lubridate)
#library(DT)
#library(leaflet)
library(R.utils)

#---------------------------
#Call required functions
#--------------------------
#---------------------------
#Forcing
source("Forcing_params_efi.R")
#---------------------------
#Models resolution and multi-model tool
source("Holisoils_multimodel_v1_i1.R")
source("Holisoils_RothC_v1_i1.R")
source("Holisoils_ICBM_v1_i1.R")
source("Holisoils_Century_v1_i1.R")
source("Holisoils_Yasso07_v1_i1.R")
source("Holisoils_Yasso20_v1_i1.R")
source("Holisoils_SG_v1.R")
#---------------------------
#--coded models
source("Yasso07Model_fixed.R")
#source("Millennial2Model.R")
source("Yasso20Model.R")
source("ICBM_environmental_functions.R")
#---------------------------
#Load databases
source("Upload_databases_efi.R")
#-----------------------------------------------------------------------
#--------Function to retrieve the C input from ISIMIP data-------------
#-----------------------------------------------------------------------
retreive_Cinput<-function(plot_lon,plot_lat,lon_rcp,lat_rcp,time_rcp,ncvar_rcp){
  
  #lon_rcp = lon_litter_rcp26
  #lat_rcp = lat_litter_rcp26
  #time_rcp = litter_time_rcp26
  #ncvar_rcp = litter_rcp26
  
  simulation_length<-min(plot_simulation_length,80)
  
  #Define lon-lat values to retrieve, based on user input
  lon_index_litter <- which.min(abs(lon_rcp - plot_lon))
  lat_index_litter <- which.min(abs(lat_rcp - plot_lat))
  #----------------------------------------------------
  #No spinup for (rcp26)
  #----------------------------------------------------
  #----------------------------------------------------
  #Retrieve forward (rcp26)
  #----------------------------------------------------
  #Convert time to dates
  litter_time_plot <- as.Date.character(format(time_rcp, "%Y-%m-%d"))
  print("!!!!!time plot!!!!!")
  #Select variable values from 2020 to 2099 for forward
  time_index_fwd_litter <- which(format(time_rcp, "%Y-%m-%d") >= "2020-01-01"
                                 & format(time_rcp, "%Y-%m-%d")< "2100-01-01")
  print("!!!!!time plot forward!!!!!")
  #For X timesteps
  litter_sel_lonlat_timestep_fwd = ncvar_rcp[lon_index_litter, lat_index_litter, time_index_fwd_litter]
  print("!!!!!time plot timestep!!!!!")
  #Convert variable from kg/m2/sec to t/ha/yr
  litter_sel_lonlat_timestep_fwd<-litter_sel_lonlat_timestep_fwd*(60*60*24*365.25)*10 #kg/m2/sec =>t/ha/yr
  print("!!!!!time plot time step forward!!!!!")
  #Select time series for fwd
  litter_time_sel_fwd = litter_time_plot[time_index_fwd_litter]
  litter_time_sel_plot_fwd = as.Date.character(format(litter_time_sel_fwd, "%Y-%m-%d"))
  print("!!!!!time plot %Y-%m-%d!!!!!")
  Litter_time_fwd = data.frame("Date"=litter_time_sel_plot_fwd,"Litter"=litter_sel_lonlat_timestep_fwd)
  print(length(Litter_time_fwd))
  print("!!!!!time plot litter time!!!!!")
  #Select intial date
  litter_sel_lonlat_timestep_fwd_sel = subset(Litter_time_fwd,Litter_time_fwd$Date>=start_date_simulations)
  max_length = min(simulation_length,length(Litter_time_fwd$Date))
  print("!!!!!time plot maxlenth!!!!!")
  #print("max_length")
  #print(max_length)
  #Select number of years
  litter_sel_lonlat_timestep_fwd_sel = litter_sel_lonlat_timestep_fwd_sel[1:max_length,]
  print("!!!!!time plot num years!!!!!")
  #Cinput_site=mean(litter_sel_lonlat_timestep_fwd_sel$Litter)
  Cinput_site=litter_sel_lonlat_timestep_fwd_sel$Litter
  #print(Cinput_site)
  
  
  return(Cinput_site)
}
#-----------------------------------------------------------------------
#--------Function to retrieve the Soil properties from Lucas data-------
#-----------------------------------------------------------------------
get_Soil<-function(plot_lon,plot_lat,lon_vector,lat_vector,var_matrix,ifna){
  lon_index_site <- which.min(abs(lon_vector - plot_lon))
  lat_index_site <- which.min(abs(lat_vector - plot_lat))
  retval <- var_matrix[lon_index_site,lat_index_site]
  if(is.na(retval)) retval <- ifna
  return(retval)
}
#--------------------------------------------------------------------------------
##This function retrieves ALL the climate data from ISIMIP for the selected site
#original version in Deliverable
#---------------------------------------------------------------------------------
retreive_clim_data_site<-function(plot_lon,plot_lat){
  ptm <- proc.time()
  
  
  #---------------------
  #Define forward simulation length [years]
  #as number of data/time_step
  #####################
  
  simulation_length<-min(plot_simulation_length,80)
  
  #####################
  
  #Choose how many plots to draw in the same figure (num of rows, num of columns)
  #par(mfrow=c(2,1),oma=c(2, 0, 0, 5))
  #---------------------
  ###########################
  #Read forcing climate data
  ###########################
  
  #ISIMIP data
  #----------------
  #Potential evapotransp data
  #----------------
  #Potevapotranspiration RCP26
  #----------------------------------------------------
  #Define lon-lat values to retrieve, based on user input
  lon_index_potevap_rcp26 <- which.min(abs(lon_potevap_rcp26 - plot_lon))
  lat_index_potevap_rcp26 <- which.min(abs(lat_potevap_rcp26 - plot_lat))
  #----------------------------------------------------
  #Retrieve spinup (rcp26)
  #----------------------------------------------------
  #Select variable values from 2006 to 2020 for spinup
  time_index_spinup_potevap_rcp26 <- which(format(potevap_time_rcp26, "%Y-%m-%d") >= "2006-01-01" &
                                             format(potevap_time_rcp26, "%Y-%m-%d") <= "2020-12-31")
  #For X timesteps
  potevap_sel_lonlat_timestep_rcp26_spinup = potevap_rcp26[lon_index_potevap_rcp26, lat_index_potevap_rcp26, time_index_spinup_potevap_rcp26]
  
  #Convert variable from kg m-2 s-1 to mm/month
  potevap_sel_lonlat_timestep_rcp26_spinup<-potevap_sel_lonlat_timestep_rcp26_spinup*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(potevap_sel_lonlat_timestep_rcp26_spinup)/12)
  
  #Select time series for spinup
  potevap_time_sel_spinup_rcp26 = potevap_time_plot_rcp26[time_index_spinup_potevap_rcp26]
  potevap_time_sel_plot_spinup_rcp26 = as.Date.character(format(potevap_time_sel_spinup_rcp26, "%Y-%m-01"))
  
  Potevap_month_spinup_rcp26 = data.frame("Date"=potevap_time_sel_plot_spinup_rcp26,"Potevap"=potevap_sel_lonlat_timestep_rcp26_spinup)
  
  #----------------------------------------------------
  #Retrieve forward (rcp26)
  #----------------------------------------------------
  #Select variable values from 2020 to 2099 for forward
  time_index_fwd_potevap_rcp26 <- which(format(potevap_time_rcp26, "%Y-%m-%d") >= "2020-01-01"
                                        & format(potevap_time_rcp26, "%Y-%m-%d")< "2100-01-01")
  #For X timesteps
  potevap_sel_lonlat_timestep_rcp26_fwd = potevap_rcp26[lon_index_potevap_rcp26, lat_index_potevap_rcp26, time_index_fwd_potevap_rcp26]
  
  #Convert variable from kg m-2 s-1 to mm/month
  potevap_sel_lonlat_timestep_rcp26_fwd<-potevap_sel_lonlat_timestep_rcp26_fwd*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(potevap_sel_lonlat_timestep_rcp26_fwd)/12)
  
  
  #Select time series for fwd
  potevap_time_sel_fwd_rcp26 = potevap_time_plot_rcp26[time_index_fwd_potevap_rcp26]
  potevap_time_sel_plot_fwd_rcp26 = as.Date.character(format(potevap_time_sel_fwd_rcp26, "%Y-%m-01"))
  
  Potevap_month_fwd_rcp26 = data.frame("Date"=potevap_time_sel_plot_fwd_rcp26,"Potevap"=potevap_sel_lonlat_timestep_rcp26_fwd)
  
  
  #Select intial date
  Potevap_month_fwd_rcp26_sel = subset(Potevap_month_fwd_rcp26,Potevap_month_fwd_rcp26$Date>=start_date_simulations)
  #Select number of years
  max_length = min(simulation_length*12,length(Potevap_month_fwd_rcp26_sel$Date))
  Potevap_month_fwd_rcp26_sel = Potevap_month_fwd_rcp26_sel[1:max_length,]
  
  # print("Potevap")
  # print(head(Potevap_month_fwd_rcp26_sel))
  # print(tail(Potevap_month_fwd_rcp26_sel))
  print("POTEVAP 26 OK")
  
  #----------------
  #Temperature data
  #----------------
  #Temperature RCP26
  
  #----------------------------------------------------
  #Define lon-lat values to retrieve, based on user input
  lon_index_temp_rcp26 <- which.min(abs(lon_temp_rcp26 - plot_lon))
  lat_index_temp_rcp26 <- which.min(abs(lat_temp_rcp26 - plot_lat))
  #----------------------------------------------------
  #Retrieve spinup (rcp26)
  #----------------------------------------------------
  print("??????temp index?????????")
  #Select variable values from 2006 to 2020 for spinup
  time_index_spinup_temp_rcp26 <- which(format(temp_time_plot_rcp26, "%Y-%m-%d") >= "2006-01-01" &
                                          format(temp_time_plot_rcp26, "%Y-%m-%d") <= "2020-12-31")
  print("??????temp spinup1?????????")
  #For X timesteps
  temp_sel_lonlat_timestep_rcp26_spinup = temp_rcp26[lon_index_temp_rcp26, lat_index_temp_rcp26, time_index_spinup_temp_rcp26]
  print("??????temp spinup2?????????")
  #Convert variable from K to C
  temp_sel_lonlat_timestep_rcp26_spinup<-temp_sel_lonlat_timestep_rcp26_spinup-273.15 #K to C
  print("??????temp convert?????????")
  #Select time series for spinup
  temp_time_sel_spinup_rcp26 = temp_time_plot_rcp26[time_index_spinup_temp_rcp26]
  temp_time_sel_plot_spinup_rcp26 = as.Date.character(format(temp_time_sel_spinup_rcp26, "%Y-%m-01"))
  print("??????temp almost spinup?????????")
  Temp_month_spinup_rcp26 = data.frame("Date"=temp_time_sel_plot_spinup_rcp26,"Temp"=temp_sel_lonlat_timestep_rcp26_spinup)
  print("??????temp spinup?????????")
  #----------------------------------------------------
  #Retrieve forward (rcp26)
  #----------------------------------------------------
  #Select variable values from 2020 to 2099 for forward
  time_index_fwd_temp_rcp26 <- which(format(temp_time_plot_rcp26, "%Y-%m-%d") >= "2020-01-01" 
                                     & format(temp_time_plot_rcp26, "%Y-%m-%d")< "2100-01-01")
  #For X timesteps
  temp_sel_lonlat_timestep_rcp26_fwd = temp_rcp26[lon_index_temp_rcp26, lat_index_temp_rcp26, time_index_fwd_temp_rcp26]
  
  #Convert variable from K to C
  temp_sel_lonlat_timestep_rcp26_fwd<-temp_sel_lonlat_timestep_rcp26_fwd-273.15 #K to C
  
  
  #Select time series for fwd
  temp_time_sel_fwd_rcp26 = temp_time_plot_rcp26[time_index_fwd_temp_rcp26]
  temp_time_sel_plot_fwd_rcp26 = as.Date.character(format(temp_time_sel_fwd_rcp26, "%Y-%m-01"))
  
  Temp_month_fwd_rcp26 = data.frame("Date"=temp_time_sel_plot_fwd_rcp26,"Temp"=temp_sel_lonlat_timestep_rcp26_fwd)
  
  
  #Select intial date
  Temp_month_fwd_rcp26_sel = subset(Temp_month_fwd_rcp26,Temp_month_fwd_rcp26$Date>=start_date_simulations)
  max_length = min(simulation_length*12,length(Temp_month_fwd_rcp26_sel$Date))
  #Select number of years
  Temp_month_fwd_rcp26_sel = Temp_month_fwd_rcp26_sel[1:max_length,]
  
  
  print("TEMP OK")
  
  #----------------
  #Precipitation data 
  #----------------
  #Precipitation RCP26
  #----------------------------------------------------
  #Define lon-lat values to retrieve, based on user input
  lon_index_prec_rcp26 <- which.min(abs(lon_prec_rcp26 - plot_lon))
  lat_index_prec_rcp26 <- which.min(abs(lat_prec_rcp26 - plot_lat))
  #----------------------------------------------------
  #Retrieve spinup (rcp26)
  #----------------------------------------------------
  #Select variable values from 2006 to 2020 for spinup
  time_index_spinup_prec_rcp26 <- which(format(prec_time_plot_rcp26, "%Y-%m-%d") >= "2006-01-01" &
                                          format(prec_time_plot_rcp26, "%Y-%m-%d") <= "2020-12-31")
  #For X timesteps
  prec_sel_lonlat_timestep_rcp26_spinup = prec_rcp26[lon_index_prec_rcp26, lat_index_prec_rcp26, time_index_spinup_prec_rcp26]
  
  #Convert variable from kg/m2/s (monthly mean) to mm/month
  prec_sel_lonlat_timestep_rcp26_spinup<-prec_sel_lonlat_timestep_rcp26_spinup*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(prec_sel_lonlat_timestep_rcp26_spinup)/12)
  
  #----------
  #Select time series for spinup
  prec_time_sel_spinup_rcp26 = prec_time_plot_rcp26[time_index_spinup_prec_rcp26]
  #prec_time_sel_plot_spinup_rcp26 = as.Date.character(format(prec_time_sel_spinup_rcp26, "%Y-%m-%d"))
  prec_time_sel_plot_spinup_rcp26 = as.Date.character(format(prec_time_sel_spinup_rcp26, "%Y-%m-01"))
  
  Prec_month_spinup_rcp26 = data.frame("Date"=prec_time_sel_plot_spinup_rcp26,"Precip"=prec_sel_lonlat_timestep_rcp26_spinup)
  print("??????prec spinup?????????")
  
  #----------------------------------------------------
  #Retrieve forward (rcp26)
  #----------------------------------------------------
  #Select variable values from 2020 to 2099 for forward
  time_index_fwd_prec_rcp26 <- which(format(prec_time_plot_rcp26, "%Y-%m-%d") >= "2020-01-01"
                                     & format(prec_time_plot_rcp26, "%Y-%m-%d")< "2100-01-01")
  #For X timesteps
  prec_sel_lonlat_timestep_rcp26_fwd = prec_rcp26[lon_index_prec_rcp26, lat_index_prec_rcp26, time_index_fwd_prec_rcp26]
  
  #Convert variable from kg/m2/s (monthly mean) to mm/month
  prec_sel_lonlat_timestep_rcp26_fwd<-prec_sel_lonlat_timestep_rcp26_fwd*(60*60*24)*rep(c(31,28,31,30,31,30,31,31,30,31,30,31),length(prec_sel_lonlat_timestep_rcp26_fwd)/12)
  
  
  #Select time series for fwd
  prec_time_sel_fwd_rcp26 = prec_time_plot_rcp26[time_index_fwd_prec_rcp26]
  prec_time_sel_plot_fwd_rcp26 = as.Date.character(format(prec_time_sel_fwd_rcp26, "%Y-%m-01"))
  
  Prec_month_fwd_rcp26 = data.frame("Date"=prec_time_sel_plot_fwd_rcp26,"Precip"=prec_sel_lonlat_timestep_rcp26_fwd)
  
  
  #Select intial date
  Prec_month_fwd_rcp26_sel = subset(Prec_month_fwd_rcp26,Prec_month_fwd_rcp26$Date>=start_date_simulations)
  max_length = min(simulation_length*12,length(Prec_month_fwd_rcp26_sel$Date))
  #Select number of years
  Prec_month_fwd_rcp26_sel = Prec_month_fwd_rcp26_sel[1:max_length,]
  
  
  print("PREC OK")
  
  #----------------
  #Soil moisture (top 18 cm)=>(m3/m3)
  #----------------
  
  #Soil moisture RCP26
  
  #----------------------------------------------------
  #Define lon-lat values to retrieve, based on user input
  lon_index_vswc_rcp26 <- which.min(abs(lon_vswc_rcp26 - plot_lon))
  lat_index_vswc_rcp26 <- which.min(abs(lat_vswc_rcp26 - plot_lat))
  #----------------------------------------------------
  #Retrieve spinup (rcp26)
  #----------------------------------------------------
  #Select variable values from 2006 to 2020 for spinup
  time_index_spinup_vswc_rcp26 <- which(format(vswc_time_rcp26, "%Y-%m-%d") >= "2006-01-01" &
                                          format(vswc_time_rcp26, "%Y-%m-%d") <= "2020-12-31")
  #For X timesteps
  vswc_sel_lonlat_timestep_rcp26_spinup = vswc_rcp26[lon_index_vswc_rcp26, lat_index_vswc_rcp26, time_index_spinup_vswc_rcp26]#mm3/mm3 top 18cm
  
  
  #Select time series for spinup
  vswc_time_sel_spinup_rcp26 = vswc_time_plot_rcp26[time_index_spinup_vswc_rcp26]
  vswc_time_sel_plot_spinup_rcp26 = as.Date.character(format(vswc_time_sel_spinup_rcp26, "%Y-%m-01"))
  
  Vswc_month_spinup_rcp26 = data.frame("Date"=vswc_time_sel_plot_spinup_rcp26,"Vswc"=vswc_sel_lonlat_timestep_rcp26_spinup)
  
  #----------------------------------------------------
  #Retrieve forward (rcp26)
  #----------------------------------------------------
  #Select variable values from 2020 to 2099 for forward
  time_index_fwd_vswc_rcp26 <- which(format(vswc_time_rcp26, "%Y-%m-%d") >= "2020-01-01"
                                     & format(vswc_time_rcp26, "%Y-%m-%d")< "2100-01-01")
  #For X timesteps
  vswc_sel_lonlat_timestep_rcp26_fwd = vswc_rcp26[lon_index_vswc_rcp26, lat_index_vswc_rcp26, time_index_fwd_vswc_rcp26] #mm3/mm3 top 18cm
  
  
  #Select time series for fwd
  vswc_time_sel_fwd_rcp26 = vswc_time_plot_rcp26[time_index_fwd_vswc_rcp26]
  vswc_time_sel_plot_fwd_rcp26 = as.Date.character(format(vswc_time_sel_fwd_rcp26, "%Y-%m-01"))
  
  Vswc_month_fwd_rcp26 = data.frame("Date"=vswc_time_sel_plot_fwd_rcp26,"Vswc"=vswc_sel_lonlat_timestep_rcp26_fwd)
  
  
  #Select intial date
  Vswc_month_fwd_rcp26_sel = subset(Vswc_month_fwd_rcp26,Vswc_month_fwd_rcp26$Date>=start_date_simulations)
  max_length = min(simulation_length*12,length(Vswc_month_fwd_rcp26_sel$Date))
  #Select number of years
  Vswc_month_fwd_rcp26_sel = Vswc_month_fwd_rcp26_sel[1:max_length,]
  
 
  
  
  print("SOIL MOIST OK")
  
  #====
  #Transform monthly to daily data
  #====
  #This function repeats the month's value for each day of the month
  #That way we have daily data, although with a monthly variability
  
  #Spinup
  Temp_day_spinup<-do.call("rbind", lapply(1:nrow(Temp_month_spinup_rcp26), function(i) 
    data.frame(Date = seq(Temp_month_spinup_rcp26$Date[i], 
                          (seq(Temp_month_spinup_rcp26$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
               Temp = Temp_month_spinup_rcp26$Temp[i])))
  
  # print("temp day ok")
  # print(tail(Temp_day_spinup))
  
  mean_days_in_months <- mean(c(31,28,31,30,31,30,31,31,30,31,30,31))
  
  Precip_day_spinup<-do.call("rbind", lapply(1:nrow(Prec_month_spinup_rcp26), function(i) 
    data.frame(Date = seq(Prec_month_spinup_rcp26$Date[i], 
                          (seq(Prec_month_spinup_rcp26$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
               Precip = Prec_month_spinup_rcp26$Precip[i]/mean_days_in_months)))
  
  Vswc_day_spinup<-do.call("rbind", lapply(1:nrow(Vswc_month_spinup_rcp26), function(i) 
    data.frame(Date = seq(Vswc_month_spinup_rcp26$Date[i], 
                          (seq(Vswc_month_spinup_rcp26$Date[i],length=2,by="months") - 1)[2], by = "1 days"), 
               Vswc = Vswc_month_spinup_rcp26$Vswc[i])))
  
  
  # print("prec day ok")
  # print(tail(Precip_day_spinup))
  
  #Fwd slected
  
  Temp_day_fwd_rcp26_sel<-do.call("rbind", lapply(1:nrow(Temp_month_fwd_rcp26_sel), function(i) 
    data.frame(Date = seq(Temp_month_fwd_rcp26_sel$Date[i]- day(Temp_month_fwd_rcp26_sel$Date[i]) +1,#from day 1
                          Temp_month_fwd_rcp26_sel$Date[i]- day(Temp_month_fwd_rcp26_sel$Date[i]) +as.numeric(days_in_month(Temp_month_fwd_rcp26_sel$Date[i])),#to last day month
                          by = "1 days"), 
               Temp = Temp_month_fwd_rcp26_sel$Temp[i])))
  # print("Temp")
  # print(head(Temp_day_fwd_rcp26_sel))
  # print(tail(Temp_day_fwd_rcp26_sel))
  
  Precip_day_fwd_rcp26_sel<-do.call("rbind", lapply(1:nrow(Prec_month_fwd_rcp26_sel), function(i) 
    data.frame(Date = seq(Prec_month_fwd_rcp26_sel$Date[i]- day(Prec_month_fwd_rcp26_sel$Date[i]) +1,#from day 1
                          Prec_month_fwd_rcp26_sel$Date[i]- day(Prec_month_fwd_rcp26_sel$Date[i]) +as.numeric(days_in_month(Prec_month_fwd_rcp26_sel$Date[i])),#to last day month
                          by = "1 days"), 
               Precip = Prec_month_fwd_rcp26_sel$Precip[i]/mean_days_in_months)))
  
  Vswc_day_fwd_rcp26_sel<-do.call("rbind", lapply(1:nrow(Vswc_month_fwd_rcp26_sel), function(i) 
    data.frame(Date = seq(Vswc_month_fwd_rcp26_sel$Date[i]- day(Vswc_month_fwd_rcp26_sel$Date[i]) +1,#from day 1
                          Vswc_month_fwd_rcp26_sel$Date[i]- day(Vswc_month_fwd_rcp26_sel$Date[i]) +as.numeric(days_in_month(Vswc_month_fwd_rcp26_sel$Date[i])),#to last day month
                          by = "1 days"), 
               Vswc = Vswc_month_fwd_rcp26_sel$Vswc[i])))
  
  
  #Fwd climate change scenario RCP26
  
  Temp_day_fwd_rcp26<-do.call("rbind", lapply(1:nrow(Temp_month_fwd_rcp26), function(i)
    data.frame(Date = seq(Temp_month_fwd_rcp26$Date[i]- day(Temp_month_fwd_rcp26$Date[i]) +1,#from day 1
                          Temp_month_fwd_rcp26$Date[i]- day(Temp_month_fwd_rcp26$Date[i]) +as.numeric(days_in_month(Temp_month_fwd_rcp26$Date[i])),#to last day month
                          by = "1 days"),
               Temp = Temp_month_fwd_rcp26$Temp[i])))
  
  
  # print("temp day ")
  # print(tail(Temp_day_fwd_rcp26,20))
  
  Precip_day_fwd_rcp26<-do.call("rbind", lapply(1:nrow(Prec_month_fwd_rcp26), function(i) 
    data.frame(Date = seq(Prec_month_fwd_rcp26$Date[i]- day(Prec_month_fwd_rcp26$Date[i]) +1,#from day 1
                          Prec_month_fwd_rcp26$Date[i]- day(Prec_month_fwd_rcp26$Date[i]) +as.numeric(days_in_month(Prec_month_fwd_rcp26$Date[i])),#to last day month
                          by = "1 days"), 
               Precip = Prec_month_fwd_rcp26$Precip[i]/mean_days_in_months)))
  # print("precip day ")
  # print(tail(Precip_day_fwd_rcp26,20))
  
  Vswc_day_fwd_rcp26<-do.call("rbind", lapply(1:nrow(Vswc_month_fwd_rcp26), function(i) 
    data.frame(Date = seq(Vswc_month_fwd_rcp26$Date[i]- day(Vswc_month_fwd_rcp26$Date[i]) +1,#from day 1
                          Vswc_month_fwd_rcp26$Date[i]- day(Vswc_month_fwd_rcp26$Date[i]) +as.numeric(days_in_month(Vswc_month_fwd_rcp26$Date[i])),#to last day month
                          by = "1 days"), 
               Vswc = Vswc_month_fwd_rcp26$Vswc[i])))
  
  
  
  print("data read ok")
  
  print("Elapsed time (minutes)")
  print((proc.time() - ptm)/60) #minutes
  
  
  ###############################
  #Return a list of 20 elements
  ##############################
  #spinup - fwd rcp2.6 - fwd rcp 6.0
  #1-3: potential evapotranspiration (mm/month) - sel
  #4-6: precipitation (mm), daily -sel
  #7-9: temperature (C), daily -sel
  #10-12: soil moisture (mm3/mm3) top 10cm, daily -sel
  
  #13-14: potential evapotranspiration (mm/month) - all
  #15-16: precipitation (mm), daily -all
  #17-18: temperature (C), daily -all
  #19-20: soil moisture (mm3/mm3) top 10cm, daily -all
  
  return(list(Potevap_month_spinup_rcp26,Potevap_month_fwd_rcp26_sel,Potevap_month_fwd_rcp26_sel,#Potevap_month_fwd_rcp60_sel,
              Precip_day_spinup,Precip_day_fwd_rcp26_sel,Precip_day_fwd_rcp26_sel,#Precip_day_fwd_rcp60_sel,
              Temp_day_spinup,Temp_day_fwd_rcp26_sel,Temp_day_fwd_rcp26_sel,#Temp_day_fwd_rcp60_sel,
              Vswc_day_spinup,Vswc_day_fwd_rcp26_sel,Vswc_day_fwd_rcp26_sel,#Vswc_day_fwd_rcp60_sel,
              Potevap_month_fwd_rcp26,Potevap_month_fwd_rcp26,#Potevap_month_fwd_rcp60,
              Precip_day_fwd_rcp26,Precip_day_fwd_rcp26,#Precip_day_fwd_rcp60,
              Temp_day_fwd_rcp26,Temp_day_fwd_rcp26,#Temp_day_fwd_rcp60,
              Vswc_day_fwd_rcp26,Vswc_day_fwd_rcp26))#Vswc_day_fwd_rcp60))
}
#------------------------------------------------------------
#--------Function to calculate C input from mortality event-------------
#-----------------------------------------------------------
calculate_Cin_mortality_event<-function(Cinput_ag,Cinput_bg,Cbiomass_cveg_ag,Cbiomass_cveg_bg,mortality_index, harvest_index,year_disturbance){
  
  #===MORTALITY EVENT====
  #print("Retrieving C input for disturbance scenario")
  ###Supposing mortality event at year 1
  mort_year = year_disturbance
  ##########################################
  #Carbon input from survival vegetation
  ##########################################
  #AGin (C input from aboveground survival vegetation) after mortality event
  #retreive_Cinput_ag_rcp26_fix_NOTDIED = Cinput_ag[c(mort_year:length(Cinput_ag))]*(mortality_index/100.) #Amount of abouveground litter input that remains after mortality
  retreive_Cinput_ag_rcp26_fix_NOTDIED = Cinput_ag[c(mort_year:length(Cinput_ag))]*(1-mortality_index/100.) #Amount of abouveground litter input that remains after mortality
  #BGin_notdied (C input from belowground survival vegetation) after mortality event
  retreive_Cinput_bg_rcp26_fix_NOTDIED = 0.333*(1.92*(100.*retreive_Cinput_ag_rcp26_fix_NOTDIED)+130.)*0.01 #belowground C input of the NON died vegetation #where *100 is the conversion factor from tC/ha to (g/m2)
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
  #>
  
  #AG_cveg
  Cbiomass_rcp26_ag = Cbiomass_cveg_ag
  #AGdied
  #Cbiomass_rcp26_ag_DIED = Cbiomass_rcp26_ag*(input$Mort_rate/100.) #Amount of aboveground biomass died after disturbance event
  Cbiomass_rcp26_ag_DIED = Cbiomass_rcp26_ag*(mortality_index/100.) #Amount of aboveground biomass died after disturbance event
  #AGkept
  #Cbiomass_rcp26_ag_KEPT = Cbiomass_rcp26_ag_DIED*(1-input$Harv_rate/100.) #Amount of aboveground biomass kept on the soil after harvest
  Cbiomass_rcp26_ag_KEPT = Cbiomass_rcp26_ag_DIED*(1-harvest_index/100.) #Amount of aboveground biomass kept on the soil after harvest
  
  #FOL_TOTbiomass
  Cbiomas_rcp26_ag_FOLIAR = (Cbiomass_cveg_ag+Cbiomass_cveg_bg)*0.15 #foliar biomass as a fraction of TOT biomass, see Konôpka et al., 2021 for different tree species
  #FOLdied
  #Cbiomas_rcp26_ag_DIED_FOLIAR = Cbiomas_rcp26_ag_FOLIAR*(input$Mort_rate/100.) #foliar biomass of the died vegetation
  Cbiomas_rcp26_ag_DIED_FOLIAR = Cbiomas_rcp26_ag_FOLIAR*(mortality_index/100.) #foliar biomass of the died vegetation
  Cbiomass_rcp26_bg_DIED = 0.333*(1.92*(100.*Cbiomas_rcp26_ag_DIED_FOLIAR)+130.)*0.01 #belowground biomass of the died vegetation #where *100 is the conversion factor from tC/ha to (g/m2)
  
  
  #if(input$Mort_rate>0){
  if(mortality_index>0){
    #AGinTOT (time series)
    Cinput_disturbance_ag = Cinput_ag
    #After disturbance
    Cinput_disturbance_ag[c(mort_year:length(Cinput_disturbance_ag))] = retreive_Cinput_ag_rcp26_fix_NOTDIED
    #Add additional AG biomass at year of disturbance
    Cinput_disturbance_ag[mort_year] = Cinput_disturbance_ag[mort_year]+Cbiomass_rcp26_ag_KEPT
    
    #BGinTOT (time series)
    Cinput_disturbance_bg = Cinput_bg
    
    #After disturbance
    Cinput_disturbance_bg[c(mort_year:length(Cinput_disturbance_bg))] =  retreive_Cinput_bg_rcp26_fix_NOTDIED
    #Add additional BG biomass at year of disturbance
    Cinput_disturbance_bg[mort_year] = Cinput_disturbance_bg[mort_year]+Cbiomass_rcp26_bg_DIED
    
  }else{
    #AGinTOT (time series)
    Cinput_disturbance_ag = Cinput_ag
    #BGinTOT (time series)
    Cinput_disturbance_bg = Cinput_bg
  }
  
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  #Cin TOT in case of disturbance (time series)
  #Cinput_disturbance = Cinput_disturbance_ag+Cinput_disturbance_bg
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  return(list(Cinput_disturbance_ag,Cinput_disturbance_bg))
}

plot_simulation_length<-20
computation_time_step_fwd<-1
start_date_simulations<-as.Date("2021-01-01")
#read input from file; filename fixed in this version:input/plots_input.csv 
infile<-paste("input/plots_input.csv",sep="")
dfin<-read.csv(infile)
for (i in rownames(dfin)) {
  plot_name<-dfin[i,"plotname"]
  lon<-dfin[i,"lon"]
  lat<-dfin[i,"lat"]
  dom_species<-dfin[i,"dom_species"]
  step<-dfin[i,"step"]
  AWEN_input<-list_AWEN[dom_species]
  DR_in<-(AWEN_input[[1]][[2]]+AWEN_input[[1]][[3]])/(AWEN_input[[1]][[1]]+AWEN_input[[1]][[4]])
  plotBiomass_ag<-dfin[i,"Biomass_ag"]
  plotBiomass_bg<-dfin[i,"Biomass_bg"]
  #Year of disturbace/harvest
  year_dist<-1
  #Soil parameters
  soil_plot<-get_Soil(lon,lat,lon_SOC,lat_SOC,SOC_data,49)
  soil_plot
  clay_plot<-get_Soil(lon,lat,lon_clay,lat_clay,clay_data,49)
  silt_plot<-get_Soil(lon,lat,lon_silt,lat_silt,silt_data,49)
  CN_plot<-get_Soil(lon,lat,lon_CN,lat_CN,CN_data,10)
  BD_plot<-get_Soil(lon,lat,lon_BD,lat_BD,BD_data,1.5)
  
  retreive_Cinput_ag_rcp26_fix = retreive_Cinput(lon,lat,lon_litter_rcp26,lat_litter_rcp26,litter_time_rcp26,litter_rcp26_ag)
  retreive_Cinput_bg_rcp26_fix = retreive_Cinput(lon,lat,lon_litter_rcp26,lat_litter_rcp26,litter_time_rcp26,litter_rcp26_bg)
  plot_clim_data<-retreive_clim_data_site(lon,lat)
  res_columns=c("harvest_ratio","mortality","year1","year2","year3","year4","year5","year6","year7","year8","year9","year10","year11","year12","year13","year14","year15","year16","year17","year18","year19","year20")
  res_data<-data.frame(matrix(nrow=0,ncol=length(res_columns)))
  colnames(res_data) = res_columns
  start_time <- Sys.time()
  #processing
  for (harvest_ratio in seq(0,90,step)) {
    for (mortality_ratio in seq(0,90,step)) { 
  
      Cinput_disturbance_vec<-calculate_Cin_mortality_event(retreive_Cinput_ag_rcp26_fix,retreive_Cinput_bg_rcp26_fix,plotBiomass_ag,plotBiomass_bg,mortality_ratio,harvest_ratio,year_dist)
      #Call the multimodel with comment of params
      test_efi <-Call_MULTIMODEL_i1(plot_figures='false',#plot_figures, #set to true in Forcing_params.R
                                    simulation_length=plot_simulation_length, #user defined
                                    spinup_length=spinup_length, #set to 5000 in Forsing_params.R
                                    computation_time_step_fwd=computation_time_step_fwd, #user defined set to 1 in app.R
                                    start_date_simulations=start_date_simulations, #User defined
                                    temperature_spinup=plot_clim_data[[7]],#Temp_day_spinup, #Temp_day_spinup=data_clim[[7]]; data_clim - returns by retreive_clim_data_site()
                                    precipitation_spinup=plot_clim_data[[4]],#Precip_day_spinup,#data_clim[[4]]
                                    potential_evapotranspiration_spinup=plot_clim_data[[1]],#Potevap_month_spinup_rcp26,#data_clim[[1]]
                                    soilmoisture_spinup=plot_clim_data[[10]]$Vswc,#as.numeric(Vswc_day_spinup$Vswc),#data_clim[[10]]
                                    temperature_fwd=plot_clim_data[[8]],#Temp_day_fwd_rcp26_sel,#data_clim[[8]]
                                    precipitation_fwd=plot_clim_data[[5]],#Precip_day_fwd_rcp26_sel,#data_clim[[5]]
                                    potential_evapotranspiration_fwd=plot_clim_data[[2]],#Potevap_month_fwd_rcp26_sel,#data_clim[[2]]
                                    soilmoisture_fwd=plot_clim_data[[11]]$Vswc,#as.numeric(Vswc_day_fwd_rcp26_sel$Vswc),#data_clim[[11]]
                                    SOC_0=soil_plot,#input$SOC,#get_soil(,,,,SOC_data,)
                                    C_input_ag_spinup=as.numeric(mean(retreive_Cinput_ag_rcp26_fix)),#retreive_Cinput<-function(lon_rcp,lat_rcp,time_rcp,ncvar_rcp)
                                    C_input_bg_spinup=as.numeric(mean(retreive_Cinput_bg_rcp26_fix)),#retreive_Cinput<-function(lon_rcp,lat_rcp,time_rcp,ncvar_rcp)
                                    C_input_ag_fwd=as.numeric(Cinput_disturbance_vec[[1]]),#as.numeric(retreive_Cinput_ag_rcp26_fix),#uses Landuse
                                    C_input_bg_fwd=as.numeric(Cinput_disturbance_vec[[2]]),#as.numeric(retreive_Cinput_bg_rcp26_fix),#uses Landuse
                                    clay_p=clay_plot,#input$clay_slider,#get_soil(,,,,clay_data,)
                                    silt_p=silt_plot,#input$silt,#get_soil(,,,,silt_data,49)
                                    soil_thickness=soil_OM_thick,#input$soilthick,#from user input; default 25 in Forcing_params.R
                                    pH_p=ph_site,#set to 6 in Forcing_params.R
                                    lignin_to_nitrogen=lignin_to_nitrogen_ratio,#input$LNratio,#lignin_to_nitrogen_ratio=0.5 set in Forcing_params.R
                                    structural_in_lignin=structural_in_lignin_ratio,#input$SLratio,#structural_in_lignin_ratio=0.1 set in Forcing_params.R
                                    woodylittersize=woodylittersize_scalar,#input$WLS,#set to 2 in Forsing_params.R
                                    AWEN_in=AWEN_input[[1]],#df_AWEN from upload_databases.R by species
                                    decomp_to_resist_ratio=DR_in,
                                    CN_Ratio=CN_plot,#input$CNratio,#get_soil(,,,,CN_data,)
                                    Bulk_Density=BD_plot,#input$bulkdensity,#get_soil(,,,,BD_data,)
                                    WFPS=water_filled_pore_space,#set to -1 in Forcing_params.R
                                    CH4_Conc=1777,#input$CH4_data,#df_CH4 by year
                                    decomposition_param_RothC=ksRothC,#set in Forcing_params.R
                                    decomposition_param_ICBM=param_ICBM,#set in Forcing_params.R
                                    decomposition_param_Century=ksCent,#set in Forcing_params.R
                                    decomposition_param_Yasso07=paramYasso07,#set in Forcing_params.R
                                    decomposition_param_Yasso20=paramYasso20)#set in Forcing_params.R
      
      print("End of multimodel simulation RCP2.6")
      res_data<-rbind(res_data,c(harvest_ratio,mortality_ratio,test_efi[[7]]))
      colnames(res_data) = res_columns
      
    }
  }
  colnames(res_data) = res_columns
  outfile<-paste("out/",plot_name,".csv",sep="")
  write.csv(res_data,outfile)
  end_time <- Sys.time()
  extime<-end_time - start_time
  print("executied in:")
  print(extime)
}