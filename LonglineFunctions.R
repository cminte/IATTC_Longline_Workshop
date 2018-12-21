# ######################################################
# Load libraries and functions for
# Data preparation and mapping
# C.M. Minte-Vera IATTC cminte@iattc.org
# Version Dec.19.2018 # This version is almost final 
# Need to have coastline and EEZ shape files to use it
# Files: "ne_110m_coastline.shp"
# "eez_boundaries_v10.shp"
# Need to be in the correct folders
# #put here the directories where the shp files are 
EEZ.path<-"C:/Users/cminte/Documents/Carolina2018/Meetings/LONGLINE workshop/@Working Code/World_EEZ_v10_20180221"
coast.path<-"C:/Users/cminte/Documents/Carolina2018/Meetings/LONGLINE workshop/@Working Code/ne_110m_coastline"
#  
##############################################

require(tidyverse)
require("RColorBrewer")
require(rgdal) #Read OGR vector maps into Spatial objects
require(tools) #Utilities for listing files, and manipulating file paths.

#######################FUNCTIONS
fortify.shape <- function(x){
  x@data$id <- rownames(x@data)
  x.f <- fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
}
setup_PO_regions<-function(dat,reg_configuration=1)
{
  require(tidyverse)
  dat$region<- NA
  if(reg_configuration==1) #Initial configuration based on the Tree analsyis
  {
    dat<-mutate(dat,region=replace(region, which(dat$lat5>10 & dat$lon5< I(-110) & dat$lon5>-150),"EPO0")) %>%
      mutate(region=replace(region, which(dat$lat5> -10 & dat$lat5<10 & dat$lon5< I(-110) & dat$lon5>-150),"EPO1")) %>%
      mutate(region=replace(region, which(dat$lat5> -10 & dat$lon5>-110),"EPO2")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-10) & dat$lon5< I(-110) & dat$lon5>-150),"EPO3")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-10) & dat$lon5>-110),"EPO4")) %>%
      mutate(region=replace(region, which(dat$lat5>10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO2")) %>%
      mutate(region=replace(region, which(dat$lat5<10 &dat$lat5>-10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO4")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-10) & dat$lon5< I(-150) & dat$lon5>-180),"WCPO6"))
    dat$region<- as.factor(dat$region)
  }
  if(reg_configuration==2) #split EPO at 0 west of 110W, isolate the Galapagos area
  {
    dat<-mutate(dat,region=replace(region, which(dat$lat5>5 & dat$lon5> I(-150)),"EPO1")) %>%
      mutate(region=replace(region, which(dat$lat5<5 & dat$lat5> 0 & dat$lon5< I(-110) & dat$lon5> -150),"EPO1")) %>% 
      mutate(region=replace(region, which(dat$lat5<0 & dat$lon5< I(-110) & dat$lon5>-150), "EPO2")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-5) & dat$lon5> -110),"EPO2")) %>%                   
      mutate(region=replace(region, which(dat$lat5<5 & dat$lat5> -5 & dat$lon5>-110),"Galapagos")) %>%
      mutate(region=replace(region, which(dat$lat5>10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO2")) %>%
      mutate(region=replace(region, which(dat$lat5<10 & dat$lat5>-10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO4")) %>%
      mutate(region=replace(region, which(dat$lat5<I(-10) & dat$lon5< I(-150) & dat$lon5>-180),"WCPO6"))
    dat$region<- as.factor(dat$region)
  }
  if(reg_configuration==3) #split EPO at 0 west of 110W, isolate east of 110 W
  {
    dat<-mutate(dat,region=replace(region, which(dat$lat5>5 & dat$lon5> I(-150)),"EPO1")) %>%
      mutate(region=replace(region, which(dat$lat5> 0 & dat$lon5< I(-110) & dat$lon5>-150),"EPO1")) %>%
      mutate(region=replace(region, which(dat$lat5< 0 & dat$lon5< I(-110)& dat$lon5>-150),"EPO2")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-5) & dat$lon5> -110),"EPO2")) %>%
      mutate(region=replace(region, which(dat$lat5< 5 & dat$lat5& dat$lon5>-110),"East110W")) %>%
      mutate(region=replace(region, which(dat$lat5> -5 & dat$lat5& dat$lon5>-110),"East110W")) %>%
      mutate(region=replace(region, which(dat$lat5> 10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO2")) %>%
      mutate(region=replace(region, which(dat$lat5<10 &dat$lat5>-10 & dat$lon5< I(-150) & dat$lon5>-180),"WCPO4")) %>%
      mutate(region=replace(region, which(dat$lat5< I(-10) & dat$lon5< I(-150) & dat$lon5>-180),"WCPO6"))
    dat$region<- as.factor(dat$region)
  }
  
  if(reg_configuration==4) #split EPO and WCPO
  {
    dat<-mutate(dat,region=replace(region, which(dat$lon5> -150),"EPO")) %>%
      mutate(region=replace(region, which(dat$lon5 < I(-150)),"WCPO")) 
    dat$region<- as.factor(dat$region)
  }
  
  if(reg_configuration==5) #split jurisdiction (EPO and WCPO)
  {
    dat<-mutate(dat,region=replace(region, which(dat$lon5> I(-150)),"EPO - IATTC")) %>%
      mutate(region=replace(region, which(dat$lon5 < I(-130) & dat$lon5> I(-150) & dat$lat5 < I(-4)),"EPO - Overlap")) %>%
      mutate(region=replace(region, which(dat$lon5< I(-150)),"WCPO - WCPFC")) 
    dat$region<- as.factor(dat$region)
  }
  
  return(dat)
}


#numbers of hooks per set by 5yr interval and region
setup_PO_time_blocks<-function(dat,time_step=5)
{
  #need year variable as op_yr
  
  dat$time_block<- NA
  min_year<-min(dat$op_yr,na.rm=T)
  
  if(time_step==5) # 5 yrs time blocks
  { 
    Mylabels=c("1975-79","1980-84","1985-89","1990-94","1995-99","2000-04","2005-09","2010-14","2015-17")
  }
  dat<- mutate(dat,time_block=(1 + floor((op_yr-min_year)/time_step))) %>%
    mutate(time_block=parse_factor(time_block,levels=unique(sort(time_block))))
  dat$time_block=factor(dat$time_block,levels=levels(dat$time_block),labels=Mylabels)
  return(dat)
}

## Set Ut HBF categories
#hbf 
#Sung Il:
#<=8 (8> <=15), (15>,<=20 ) >20
#JPN figure:
#<=5 , 5-7, 8-10, 11-13, 14-16, 17-19, 20-21, >21

#Hooks between floats * Number of Hooks: EPO and WCPO and  EPO1 and WCPO4 -  by categories
#Same but as % - with Sung Il categories and Cleridy's colors
setup_hbf_categories<-function(dat,configuration=1)
{
  #need hooks between floats variable as hbf
  dat$hbf_cate<- NA
  if(configuration==1) # JPN figure
  { 
    Mylabels=c("<5","5-7","8-10","11-13","14-16","17-19","20-21",">21")  
    dat<- mutate(dat,hbf_cate=replace(hbf_cate, which(hbf<5),Mylabels[1])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>4 & hbf <8), Mylabels[2])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>7 & hbf <11),Mylabels[3])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>10 & hbf <14),Mylabels[4])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>13 & hbf <17),Mylabels[5])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>16 & hbf <20),Mylabels[6])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>19 & hbf <22),Mylabels[7])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>21),Mylabels[8])) 
    dat$hbf_cate<- as.factor(dat$hbf_cate)
    
  }
  if(configuration==2) # Sung Il
  { 
    #<=8 (8> <=15), (15>,<=20 ) >20
    Mylabels=c("<9","9-15","16-20",">20")  
    dat<- mutate(dat,hbf_cate=replace(hbf_cate, which(hbf<9),Mylabels[1])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>8 & hbf <16), Mylabels[2])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>15 & hbf <21),Mylabels[3])) %>%
      mutate(hbf_cate=replace(hbf_cate, which(hbf>20),Mylabels[4])) 
    dat$hbf_cate<- as.factor(dat$hbf_cate)
    
  }
  dat<-mutate(dat, hbf_cate=parse_factor(hbf_cate,levels=Mylabels))
  return(dat)
}
##########FUNCTION FOR MAPS#############
#with tiles, EEZ and continents (12.19.2018)
Do_Map<-function(lat,lon,variable,by,label,MyScale,cols,Mylegend,MyTitle="Korea",
                    MyCaption="Data source: Operational Level Data",save=TRUE){
  
  MyDat<-as.data.frame(cbind(lat,lon,variable,by))
  p1<- ggplot()+
    geom_tile(data=MyDat,mapping=aes(x=lon,y=lat,fill=cut(variable,breaks=MyScale))) +
    coord_quickmap(ylim = c(-50,50),xlim = c(-250,-70))+
    scale_fill_manual(values = cols,labels=levels(cut(variable,breaks=MyScale)))+  
    geom_path(data = EEZ, aes(x = long, y = lat, group = group), 
              color = "#606060", size = 0.25) +
    geom_path(data = dat.coast, aes(x = long, y = lat, group = group), 
              color = "black", size = 0.5) +
    #overlap WCPFC and IATTC
    geom_polygon(data=as.data.frame(cbind("long"=c(-150,-150,-130,-130),"lat"=c(-50,-4,-4,-50))),
                 mapping=aes(x=long,y=lat),fill="cyan",alpha=0.1,color="black",linetype = 2,size=0.5)+
    scale_x_continuous(limits = c(-250,-70), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-50,50), expand = c(0, 0)) +
    labs( title =MyTitle,
          subtitle = str_c("Map of",label),
          caption = MyCaption  ) +
    guides(fill=guide_legend(title=Mylegend)) +
    labs(list( x = "Longitude", y = "Latitude")) +
    theme_bw(8)+
    facet_wrap(~by)
  
  if(save) ggsave(str_c("Figure_Maps_",label,".pdf"),plot=p1,scale=2)
  else print(p1)
}

####################### END FUNCTIONS 

## COASTLINE SHAPE file
print("Attention: **Need to download Coastline shape from https://www.naturalearthdata.com/downloads/110m-physical-vectors/ **")
#https://www.naturalearthdata.com/downloads/110m-physical-vectors/
file.coast <- "ne_110m_coastline.shp"
dat.coast <- readOGR(dsn = coast.path,layer = file_path_sans_ext(file.coast))
#Fortify the shapefile data using `fortify.shape()`:
dat.coast <- fortify.shape(dat.coast) # a 410951x8 dataframe
# Specify the desired domain (Pacific Ocean :
dat.coast$long[dat.coast$long >0]<-  dat.coast$long[dat.coast$long >0] - 360 

#### EEZ SHAPE file
print("Attention: **Need to download EEZ shape from http://www.marineregions.org/downloads.php  **")
#DOWNLOAD FROM HERE:
#http://www.marineregions.org/download_file.php?name=World_EEZ_v10_20180221.zip
# Reference
#Flanders Marine Institute (2018). Maritime Boundaries Geodatabase: Maritime Boundaries and Exclusive Economic Zones (200NM), version 10. Available online at http://www.marineregions.org/ https://doi.org/10.14284/312
EEZ.file<-"eez_boundaries_v10"
dat.eez <- readOGR(dsn = EEZ.path,layer = file_path_sans_ext(EEZ.file))
#summary(dat.eez)
EEZ<-fortify.shape(dat.eez)
EEZ$long[EEZ$long >0]<-  EEZ$long[EEZ$long >0] - 360 # make all latitudes negative

