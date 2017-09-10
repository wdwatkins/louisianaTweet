#clean up
rm(list=ls())

#get cumulative precipitation by county
getPrecip <- function(states, startDate, endDate){
  
  # wg_s <- webgeom(geom = 'derivative:US_Counties', attribute = 'STATE')
  # wg_c <- webgeom(geom = 'derivative:US_Counties', attribute = 'COUNTY')
  # wg_f <- webgeom(geom = 'derivative:US_Counties', attribute = 'FIPS')
  # wi
  # county_info <- data.frame(state = query(wg_s, 'values'), county = query(wg_c, 'values'), 
  #                           fips = query(wg_f, 'values'), stringsAsFactors = FALSE) %>% 
  #   unique() 
  
  # counties_fips <- county_info %>% filter(state %in% states) %>%
  #   mutate(state_fullname = tolower(state.name[match(state, state.abb)])) %>%
  #   mutate(county_mapname = paste(state_fullname, tolower(county), sep=",")) %>%
  #   mutate(county_mapname = unlist(strsplit(county_mapname, split = " county")))
  
  counties_fips <- maps::county.fips %>% 
        mutate(statecounty=as.character(polyname)) %>% # character to split into state & county
        tidyr::separate(polyname, c('statename', 'county'), ',') %>%
        mutate(fips = sprintf('%05d', fips)) %>% # fips need 5 digits to join w/ geoknife result
        filter(statename %in% states) 
  
  stencil <- webgeom(geom = 'derivative:US_Counties',
                     attribute = 'FIPS',
                     values = counties_fips$fips)
  
  fabric <- webdata(url = 'http://cida.usgs.gov/thredds/dodsC/stageiv_combined', 
                    variables = "Total_precipitation_surface_1_Hour_Accumulation", 
                    times = c(as.POSIXct(startDate), 
                              as.POSIXct(endDate)))
  
  job <- geoknife(stencil, fabric, wait = TRUE, REQUIRE_FULL_COVERAGE=FALSE)
  check(job)
  precipData <- result(job, with.units=TRUE)
  precipData2 <- precipData %>% 
    select(-variable, -statistic, -units) %>% 
    gather(key = fips, value = precipVal, -DateTime) %>% 
    left_join(counties_fips, by="fips")
  
  return(precipData2)
  
}

#get a daily mean discharge value classified by site history percentile
getClassifiedSites <- function(states, storm.date){
  #download each state individually
  for(st in states){
    
    stDV <- renameNWISColumns(readNWISdata(service="dv",
                                           parameterCd="00060",
                                           stateCd = st,
                                           startDate = storm.date,
                                           endDate = storm.date))
    if(st != states[1]){
      storm.data <- full_join(storm.data,stDV)
      sites <- full_join(sites, attr(stDV, "siteInfo"))
    } else {
      storm.data <- stDV
      sites <- attr(stDV, "siteInfo")
    }
  }
  
  #retrieve stats data, dealing with 10 site limit to stat service requests
  reqBks <- seq(1,nrow(sites),by=10)
  statData <- data.frame()
  for(i in reqBks) {
    getSites <- sites$site_no[i:(i+9)]
    currentSites <- readNWISstat(siteNumbers = getSites,
                                 parameterCd = "00060", 
                                 statReportType="daily",
                                 statType=c("p10","p25","p50","p75","p90","mean"))
    statData <- rbind(statData,currentSites)
  }
  
  statData.storm <- statData[statData$month_nu == month(storm.date) & 
                               statData$day_nu == day(storm.date),]
  
  finalJoin <- left_join(storm.data,statData.storm)
  finalJoin <- left_join(finalJoin,sites) 
  
  #remove sites without current data 
  finalJoin <- finalJoin[!is.na(finalJoin$Flow),] 
  
  
  #classify current discharge values
  finalJoin$class <- NA
  finalJoin$class <- ifelse(is.na(finalJoin$p25), 
                            ifelse(finalJoin$Flow > finalJoin$p50_va, "greenyellow","darkorange1"),
                            ifelse(finalJoin$Flow < finalJoin$p25_va, "cyan",
                                   ifelse(finalJoin$Flow > finalJoin$p75_va, "red2","green4")))
  return(finalJoin)
}


# function to map cumulative precip data using R package maps
precipMapWithStats <- function(precipData, date, plotTitle, maxPrecip, step, 
                               classifiedSites, colSteps, fileName){
  par(mar = c(0,0,3,0))
  
  png(fileName, width = 7, height = 5, res = 150, units = 'in')
  m1 <- map('county', regions = precipData$state_fullname, col = "lightgrey")
  m2 <- map('state', regions = precipData$state_fullname, 
            add = TRUE, lwd = 1.5, col = "darkgrey")
  
  #St. Martin parish is discontinuous, so the maps package treats it as two things
  #need to duplicate the data entry
  stMartinRows <- as.data.frame(precipData[precipData$statecounty=="louisiana,saint martin parish",])
  stMartinRows <- rbind(stMartinRows, stMartinRows)
  stMartinRows$statecounty <- paste0(stMartinRows$statecounty, c(":north",":south"))
 precipData <- bind_rows(precipData, stMartinRows)
  
  # some county names are mismatched, only plot the ones that maps library 
  # knows about and then order them the same as the map
  precipData <- precipData %>%
    mutate(statecounty = gsub(x = statecounty, pattern = 'saint', replacement = 'st')) %>%  
    mutate(statecounty = gsub(x = statecounty, pattern = " parish", replacement = "")) %>%
    filter(statecounty %in% m1$names)
  
  m3 <- map('county', regions = precipData$statecounty, 
            add = TRUE, fill = TRUE, col = precipData$cols)
  par(xpd=TRUE)
  legend(x = "bottomright", inset=c(-0.3,0), fill = colSteps, cex = 0.7, bty = 'n', 
         title = "Cumulative\nPrecipitation (in)",
         legend = c(paste('<', precip_breaks[-c(1,length(precip_breaks))]), 
                    paste('>', tail(precip_breaks,2)[1]))) # greater
  # legend("topright",inset=c(-0.1,.00),
  #        legend=c("Q > P50*","Q < P50*","Q < P25","P25 < Q < P75","Q > P75"),
  #        pch=19,cex = 0.75,pt.cex = 1.2,
  #        col = c("darkorange1","greenyellow","cyan","green4","red2"),
  #        ncol = 1)
  map.scale(ratio=FALSE,cex = 0.75,
            grconvertX(.07,"npc"), 
            grconvertY(.1, "npc"))
  text("*Other percentiles not\navailable for these sites", cex=0.75,
       x=grconvertX(.97,"npc"), 
       y=grconvertY(0.68, "npc"))
  graphics::title(plotTitle,
                  line = 2, cex.main=1.2)  #title was being masked by geoknife
  mtext(side = 3, line = 0, cex = 0.9, 
        text= paste("Cumulative precipitation and daily mean discharge percentile rank\n",date))
  points(classifiedSites$dec_lon_va,
          classifiedSites$dec_lat_va, 
          col=classifiedSites$class, pch=19)
  dev.off()
}

library(dplyr)
library(tidyr)
library(geoknife) #order matters because 'query' is masked by a function in dplyr
library(RColorBrewer)
library(maps)
library(dataRetrieval)
library(lubridate)

states <- c('FL')
start <- as.Date("2017-09-07")
end <- as.Date("2017-09-10")

days <- seq.Date(start, end, by = "date")
step=.5
maxPrecip <- 3
precip_breaks <- c(seq(0, maxPrecip-step, by = step), maxPrecip)
ncolors <- length(precip_breaks) - 1
colSteps <- colorRampPalette(brewer.pal(ncolors,'Blues'))(ncolors)

totalPrecip <- NULL
for(d in as.list(days)){ #for will convert dates to numeric if not in list
  
  precipData <- getPrecip(states = states, 
                          startDate = as.character(d), 
                          endDate = as.character(d+1))
  precipData$precipVal <- precipData$precipVal/25.4 #convert mm to inches
  precipData_cols <- precipData %>% 
    group_by(statename, statecounty) %>% 
    summarize(cumprecip = sum(precipVal)) 
  if(is.null(totalPrecip)){
    totalPrecip <- precipData_cols
  }else{
    totalPrecip$cumprecip <- totalPrecip$cumprecip + precipData_cols$cumprecip
  }
  totalPrecip <- totalPrecip %>% mutate(cols = cut(cumprecip, breaks = precip_breaks, labels = colSteps, right=FALSE)) %>% 
    mutate(cols = as.character(cols))
  statsData <- getClassifiedSites(states = "LA", storm.date = d)
  statsData <- select(statsData, dec_lon_va, dec_lat_va, class)
  
  precipMapWithStats(precipData = totalPrecip, 
                     date=d,plotTitle = 'August 2016 Louisiana Storm', maxPrecip = maxPrecip, step=step, colSteps=colSteps, 
                     fileName = paste0("la_precip", d,".png"), classifiedSites = statsData)
}

#animate
system("convert -delay 80 la_precip2016* la_flood.gif")