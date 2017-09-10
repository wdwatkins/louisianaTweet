getPrecip <- function(states, startDate, endDate){
  
  wg_s <- webgeom(geom = 'derivative:US_Counties', attribute = 'STATE')
  wg_c <- webgeom(geom = 'derivative:US_Counties', attribute = 'COUNTY')
  wg_f <- webgeom(geom = 'derivative:US_Counties', attribute = 'FIPS')
  county_info <- data.frame(state = query(wg_s, 'values'), county = query(wg_c, 'values'), 
                            fips = query(wg_f, 'values'), stringsAsFactors = FALSE) %>% 
    unique() 
  
  counties_fips <- county_info %>% filter(state %in% states) %>%
    mutate(state_fullname = tolower(state.name[match(state, state.abb)])) %>%
    mutate(county_mapname = paste(state_fullname, tolower(county), sep=",")) %>%
    mutate(county_mapname = unlist(strsplit(county_mapname, split = " county")))
  
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

# function to map cumulative precip data using R package maps
precipMap <- function(precipData, startDate, endDate, plotTitle, maxPrecip, step, fileName){
  precip_breaks <- c(seq(0, maxPrecip-step, by = step), maxPrecip)
  ncolors <- length(precip_breaks) - 1
  cols <- colorRampPalette(brewer.pal(ncolors,'Blues'))(ncolors)
  
  precipData_cols <- precipData %>% 
    group_by(state_fullname, county_mapname) %>% 
    summarize(cumprecip = sum(precipVal)) %>% 
    mutate(cols = cut(cumprecip, breaks = precip_breaks, labels = cols, right=FALSE)) %>% 
    mutate(cols = as.character(cols))
  
  par(mar = c(0,0,3,0))
  
  png(fileName, width = 7, height = 5, res = 150, units = 'in')
  m1 <- map('county', regions = precipData_cols$state_fullname, col = "lightgrey")
  m2 <- map('state', regions = precipData_cols$state_fullname, 
            add = TRUE, lwd = 1.5, col = "darkgrey")
  
  #St. Martin parish is discontinuous, so the maps package treats it as two things
  #need to duplicate the data entry
  stMartinRows <- as.data.frame(precipData_cols[precipData_cols$county_mapname=="louisiana,saint martin parish",])
  stMartinRows <- rbind(stMartinRows, stMartinRows)
  stMartinRows$county_mapname <- paste0(stMartinRows$county_mapname, c(":north",":south"))
  precipData_cols <- bind_rows(precipData_cols, stMartinRows)
  
  # some county names are mismatched, only plot the ones that maps library 
  # knows about and then order them the same as the map
  precipData_cols <- precipData_cols %>%
    mutate(county_mapname = gsub(x = county_mapname, pattern = 'saint', replacement = 'st')) %>%  
    mutate(county_mapname = gsub(x = county_mapname, pattern = " parish", replacement = "")) %>%
    filter(county_mapname %in% m1$names)
  precipData_cols <- precipData_cols[na.omit(match(m1$names, precipData_cols$county_mapname)),]
  
  m3 <- map('county', regions = precipData_cols$county_mapname, 
            add = TRUE, fill = TRUE, col = precipData_cols$cols)
  par(xpd=TRUE)
  legend(x = "bottomright", inset=c(-0.2,0), fill = cols, cex = 0.7, bty = 'n', 
         title = "Cumulative\nPrecip (in)",
         legend = c(paste('<', precip_breaks[-c(1,length(precip_breaks))]), 
                    paste('>', tail(precip_breaks,2)[1]))) # greater
  graphics::title(plotTitle,
                  line = 2, cex.main=1.2)  #title was being masked by geoknife
  mtext(side = 3, line = 1, cex = 0.9, 
        text= paste("By county from", startDate, "to", endDate))
  dev.off()
}

library(dplyr)
library(tidyr)
library(geoknife) #order matters because 'query' is masked by a function in dplyr
library(RColorBrewer)
library(maps)

states <- c('LA')
start <- "2016-08-11 00:00:00"
end <- "2016-08-14 12:00:00"

precipData <- getPrecip(states = states, 
                        startDate = start, 
                        endDate = end)
precipData$precipVal <- precipData$precipVal/25.4 #convert mm to inches
precipMap(precipData, 
          startDate = start, 
          endDate = end,plotTitle = 'August 2016 Louisiana Storm', maxPrecip = 20, step=2, fileName = "la_precip.png")