sbbg_precip <- data.frame(readxl::read_xls('./data_files/Stachys/Setup/321dailys.xls'))

names(sbbg_precip)

annual_precip <- aggregate(daily_rain ~ water_year, sbbg_precip, FUN = sum)

mean(annual_precip$daily_rain)

annual_precip[annual_precip$water_year==2017,]