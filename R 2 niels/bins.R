
#### Divide data in 10 bins with same proportion of original set of data
bins_crossvalidation <- function(data, level1, level2 ){

data$id <- 1:nrow(data)

deviated = data %>%  
  filter(Polarity==level1)%>% 
  dplyr::select(id)
deviated$id <- sample(deviated$id) 

straight = data %>% 
  filter(Polarity==level2)%>%
  dplyr::select(id)
straight$id <- sample(straight$id) 

bins  <- rep(1:10, nrow(straight) / 10)
bins_straight <- split(straight, bins)

bins  <- rep(1:10, nrow(deviated) / 10)
bins_deviated <- split(deviated, bins)


bin1 <- rbind(bins_straight$`1`, bins_deviated$`1`)
bin2 <- rbind(bins_straight$`2`, bins_deviated$`2`)
bin3 <- rbind(bins_straight$`3`, bins_deviated$`3`)
bin4 <- rbind(bins_straight$`4`, bins_deviated$`4`)
bin5 <- rbind(bins_straight$`5`, bins_deviated$`5`)
bin6 <- rbind(bins_straight$`6`, bins_deviated$`6`)
bin7 <- rbind(bins_straight$`7`, bins_deviated$`7`)
bin8 <- rbind(bins_straight$`8`, bins_deviated$`8`)
bin9 <- rbind(bins_straight$`9`, bins_deviated$`9`)
bin10 <- rbind(bins_straight$`10`, bins_deviated$`10`)

rm(bins_straight,bins_deviated)

bins <<- list(bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10)

rm(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10)

}

