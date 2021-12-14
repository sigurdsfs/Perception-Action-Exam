# PACKAGES ####

library(reshape2)
library(devtools)
#require(rocauc)
require(pROC)
require(tidyverse)
require(plyr)
require(Hmisc)
require(ggplot2)
library(ellipse)
library(boot)
library(reshape)
library(ggpubr)
library(GGally)
require(lme4)
library(MASS) # NB: this will mask dplyr::select


## palette ####
cbPalette <- c("#000000", "#199958", "#CC1919", 
               "#00CDCD", "#7E10AF", "#AF7E10", "#1041AF", "#808080")
myPalette <- c("#000099",  "#FF3333", "#CC0000", "#990000")


### FUNCTIONS ####

# Standard error of the mean
se <- function(x, ...) {
  n <- length(x)
  return(sd(x, ...)/sqrt(n))
}

# Area under the curve
auc <- function(x,y) {
  # order x and y
  ord.x = sort(x)
  ord.y = y[rank(x)]
  
  sum(diff(ord.x)*(ord.y[-1]+ord.y[-length(ord.y)]))/2 
}

# Beta distribution parameters
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

##Plots
source("plot_measures.R")
source("plot_single_item.R")

# area under the ROC curve
auc_roc <- function(data, indices, score, label, n){
  n_data <- if(n!=FALSE) subset(data, Subject %in% sample(data$Subject, n)) else data 
  data <- n_data[indices,]
  score.te <- unlist(data[,score])  #Measureen
  label.te <- unlist(data[,label]) #Polarity
  roc.te <- roc(label.te, score.te)
  return(roc.te$auc)
}

#randomization within
random.within <- function (data, class) {
  data$Subject <- factor(data$Subject)
  vector.subject <- c('0')
  for (s in levels(data$Subject)){
    data.sub <- subset(data, Subject==s)
    tr <- unlist(data.sub[,class])
    tr.random <- as.character(sample(tr, length(tr), replace=FALSE))
    vector.subject <- c(vector.subject,tr.random)}
  vector.subject <- vector.subject[2:length(vector.subject)]
  data$Random <- vector.subject 
  return(data)
}

# power
power <-function(value){length(which(value> .5))/1000}

#Random classifier function
random_classifier <- function(training, test, iterations){
  # 1. Extract empirical frequency in training set for each class
  class1 <<- sum(training$Polarity=="straight")/length(training$Polarity)
  class2 <<- sum(training$Polarity=="deviated")/length(training$Polarity)
  
  # 2. Assign labels randomly based on probability
  labels <- c('straight', 'deviated')
  #test$random_classifier <- sample(labels, length(test$Polarity), replace=TRUE, prob=c(class1,class2))
  parameters <- estBetaParams(class1, 0.1)
  #test$random_classifier <- rbeta(length(test$Polarity), parameters$alpha, parameters$beta)
  my <<- rbeta(length(test$Polarity)*iterations, parameters$alpha, parameters$beta)
  calibrationTest <<- test
}

# Divide data in 10 bins with same proportion of original set of data
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

multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL, title = NULL, 
                      fontsize = 14, fontfamily = "Helvetica", fontface = "bold") {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (length(title)>0){
    layout <- rbind(rep(0, ncol(layout)), layout)
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                               ncol(layout), 
                                               heights = if (length(title)>0) {unit(c(0.5, rep(5,nrow(layout)-1)), "null")}
                                               else {unit(c(rep(5, nrow(layout))), "null")})))
    if(length(title) > 1){
      ncols <- 1:ncol(layout)
      for(i in seq(ncols)){
        grid.text(title[i], 
                  vp = viewport(layout.pos.row = 1, layout.pos.col = i),
                  gp = gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fontface))
      }
    } else {
      grid.text(title, 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)),
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fontface))
    }
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

