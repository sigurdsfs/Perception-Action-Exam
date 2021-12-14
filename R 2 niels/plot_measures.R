
##PLOTTING different measures
#Input: data= data, measure = 'measure', division = 'division'
#Output: multiplot with plots for histogram, density and bar graph taken from means + print of means + SE
require(ggplot2)

plot_measure <- function(data, measure, division){
  histogram <- ggplot(data, aes_string(x=measure, fill=division), environment = environment()) +
    geom_histogram(bins=12,  position="dodge")+
    scale_fill_brewer(palette="Set1")+
    theme_minimal() + 
    theme(legend.position = "none")
  
  density <-ggplot(data, aes_string(x=measure, fill=division), environment = environment()) +
    geom_density(alpha=.5)+
    scale_fill_brewer(palette="Set1")+
    theme_minimal() +theme(legend.position = "none")
  
  
  
  mydata.agreggated <- ddply(data, c(division, "Subject"),
                             function(x,ind){mean(x[,ind])},measure)
  
  mydata.agreggated.overall <- ddply(mydata.agreggated, c(division),
                                     function(mydata.agreggated)c(mean=mean(mydata.agreggated$V1, na.rm=T), se=se(mydata.agreggated$V1, na.rm=T) ))
  
 

  mean.plot <- ggplot(mydata.agreggated.overall, aes(x=mydata.agreggated.overall[,1], y=mean, fill=mydata.agreggated.overall[,1])) +
    geom_bar(position=position_dodge(), stat="identity") +
    scale_fill_brewer(palette="Set1")+
    theme_minimal()+
    xlab(' ')  +
    theme(legend.position = "none") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))
 
   print(mydata.agreggated.overall) 
  
   return(multiplot(density, mean.plot,
                   cols=2))  }


