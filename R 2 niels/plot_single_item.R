#Comparing some trials to test

plot_single_subject <- function (subject, item) {
  
  p1 <- ggplot(subset(normalized_positions.plot, Subject==subject & Item.number ==item), aes(x=Time.Step, y=Acceleration_Smooth, color=Polarity, group=grp)) +
    geom_line(alpha=.5, color='green') +
    geom_point(color='green', alpha=.5)+
    ggtitle('')+
    theme_minimal() +
    expand_limits(y=c(-.1,.1)) +
    theme(legend.position = "none")
  
  p2 <- ggplot(subset(normalized_positions.plot, Subject==subject & Item.number ==item), aes(x=Time.Step, y=Acceleration, group=grp)) +
    geom_line(alpha=.5, color='blue') +
    geom_point(color='blue', alpha=.5)+
    ggtitle('')+
    theme_minimal() +
    expand_limits(y=c(-.1,.1)) +
    theme(legend.position = "none")
  
  p3 <- ggplot(subset(normalized_positions.plot, Subject==subject & Item.number ==item), aes(x=Time.Step, y=X.Position, color=Polarity, group=grp)) +
    geom_line(alpha=.5, color='red') +
    geom_point(color='red', alpha=.5)+
    ggtitle('')+
    theme_minimal() +
    expand_limits(y=c(-1.5,1.5)) +
    theme(legend.position = "none")
  
  p4 <- ggplot(subset(normalized_positions.plot, Subject==subject & Item.number ==item), aes(x=Time.Step, y=LogRatio, color=Polarity, group=grp)) +
    geom_line(alpha=.5, color='black') +
    geom_point(color='black', alpha=.5)+
    ggtitle('')+
    theme_minimal() +
    theme(legend.position = "none")
  
  return (multiplot(p1,p2,p3,p4))
}


