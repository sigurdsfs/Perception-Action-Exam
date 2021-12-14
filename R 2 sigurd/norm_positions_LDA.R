### Testing LDA on negation data
#Create data frame for LDA
x <- paste0('x', sprintf("%03d", c(1:101)))
y <- paste0('y', sprintf("%03d", c(1:101)))
a <- paste0('a', sprintf("%03d", c(1:101)))
v <- paste0('v', sprintf("%03d", c(1:101)))
t <- paste0('t', sprintf("%03d", c(1:101)))

# Each x and y coordenate into two columns (101 coordenates per trial) 
normalized_positions = negation_data %>%
  dplyr::select(Subject, Item.number, Polarity, Response, Normalized.positions.X,Normalized.positions.Y, Velocity, Acceleration, RawTime) %>%
  separate(Normalized.positions.Y, into= y, sep = ",") %>%
  separate(Normalized.positions.X, into= x, sep = ",") %>%
  separate(Velocity, into= v, sep = ",") %>%
  separate(RawTime, into= t, sep = ",") %>%
  separate(Acceleration, into= a, sep = ",")

normalized_positions[y] <- sapply(normalized_positions[y],as.numeric)
normalized_positions[x] <- sapply(normalized_positions[x],as.numeric)
normalized_positions[v] <- sapply(normalized_positions[v],as.numeric)
normalized_positions[a] <- sapply(normalized_positions[a],as.numeric)
normalized_positions[t] <- sapply(normalized_positions[t],as.numeric)

# Taking the negative of false items, to have everything in the same scale
normalized_positions_false = normalized_positions%>%
  filter(Response=='false')%>%
  dplyr::mutate_at(vars(starts_with('x')), funs('-'))
normalized_positions_true = filter(normalized_positions, Response=='true')
normalized_positions = bind_rows(normalized_positions_false,normalized_positions_true)
rm(normalized_positions_true, normalized_positions_false)


#More about classes
normalized_positions$Subject <- factor(normalized_positions$Subject)
normalized_positions$Polarity <- factor(normalized_positions$Polarity)
normalized_positions$Response <- factor(normalized_positions$Response)


