
## BASELINE LONG VS. SHORT RESPONSES 
# Just for affirmative true cases
negation_data.true_aff <- subset(negation_data, Polarity=='P' & Expected_response=='true')
negation_data.true_aff = negation_data.true_aff %>%
  dplyr::select(Subject, Normalized.positions.X, Normalized.positions.Y, Item.number, Sentence, RT, Response, Velocity, Acceleration, RawTime) 

### Division between uncertain and straight trials based on response times ####
# 1. Divide between long and short affirmative true cases. Arguably, long affirmative true cases are uncertain 
## Taking mean log response time per subject (for only those trials)
meanRT_subject <- ddply(negation_data.true_aff, c("Subject"),
                        function(negation_data.true_aff)c(RTmean=mean(log(negation_data.true_aff$RT))))

negation_data.true_aff <- merge(meanRT_subject , negation_data.true_aff, by='Subject')

## Categorizing the trials into 'straight' or 'uncertain' based on the subject mean RT
negation_data.true_aff$Class <- if_else(log(negation_data.true_aff$RT)< negation_data.true_aff$RTmean , 'straight', 'uncertain')


# 2. See how does our LDA does dividing between straight and uncertain
# LDA FULL (coordinates, velocity and euclidean-based acceleration) ####
load('LDA-Full.RData')
#Create data frame for LDA
x <- paste0('x', sprintf("%03d", c(1:101)))
y <- paste0('y', sprintf("%03d", c(1:101)))
a <- paste0('a', sprintf("%03d", c(1:101)))
v <- paste0('v', sprintf("%03d", c(1:101)))
t <- paste0('t', sprintf("%03d", c(1:101)))

# Each x and y coordenate into two columns (101 coordenates per trial) 
normalized_positions = negation_data.true_aff %>%
  dplyr::select(Subject, Item.number, Class, Response, Normalized.positions.X,Normalized.positions.Y, Velocity, Acceleration, RawTime) %>%
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

#More about classes
normalized_positions$Subject <- factor(normalized_positions$Subject)
normalized_positions$Class <- factor(normalized_positions$Class)
normalized_positions$Response <- factor(normalized_positions$Response)


normalized_positions.new <- normalized_positions %>%
  dplyr::select(Subject, Item.number, Class, Response, one_of(all_data_columns))

normalized_positions.new_pca <- bind_cols(normalized_positions.new,
                                          as.data.frame(predict(m_pca, normalized_positions.new)[,1:n_pca]))

lda_measure.new.df <- data_frame(
  lda_measure=c(as.matrix(dplyr::select(normalized_positions.new_pca, starts_with("PC"))) %*% v_lda- b_lda),
  Subject = normalized_positions.new_pca$Subject, 
  Item.number = normalized_positions.new_pca$Item.number, 
  Class = normalized_positions.new_pca$Class, 
  Response = normalized_positions.new_pca$Response)

##Including the relevant lda_measure in the data
negation_data.true_aff$Subject <- factor(negation_data.true_aff$Subject)
negation_data.true_aff$Response <- factor(negation_data.true_aff$Response)
negation_data.true_aff$Class <- factor(negation_data.true_aff$Class)
negation_data.true_aff <- dplyr::full_join(lda_measure.new.df, negation_data.true_aff, by=c("Subject", "Item.number", "Class", "Response"))

## See the trajectories
X.plot = negation_data.true_aff %>%
  dplyr::select(Subject,Class, Response, Normalized.positions.X, Item.number, Sentence) %>%
  separate(Normalized.positions.X, into= as.character(c(1:101)), sep = ",") %>%
  gather(Time.Step, X.Position, 4:104) 

Y.plot = negation_data.true_aff %>%
  dplyr::select(Subject,Class, Response, Normalized.positions.Y, Item.number, Sentence)%>%
  separate(Normalized.positions.Y, into= as.character(c(1:101)), sep = ",") %>%
  gather(Time.Step, Y.Position, 4:104) 

normalized_positions.plot_true_af <- merge(X.plot,Y.plot)
normalized_positions.plot_true_af$X.Position<- as.numeric(normalized_positions.plot_true_af$X.Position)
normalized_positions.plot_true_af$Y.Position <- as.numeric(normalized_positions.plot_true_af$Y.Position)



