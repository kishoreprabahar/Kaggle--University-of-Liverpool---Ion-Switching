setwd("C:\\Users\\anush\\OneDrive\\Desktop\\fs\\dar\\finalproject")

library(data.table) 

#read the clean data from kaggel without drift and noise 

train <- fread("traink.csv")
test <- fread("testk.csv")

#removing the garbage values to save memory of pc
gc()

time <- test$time
#extending the time frame
xT1 <- seq(from = -9.9998, to = 0, by = .0001)
xT2 <- seq(from = 500.0001, to = 510.9999, by = .0001)
#dividing the time frame into two for better approach
with(train,{
  xtend1 <<- train[time >=40.0001 & time < 50,]
  xtend2 <<- train[time >=50.0001 & time < 61,] 
})
#now order the sequence
xtend1 <- xtend1[order(xtend1$time),]
xtend2 <- xtend2[order(xtend2$time),]

xtend1$time <- xT1
xtend2$time <- xT2
train <- rbind(xtend1, train, xtend2)
#sorting the test data
xT1 <- seq(from = 480.0001, to = 500.0000, by = .0001)
xT2 <- seq(from = 700.0001, to = 720.0000, by = .0001)
with(test,{
  xtend1 <<- rbind(test[time > 500 & time <= 510,],test[time > 500 & time <= 510,])
  xtend2 <<- test[time > 680 & time <= 700,]
})
xtend1 <- xtend1
xtend2 <- xtend2[order(xtend2$time),]
xtend1$time <- xT1
xtend2$time <- xT2
test <- rbind(xtend1, test, xtend2)
rm(xtend1, xtend2, xT1, xT2)
library (RcppRoll)
lags = c(11, 26, 51, 101, 1001)
CRoll <- function(DF, lags){
  
  Features = NULL
  
  for (l in lags) {
    Start = Sys.time()
    #calculating the mean
    Mn <-
      RcppRoll::roll_mean(DF$signal,
                          n = l,
                          fill = NA,
                          align = "center")
    End = Sys.time()
    print(End - Start)
  #claculating standard deviation  
    Start = Sys.time()
    Std <-
      RcppRoll::roll_sd(DF$signal,
                        n = l,
                        fill = NA,
                        align = "center")
    End = Sys.time()
    print(End - Start)
    #calculating the maximum
    Start = Sys.time()
    Mx <-
      RcppRoll::roll_max(DF$signal,
                         n = l,
                         fill = NA,
                         align = "center")
    
    End = Sys.time()
    print(End - Start)
    #calculating the minimum
    Start = Sys.time()
    Mi <-
      RcppRoll::roll_min(DF$signal,
                         n = l,
                         fill = NA,
                         align = "center")
    End = Sys.time()
    print(End - Start)
   #calculating the weighted mean 
    Start = Sys.time()
    Wm <-
      RcppRoll::roll_mean(
        DF$signal,
        n = l,
        weights = c((1:((l-1)/2)) / ((l-1)/2 * ((l-1)/2 + 1) / 1)
                    -max((1:((l-1)/2)) / ((l-1)/2 * ((l-1)/2 + 1) / 1))/((l-1)/2), 
                    max((1:((l-1)/2)) / ((l-1)/2 * ((l-1)/2 + 1) / 1))*2,
                    (((l-1)/2):1) / ((l-1)/2 * ((l-1)/2 + 1) / 1)
                    -max((1:((l-1)/2)) / ((l-1)/2 * ((l-1)/2 + 1) / 1))/((l-1)/2)),
        fill = NA,
        align = "center"
      )
    
    End = Sys.time()
    print(End - Start)
    Features <- cbind(Features, Mn, Std, Mx, Mi, Wm)
    
    
    colnames(Features)[c(
      (ncol(Features) - 4),
      (ncol(Features) - 3),
      (ncol(Features) - 2),
      (ncol(Features) - 1),
      ncol(Features))] <-
      c(
        paste("Mn", l, sep = "_"),
        paste("Std", l, sep = "_"),
        paste("Mx", l, sep = "_"),
        paste("Mi", l, sep = "_"),
        paste("Wm", l, sep = "_")
      )
    
    gc()
    print(l)
  }
  
  return(Features)
}
LRRolling <- function(DF){
  
  LR = NULL
  
  l=100
  R_W <- 
    RcppRoll::roll_mean(
      DF$signal,
      n = l,
      weights = (1:l) / (l * (l + 1) / 2),
      fill = NA,
      align = "right"
    )
  L_W <- 
    RcppRoll::roll_mean(
      DF$signal,
      n = l,
      weights = (l:1) / (l * (l + 1) / 2),
      fill = NA,
      align = "left"
    )
  l=1000
  R_Wt <- 
    RcppRoll::roll_mean(
      DF$signal,
      n = l,
      weights = (1:l) / (l * (l + 1) / 2),
      fill = NA,
      align = "right"
    )
  L_Wt <- 
    RcppRoll::roll_mean(
      DF$signal,
      n = l,
      weights = (l:1) / (l * (l + 1) / 2),
      fill = NA,
      align = "left"
    )
  
  LR <- cbind(LR, R_W, R_Wt, L_W, L_Wt)
  
  return(LR)
}
SplFeatures <- function(DF){
  
  SFeatures = NULL
  #calculating the lags
  L1 <- c(NA, DF$signal[1:(nrow(DF)-1)])
  L2 <- c(c(NA,NA), DF$signal[1:(nrow(DF)-2)])
  L3 <- c(c(NA, NA, NA), DF$signal[1:(nrow(DF)-3)])
  
  F1 <- c(DF$signal[2:nrow(DF)], NA)
  F2 <- c(DF$signal[3:nrow(DF)], c(NA, NA))
  F3 <- c(DF$signal[4:nrow(DF)], c(NA, NA, NA))
  #calculating the special signals
  dL1 <- c(NA, diff(DF$signal))
  dF1 <- c(diff(DF$signal), NA)
  sigAb <- abs(DF$signal)
  
  SigSq <- DF$signal^2
  sigSqR <- sign(DF$signal)*abs(DF$signal)^(1/2)
  SFeatures <- cbind(SFeatures, L1, L2, L3, F1, F2, F3, dL1, dF1, 
                     sigAb, 
                      SigSq, 
                     sigSqR)
  
  return(SFeatures)
}
RollF <- CRoll(train, lags)
LRF <- LRRolling(train)
spls <- SplFeatures(train)
DF_train <- cbind(train, RollF, LRF, spls)
rm(RollF, LRF, spls)
DF_train$signal_M1001 <- DF_train$signal - DF_train$Mn_1001
DF_train$signal_M101 <- DF_train$signal - DF_train$Mn_101
DF_train$batch <- DF_train$time %/%10
#checking 75% probability
b75 <- aggregate(signal~batch, data = DF_train, FUN = quantile, probs = .75)
colnames(b75)[2] <- "signal75"
DF_train <- merge(x = DF_train, y = b75, by = "batch", all.x = T)
#checking 25% of probability
b25 <- aggregate(signal~batch, data = DF_train, FUN = quantile, probs = .25)
colnames(b25)[2] <- "signal25"
DF_train <- merge(x = DF_train, y = b25, by = "batch", all.x = T)
#maxi prob
bMax <- aggregate(signal~batch, data = DF_train, FUN = max)
colnames(bMax)[2] <- "signalMax"
DF_train <- merge(x = DF_train, y = bMax, by = "batch", all.x = T)
#min prob
bMin <- aggregate(signal~batch, data = DF_train, FUN = min)
colnames(bMin)[2] <- "signalMin"
DF_train <- merge(x = DF_train, y = bMin, by = "batch", all.x = T)
#ordering data
DF_train <- DF_train[order(DF_train$time),]
#calculating upper and lower difference
DF_train$UL <- DF_train$Mx_1001 - DF_train$Mi_1001
#grouping the channels
DF_train$DD <- 0
DF_train$DD <- ifelse(DF_train$time <= 100, 1, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 100 & DF_train$time <= 150, 1, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 150 & DF_train$time <= 200, 3, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 200 & DF_train$time <= 250, 10, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 250 & DF_train$time <= 300, 5, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 300 & DF_train$time <= 350, 1, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 350 & DF_train$time <= 400, 3, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 400 & DF_train$time <= 450, 5, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 450 & DF_train$time <= 500, 10, DF_train$DD)
DF_train$DD <- ifelse(DF_train$time > 500, 1, DF_train$DD)
#comment
DF_train <- DF_train[complete.cases(DF_train),]
rm(b75, b25, bMax, bMin, train)
gc()
#formatting test data
RollF <- CRoll(test, lags)
LRF <- LRRolling(test)
specials <- SplFeatures(test)
DF_test <- cbind(test, RollF, LRF, specials)
rm(RollF, LRF, specials)
#mean calculation
DF_test$signal_M1001 <- DF_test$signal - DF_test$Mn_1001
DF_test$signal_M101 <- DF_test$signal - DF_test$Mn_101
#dividing into batches 
DF_test$batch <- DF_test$time %/%10
# the batches 75,25,min and max
b75 <- aggregate(signal~batch, data = DF_test, FUN = quantile, probs = .75)
colnames(b75)[2] <- "signal75"
DF_test <- merge(x = DF_test, y = b75, by = "batch", all.x = T)

b25 <- aggregate(signal~batch, data = DF_test, FUN = quantile, probs = .25)
colnames(b25)[2] <- "signal25"
DF_test <- merge(x = DF_test, y = b25, by = "batch", all.x = T)

bMax <- aggregate(signal~batch, data = DF_test, FUN = max)
colnames(bMax)[2] <- "signalMax"
DF_test <- merge(x = DF_test, y = bMax, by = "batch", all.x = T)

bMin <- aggregate(signal~batch, data = DF_test, FUN = min)
colnames(bMin)[2] <- "signalMin"
DF_test <- merge(x = DF_test, y = bMin, by = "batch", all.x = T)
#order the data
DF_test <- DF_test[order(DF_test$time),]
#upper and lower difference
DF_test$UL <- DF_test$Mx_1001 - DF_test$Mi_1001

#generating batch channels
DF_test$DD <- 0
DF_test$DD <- ifelse(DF_test$time <= 500, 1, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 500 & DF_test$time <= 510, 1, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 510 & DF_test$time <= 520, 3, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 520 & DF_test$time <= 530, 5, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 530 & DF_test$time <= 540, 1, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 540 & DF_test$time <= 550, 1, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 550 & DF_test$time <= 560, 10, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 560 & DF_test$time <= 570, 5, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 570 & DF_test$time <= 580, 10, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 580 & DF_test$time <= 590, 1, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 590 & DF_test$time <= 600, 3, DF_test$DD)
DF_test$DD <- ifelse(DF_test$time > 600, 1, DF_test$DD)

DF_test <- DF_test[complete.cases(DF_test),]
rm(b75, b25, bMax, bMin, test)
gc()
#converting into data frames 
DF_train <- as.data.frame(DF_train)
DF_test <- as.data.frame(DF_test)
#modifiying train and test batches 
DF_train <- DF_train[DF_train$time > 0 & DF_train$time <= 502,-which(colnames(DF_train) %in% c("time", "batch"))]
DF_test <- DF_test[DF_test$time > 500 & DF_test$time <= 700,-which(colnames(DF_test) %in% c("time", "batch"))]
gc()
#Algorithm implementation
library(keras)
install.packages("tensorflow")
library(tensorflow)
install_tensorflow()
y_train <- to_categorical(DF_train$open_channels)
#droping open channels
DF_train <-DF_train[,-c(which(colnames(DF_train) %in% c("open_channels")))]
DF_train <- as.matrix(DF_train)
DF_test <- as.matrix(DF_test)
#normalizing the data
for(c in 1:ncol(DF_train)){
  trm <- (DF_train[,c] - mean(c(DF_train[,c], DF_test[,c]))) / sd(c(DF_train[,c], DF_test[,c]))
  tem <- (DF_test[,c] - mean(c(DF_train[,c], DF_test[,c]))) / sd(c(DF_train[,c], DF_test[,c]))
  
  DF_train[,c] <- trm
  DF_test[,c] <- tem
}
rm(tem, trm)
gc()
DF_train <- array_reshape(DF_train, c(dim(DF_train), 1))
gc()
DF_test <- array_reshape(DF_test, c(dim(DF_test), 1))
gc()
#learning the epochs
learning_scheduler = function(epoch, lr) {
  if (epoch < 25) {
    return(.001-.00002*epoch)
  } else if(epoch >= 25 & epoch < 35){
    return(.0012-.00002*epoch)
  } else if(epoch >= 35 & epoch < 40){
    return(.0013-.00002*epoch)
  } else if(epoch >= 40 & epoch < 50){
    return(.0014-.00002*epoch)
  } else {
    return(.0015-.00002*epoch)
  }
}


model <- keras_model_sequential() 
model %>%layer_conv_1d(filters = 10,kernel_size = 8,strides = 1,activation='relu',padding = "same", input_shape = c(dim(DF_train)[2],1)) %>%
  layer_conv_1d(filters = 15,kernel_size = 6,dilation_rate = 8,activation='relu',padding = "same") %>%
  layer_conv_1d(filters = 20,kernel_size = 3,dilation_rate = 6,activation='relu',padding = "same") %>%
  layer_flatten() %>% 
  
  layer_dense(units = 96,activation = 'relu'
              #,regularizer_l1_l2(l1 = .01, l2 = .0001)
  ) %>%
  layer_dropout(rate = 0.05) %>%
  layer_dense(units = 32,
              activation = 'relu'
  ) %>%
  layer_dropout(rate = 0.025) %>%
  layer_dense(units = 11, activation = 'softmax')
summary(model)
model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_adam(beta_1 = .9, beta_2 = .99),
  metrics = c('accuracy')
)
gc()
history <- model %>% fit(
  DF_train, y_train, 
  epochs = 60, 
  batch_size = 50000,
  callbacks = callback_learning_rate_scheduler(learning_scheduler)
)
gc()
plot(history)
predictKeras <- model %>% predict_classes(DF_test)
#changing to data frame
predictKeras <- as.data.frame(predictKeras)
predictKeras <- predictKeras
sample_submission <- fread("sample_submission.csv", colClasses = "character")
sample_submission$open_channels <- predictKeras
fwrite(sample_submission, "submission.csv")
gc()
#plotting a graph using test data

test <- fread("testk.csv")
library(ggplot2)
sample_submission$open_channels <- as.factor(sample_submission$open_channels)
sample_submission$time <- test$time[test$time > 500 & test$time <= 700]
sample_submission$signal <- test$signal[test$time > 500 & test$time <= 700]
options(repr.plot.width = 15, repr.plot.height = 10)
ggplot(sample_submission, aes(x=time, y=signal)) +
  geom_point(shape=16, aes(color = open_channels), alpha = 0.4) +
  theme_grey(base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=10)))
gc()
