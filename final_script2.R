library(sf)
library(tidyverse)
library(FNN)
library(lubridate)
library(ggmap)
library(AppliedPredictiveModeling)
library(caret) 
library(MASS)
library(doParallel)
library(foreach)

options(scipen = 999)
Sys.setlocale("LC_TIME", "English")

setwd("F:/SpatialAnalysisforPlanning/final")
load("final_dataset2.RData")

mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}


# read all Divvy stations
station.all <- read.csv("data/Divvy_Stations_2017_Q1Q2.csv")

# read all trips
trip.q1 <- read.csv("data/Divvy_Trips_2017_Q2.csv")

# select trips from 01/01/2017 to 01/07/2017
trip.q1$start_time <- mdy_hms(trip.q1$start_time)
train.period <- interval(mdy("05/07/2017"),mdy("05/14/2017"))

trip.in.period <- trip.q1 %>%
  filter(start_time %within% train.period)

# 150 stations closest to LOOP
# select stations where trip starts
station.start <- station.all %>%
  filter(id %in% trip.in.period$from_station_id)

# 150 nearest neighbors
station.xy <- station.start %>%
  select(longitude,latitude) %>%
  as.matrix()

loop.xy <- matrix(c(-87.625,41.8786),nrow = 1)

nn <- get.knnx(station.xy,loop.xy,k=150)

station <- station.start[nn$nn.index,]
rm(station.all,station.xy,loop.xy,nn,station.start)

# select trips from these stations
trip <- trip.in.period %>%
  filter(from_station_id %in% station$id)
rm(trip.q1,trip.in.period,train.period)

# weekday and hour of trips
trip1 <- trip %>%
  mutate(weekday=factor(wday(start_time)),
         hour=factor(hour(start_time)),
         from_station_id=factor(from_station_id)) %>%
  group_by(weekday,hour,from_station_id) %>%
  summarise(trip_count=n()) %>%
  arrange(from_station_id,weekday,hour) %>%
  as.data.frame()

# supplement all rows
for (nextStation in levels(trip1$from_station_id)){
  for (nextDay in levels(trip1$weekday)){
    for (nextHour in levels(trip1$hour)){
      if(nrow(subset(trip1,from_station_id==nextStation & weekday==nextDay & hour==nextHour))==0){
        a <- data.frame(weekday=nextDay,
                        hour=nextHour,
                        from_station_id=nextStation,
                        trip_count=0)
        trip1 <- rbind(trip1,a)
      }
    }
  }
}
trip1 <- trip1 %>%¡¡arrange(from_station_id,weekday,hour)
rm(nextDay,nextHour,nextStation,a)
ggplot()+
  geom_histogram(data = trip6,aes(trip_count),bins=30)+
  labs(x="Number of trip departures",
       y="Count",
       title="Histogram of Dependent Variable")

chicago <- get_map(location = "Chicago",
                   zoom = 13,
                   maptype = "toner-lite",
                   source = "stamen")
                   #maptype = "roadmap",
                   #source = "google")
station$id <- factor(station$id)

trip.bystation <- trip1 %>%
  group_by(from_station_id) %>%
  summarise(count_station=sum(trip_count)) %>%
  left_join(station,by=c("from_station_id"="id"))

ggmap(chicago)+
  geom_point(data = trip.bystation,aes(longitude,latitude,color=factor(ntile(count_station,5))),size=4)+
  scale_color_manual(labels=as.character(quantile(trip.bystation$count_station,
                                                  c(.1,.2,.4,.6,.8),na.rm=T)),
                     values = colorRampPalette(c("turquoise2","tomato"))(5),
                     name="Trip departures\n(Quintile breaks)")+
  labs(title="Number of Trip Departures at Each Station",
       subtitle="Since 5/7/2017 to 5/13/2017, Chicago")+
  mapTheme()

trip.byhour <- trip1 %>%
  group_by(hour) %>%
  summarise(count_hour=sum(trip_count))

ggplot(data = trip.byhour,aes(hour,count_hour,group=1))+
  geom_line(size=1.5,color="steelblue3") +
  geom_point(size=3,color="tomato")+
  labs(x="Hour in day",
       y="Trip departures",
       title="Number of Trip Departures in Each Hour")

acf(trip6$trip_count,lag.max = 24)
# trips by weekday
ggplot()+
  geom_histogram(data = trip1,aes(x=trip_count))+
  facet_wrap(~weekday)

trip.byweekday <- trip1 %>%
  group_by(weekday) %>%
  summarise(count_weekday=sum(trip_count))
trip.byweekday$weekday <- wday(seq(ymd("2017-05-07"),ymd("2017-05-13"),by="days"),label = TRUE)

ggplot(data = trip.byweekday,aes(x=weekday,y=count_weekday))+
  geom_col(width = 0.5,fill="turquoise2") +
  labs(x="Day in week",
       y="Trip departures")

# temporal & spatial autocorrelation (significant)
station.xy <- station %>%
  arrange(id) %>%
  select(longitude,latitude) %>%
  as.matrix()
nn <- get.knn(station.xy,k=10)
nn$nn.index[1:10,]
nn$nn.dist[1:10,]

trip.lag1 <- trip1 %>%
  mutate(temp.lag=0,
         spat.lag=0)
for(i in 1:150){
  for(j in 2:168){
    index <- j-168+i*168
    trip.lag1[index,]$temp.lag <- trip.lag1[index-1,]$trip_count
    
    slag.index <- nn$nn.index[i,]*168-169+j
    trip.lag1[index,]$spat.lag <- mean(trip.lag1$trip_count[slag.index])
  }
}
rm(i,j,index,slag.index)

# distance to nearest neighbors (significant)
mean.dist <- rowSums(nn$nn.dist)/10
station.neighbordist <- station %>%
  arrange(id) %>%
  mutate(near_dist=1/mean.dist) %>%
  select(id,near_dist)
rm(nn,mean.dist)
trip2 <- trip.lag1 %>%
  left_join(station.neighbordist,by=c("from_station_id"="id"))

# bus stops
bus.stops <- st_read("data/shapefile/CTA_BusStops.shp")
bus.xy <- bus.stops %>%
  as.data.frame() %>%
  select(POINT_X,POINT_Y) %>%
  as.matrix()

dist.busstop <- get.knnx(bus.xy,station.xy,k=10)
dist.busstop$nn.dist[1:10,]
mean.dist <- rowSums(dist.busstop$nn.dist)/10
near.busstop <- dist.busstop$nn.dist<mean.dist

station.busdist <- station %>%
  mutate(bus_dist=rowSums(near.busstop)) %>%
  select(id,bus_dist)
rm(bus.xy,dist.busstop,mean.dist,near.busstop)
trip4 <- trip2 %>%
  left_join(station.busdist,by=c("from_station_id"="id"))


#regression
reg4 <- glm(trip_count ~ ., family = "poisson", 
            data= trip4 %>% 
              dplyr::select(-from_station_id))
summary(reg4)

fitControl <- trainControl(method = "cv", number = 20)
set.seed(825)
lmFit4 <- train(trip_count ~ .,  
                data= trip4 %>% 
                  dplyr::select(-from_station_id),
                method = "glm", family = "poisson",
                trControl = fitControl)
lmFit4


a <- data.frame(observed=trip4$trip_count,
                predicted=reg4$fitted.values)
ggplot(data=a,aes(x=observed,y=predicted))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE)



# weather
weather.data <- data.frame(date=seq(ymd_hms("2017-05-07 0:00:00"),ymd_hms("2017-05-13 23:59:59"),by="hour"))

cl <- makeCluster(2)
registerDoParallel(cl)
results <- foreach(i=0:6,.packages="lubridate") %dopar%{
  # paste together the URL
  url.text <- paste("http://www.wunderground.com/history/airport/KORD/",
                    year(weather.data$date[i*24+1]),"/",
                    month(weather.data$date[i*24+1]),"/",
                    day(weather.data$date[i*24+1]),"/",
                    "DailyHistory.html",
                    sep="")
  # read in the webpage
  a <- scan(url.text,what="",sep="\n")
  
  # search for hourly data
  i <- grep("51 (A|P)M</td>",a)
  
  # temperature
  temp <- as.numeric(gsub("(<[^>]*>|&[^;]*;|[ F])","",a[i+2]))
  
  # precipitation
  preci <- as.numeric(gsub("(<[^>]*>|&[^;]*;|[ in]|\\t)","",a[i+25]))
  # NA to 0
  preci[is.na(preci)] <- 0
  
  return(cbind(temp,preci))
}
results <- do.call(rbind,results)
stopCluster(cl)

weather.data$temperature <- results[,1]
weather.data$precipitation <- results[,2]

trip4$temperature <- rep(weather.data$temperature,150)
trip4$precipitation <- rep(weather.data$precipitation,150)
rm(weather.data)

reg5 <- glm(trip_count ~ ., family = "poisson", 
            data= trip4 %>% 
              dplyr::select(-from_station_id))
summary(reg5)
AIC(reg4,reg5)



# remove 0am to 4am
ggplot(data = trip.byhour,aes(hour,count_hour))+
  geom_col()

trip5 <- trip4 %>%
  filter(hour %in% factor(5:23))

ggplot()+
  geom_histogram(data = trip5,aes(trip_count),binwidth = 1)

reg6 <- glm(trip_count ~ ., family = "poisson", 
            data= trip5 %>% 
              dplyr::select(-from_station_id))
summary(reg6)

reg7 <- glm.nb(trip_count ~ .,
               data= trip5 %>% 
                 dplyr::select(-from_station_id))
summary(reg7)
AIC(reg6,reg7)

# bike routes
# project the station and bike routes shapefiles using NAD_1983_StatePlane_Illinois_East_FIPS_1201_Feet projection
# use line density tool to calculate the density of bike lane and output a raster
# use extract values to points to extract cell values to station points, save as "station_bikedensity.shp"
bike.lane <- st_read("data/shapefile/station_bikedensity.shp") %>%
  as.data.frame() %>%
  dplyr::select(id,RASTERVALU) %>%
  mutate(id=factor(id))

colnames(bike.lane)[2] <- "bike_lane_density"

trip6 <- trip5 %>%
  left_join(bike.lane,by=c("from_station_id"="id"))

#reg8 <- glm.nb(trip_count ~ .,
#               data= trip6 %>% 
#                 dplyr::select(-from_station_id))
reg8 <- glm(trip_count ~ ., family = "poisson",
            data= trip6 %>% 
              dplyr::select(-from_station_id))
summary(reg8)

b <- summary(reg8)$coef[-1,1][25:31]
sx <- apply(reg8$model[4:10],2,sd)
sy <- sd(trip6$trip_count)
beta <- b * sx/sy
standardized <- data.frame(absCoef=abs(beta),
                           variable=names(beta))
ggplot(standardized, aes(x=reorder(variable,-absCoef), y=absCoef, fill=reorder(variable,-absCoef))) + 
  geom_bar(stat="identity",width = 0.6)+
  scale_fill_manual(values = heat.colors(7,1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        legend.position="none")+
  labs(x="Predictors",y="Coefficient",
       title="Absolute Standardized Coefficients for Predictors")

library(corrplot)
M <- cor(trip6 %>% dplyr::select(-trip_count,-weekday,-hour,-from_station_id))
corrplot(M, method = "number")

# RMSLE
trainset <- data.frame(observed=trip6$trip_count,
                       predicted=reg8$fitted.values) %>%
  mutate(logerror=log(predicted+1)-log(observed+1))
rmsle.train <- (mean(trainset$logerror^2))^(1/2)
testset <- data.frame(observed=trip.test$trip_count,
                      predicted=predict.glm(reg8,trip.test,type = "response"))
testset$predicted[which(testset$predicted>100)] <- 100
testset <- testset  %>%
  mutate(logerror=log(predicted+1)-log(observed+1))
rmsle.test <- (mean(testset$logerror^2))^(1/2)
data.frame(rmsle.training=rmsle.train,
           rmsle.test=rmsle.test)


# cross validation
fitControl <- trainControl(method = "cv", number = 20)
set.seed(825)
lmFit1 <- train(trip_count ~ .,  
                data= trip6 %>% 
                  dplyr::select(-from_station_id),
                method = "glm", family = "poisson",
                trControl = fitControl)
lmFit1$resample

cv <- gather(lmFit1$resample,"statistics","value",c(RMSE,MAE))
ggplot(cv, aes(value)) + 
  geom_histogram(bins = 5) +
  facet_wrap(~statistics,scales = "free_x") +
  labs(x="Value",
       y="Count",
       title="Histogram of Statistics in Cross-Validation")

# test set
test.period <- interval(mdy("05/14/2017"),mdy("05/16/2017"))
trip.q1$start_time <- mdy_hms(trip.q1$start_time)

trip.test <- trip.q1 %>%
  filter(start_time %within% test.period) %>%
  filter(from_station_id %in% station$id) %>%
  mutate(weekday=factor(wday(start_time)),
         hour=factor(hour(start_time)),
         from_station_id=factor(from_station_id)) %>%
  group_by(weekday,hour,from_station_id) %>%
  summarise(trip_count=n()) %>%
  arrange(from_station_id,weekday,hour) %>%
  as.data.frame()

for (nextStation in levels(trip.test$from_station_id)){
  for (nextDay in levels(trip.test$weekday)){
    for (nextHour in levels(trip.test$hour)){
      if(nrow(subset(trip.test,from_station_id==nextStation & weekday==nextDay & hour==nextHour))==0){
        a <- data.frame(weekday=nextDay,
                        hour=nextHour,
                        from_station_id=nextStation,
                        trip_count=0)
        trip.test <- rbind(trip.test,a)
      }
    }
  }
}
trip.test <- trip.test %>%¡¡arrange(from_station_id,weekday,hour)
rm(nextDay,nextHour,nextStation,a,trip.q1,test.period)

# time and spatial lag
trip.test <- trip.test %>%
  mutate(temp.lag=0,
         spat.lag=0)
for(i in 1:150){
  for(j in 2:48){
    index <- j-48+i*48
    trip.test[index,]$temp.lag <- trip.test[index-1,]$trip_count
    
    slag.index <- nn$nn.index[i,]*48-49+j
    trip.test[index,]$spat.lag <- mean(trip.test$trip_count[slag.index])
  }
}
rm(i,j,index,slag.index)

# distance to nearest neighbors (significant)
trip.test <- trip.test %>%
  left_join(station.neighbordist,by=c("from_station_id"="id"))

# near bus stations
trip.test <- trip.test %>%
  left_join(station.busdist,by=c("from_station_id"="id"))

# weather
weather.data <- data.frame(date=seq(ymd_hms("2017-05-14 0:00:00"),ymd_hms("2017-05-15 23:59:59"),by="hour"))

cl <- makeCluster(2)
registerDoParallel(cl)
results <- foreach(i=0:1,.packages="lubridate") %dopar%{
  # paste together the URL
  url.text <- paste("http://www.wunderground.com/history/airport/KORD/",
                    year(weather.data$date[i*24+1]),"/",
                    month(weather.data$date[i*24+1]),"/",
                    day(weather.data$date[i*24+1]),"/",
                    "DailyHistory.html",
                    sep="")
  # read in the webpage
  a <- scan(url.text,what="",sep="\n")
  
  # search for hourly data
  i <- grep("51 (A|P)M</td>",a)
  
  # temperature
  temp <- as.numeric(gsub("(<[^>]*>|&[^;]*;|[ F])","",a[i+2]))
  
  # precipitation
  preci <- as.numeric(gsub("(<[^>]*>|&[^;]*;|[ in]|\\t)","",a[i+25]))
  # NA to 0
  preci[is.na(preci)] <- 0
  
  return(cbind(temp,preci))
}
results <- do.call(rbind,results)
stopCluster(cl)

weather.data$temperature <- results[,1]
weather.data$precipitation <- results[,2]

trip.test$temperature <- rep(weather.data$temperature,150)
trip.test$precipitation <- rep(weather.data$precipitation,150)
rm(weather.data,cl)

# remove 0am to 4am
trip.test <- trip.test %>%
  filter(hour %in% factor(5:23))

# bike lane density
trip.test <- trip.test %>%
  left_join(bike.lane,by=c("from_station_id"="id"))

# predict
trip.predict <- trip.test %>%
  mutate(prediction=predict.glm(reg8,trip.test,type = "response"),
         error=prediction-trip_count)
trip.predict$prediction[which(trip.predict$prediction>100)] <- 100
acf(trip.predict$error)

# sum by hour
trip.predict.hour <- trip.predict %>%
  group_by(hour) %>%
  summarise(actualCount=sum(trip_count),
            predictedCount=sum(round(prediction))) %>%
  gather("method","value",actualCount:predictedCount)

ggplot(trip.predict.hour,aes(hour,value)) +
  geom_col(aes(fill=method),position = "dodge") +
  scale_fill_manual(values = c("turquoise2","tomato"),
                    labels = c("Observed","Predicted"))+
  labs(x="Hour in day",
       y="Sum of departures",
       title="Observed and Predicted Trip by Hour")

ggplot(trip.predict.hour,aes(hour,value,group=method)) +
  geom_line(aes(color=method),size=2) +
  scale_color_manual(values = c("turquoise2","palevioletred1"),
                     labels=c("Actual number of departures",
                              "Predicted number of departures"),
                     name="") +
  labs(x="Hour in day",
       y="Sum of departures",
       title="Observed and Predicted Departures by Hour")

# by station
trip.predict.station <- trip.predict %>%
  group_by(from_station_id) %>%
  summarise(actualCount=sum(trip_count),
            predictedCount=sum(prediction)) %>%
  mutate(abserror=abs(predictedCount-actualCount),
         percent_error=abserror/actualCount) %>%
  left_join(station,by=c("from_station_id"="id")) 

# error by station
ggmap(chicago)+
  geom_point(data = trip.predict.station,
             aes(longitude,latitude,color=factor(ntile(percent_error,5)),size=factor(ntile(actualCount,5)))) +
  scale_color_manual(labels=as.character(quantile(trip.predict.station$percent_error,
                                                  c(.1,.2,.4,.6,.8),na.rm=T)),
                     values = colorRampPalette(c("turquoise2","tomato"))(5),
                     name="Percent Error") +
  scale_size_manual(labels=as.character(quantile(trip.predict.station$actualCount,
                                                 c(.1,.2,.4,.6,.8),na.rm=T)),
                    values = c(2,3,4,5,6),
                    name="Observed Departures") +
  labs(title="Observed Departures and Percent Errors at Each Station",
       subtitle="5/14/2017 to 5/15/2017, Chicago") +
  mapTheme()

ggmap(chicago)+
  geom_point(data = trip.predict.station,
             aes(longitude,latitude,color=percent_error),size=4)+
  mapTheme()

# predicted by station
ggmap(chicago)+
  geom_point(data = trip.predict.station %>% filter(weekday=="1"),
             aes(longitude,latitude,color=factor(ntile(predictedCount,5))),size=3)+
  scale_color_manual(labels=as.character(quantile(trip.predict.station$predictedCount,
                                                  c(.1,.2,.4,.6,.8),na.rm=T)),
                     values = colorRampPalette(c("turquoise2","tomato"))(5),
                     name="Predicted Departures\n(Quintile breaks)")+
  mapTheme()

# observed by station
ggmap(chicago)+
  geom_point(data = trip.predict.station %>% filter(weekday=="1"),
             aes(longitude,latitude,color=factor(ntile(actualCount,5))),size=3)+
  scale_color_manual(labels=as.character(quantile(trip.predict.station$actualCount,
                                                  c(.1,.2,.4,.6,.8),na.rm=T)),
                     values = colorRampPalette(c("turquoise2","tomato"))(5),
                     name="Observed Departures\n(Quintile breaks)")+
  mapTheme()

# observed vs. predicted
ggplot(trip.predict %>% filter(prediction<100),aes(prediction,trip_count)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color="blue") +
  labs(x="Predicted Departures",
       y="Observed Departures",
       title="Observed Departures v.s. Predicted Values")

# error vs. predicted
ggplot(trip.predict %>% filter(prediction<100),aes(prediction,error)) +
  geom_point() +
  geom_hline(yintercept = 0,color="red") +
  labs(x="Predicted Departures",
       y="Error",
       title="Errors v.s. Predicted Values")


# compare with prediction by last week
trip.examine <- trip6 %>%
  filter(weekday %in% c("1","2"))

trip.predict <- trip.predict %>%
  mutate(predition.bylastweek=trip.examine$trip_count)

# goodness of fit
goodness.reg <- trip.predict %>%
  mutate(class=cut(percent_rank(prediction),4)) %>%
  group_by(class) %>%
  summarise(count=sum(trip_count))

goodness.lastweek <- trip.predict %>%
  mutate(class=cut(percent_rank(predition.bylastweek),4)) %>%
  group_by(class) %>%
  summarise(count=sum(trip_count))

goodness <- data.frame(group=goodness.reg$class,
                       reg_count=goodness.reg$count,
                       lastweek_count=goodness.lastweek$count) %>%
  mutate(reg_count=reg_count/sum(reg_count),
         lastweek_count=lastweek_count/sum(lastweek_count),
         group=c("0% - 25%","25% - 50%","50% - 75%","75% - 100%")) %>%
  gather(Variable, Value, reg_count:lastweek_count)

ggplot(data=goodness, aes(group,Value)) +
  geom_bar(aes(fill = Variable), position = "dodge", stat="identity") +
  scale_fill_manual(values = c("#42cbec","#f5d64e"),
                    labels=c("last week","model prediction"),
                    name="method") +
  labs(x= "Predicted Bike Demand Level",
       y="Percent of correctly predicted trip departures",
       title= "Goodness of Fit") 


