library('gstat') #kriging library
library('ggplot2') #plotting and visualization library
library('sp') #
library('sf')
library('rgdal')
library('dplyr')
library('automap')
library('raster')
library('fields')
library('stars')
library('dismo')


#set working directory
setwd("~/Documents/Upwork/Spatial_analysis/surface_wq")
#data reading
andra <- st_read('shps/Andra_pradesh.shp')
andra_wq <- read.csv('andhra_surfacewq.csv') 
#dataframes for Sodium and potassium
andra.sodium<- data.frame(andra_wq$Amount.of.Sodium..Na.,andra_wq$Surface.Water.Station.Latitude,andra_wq$Surface.Water.Station.Longitude)
andra.potassium<- data.frame(andra_wq$Amount.of.Potassium..K.,andra_wq$Surface.Water.Station.Latitude,andra_wq$Surface.Water.Station.Longitude)
#descriptive analysis
summary(andra.sodium)
summary(andra.potassium)
hist(andra.sodium$andra_wq.Amount.of.Sodium..Na., breaks=30, col = 'green', main = 'Sodium Histogram', xlab = 'Sodium')
hist(andra.potassium$andra_wq.Amount.of.Potassium..K., breaks=30, col = 'blue', main = 'Potassium Histogram', xlab = 'Potassium')
##from the histogram, data for both parameters is skewed to the right, there is therefore need to normalize > Gaussian distribution
andra.sodium['log_Na']<-log(andra.sodium$andra_wq.Amount.of.Sodium..Na.)
andra.potassium['log_K']<-log(andra.potassium$andra_wq.Amount.of.Potassium..K.)

hist(andra.sodium$log_Na, breaks=30, col = 'green', main = 'LogNa Histogram', xlab = 'Sodium')
hist(andra.potassium$log_K, breaks=30, col = 'green', main = 'LogK Histogram', xlab = 'Potassium')

#renaming dataframes for easier references
names(andra.sodium) <- c('sodium','latitude','longitude','logNa')
names(andra.potassium) <- c('potassium','latitude','longitude','logK')
#removing NAs
andra.sodium <- na.omit(andra.sodium)
andra.potassium <- na.omit(andra.potassium)
#making the dataframe spatial
coordinates(andra.sodium) <- c('longitude','latitude')  #sodium
coordinates(andra.potassium) <- c('longitude','latitude') #potassium
class(andra.potassium) #checking class type of DF

##plotting DF

andra.sodium %>% as.data.frame %>% 
  ggplot() + 
  geom_point(aes(longitude,latitude,size=sodium,color = sodium)) +
  scale_color_gradientn(colors = c("blue","green", "yellow", "red"))+
  ggtitle("Andra Pradesh Sodium Concentrations", subtitle = 'Surface water amount of Sodium in andra Pradesh state') + 
  theme_bw() + geom_sf(data = andra,fill=NA) + 
  xlab('Longitude (Degrees)') + ylab('Latitude (Degrees)') + 
  labs(size = 'Sodium conc (µg/l)')

andra.potassium %>% as.data.frame %>% 
  ggplot() + 
  geom_point(aes(longitude,latitude,size=potassium),color = rainbow(316), alpha=3/4) + 
  ggtitle("Andra Pradesh Potassium Concentrations", subtitle = 'Surface water amount of Potassium in andra Pradesh state') + 
  theme_bw() + geom_sf(data = andra,fill=NA) + 
  xlab('Longitude (Degrees)') + ylab('Latitude (Degrees)') + 
  labs(size = 'Potassium conc (µg/l)')

spplot(andra.sodium, "logNa", colorkey = TRUE, main='Sodium Plot')
bubble(andra.sodium, "logNa", key.space = "bottom")
#creating GRID for Interpolation operations
x.range <- c(76.0,84.9)
y.range <- c(12.5,19.2)

grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=0.02), y=seq(from=y.range[1], to=y.range[2], by=0.02))
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE
plot(grd, cex=1.5)
points(andra.sodium)

####INTERPOLATION####
######IDW#
Na.idw <- idw(sodium~1,andra.sodium,grd) #sodium IDW
plot(Na.idw)

K.idw <- idw(potassium~1,andra.potassium,grd) #Potassium IDW
plot(K.idw)
##cropping to AOI extents
##--Sodium
idw <- raster(Na.idw)
idw_crop <- raster::mask(idw,andra)
idw_rm <- rasterToPoints(idw_crop, spatial = TRUE)
gridded(idw_rm) <- TRUE

idw.out2 <- as.data.frame(idw_rm)
glimpse(idw.out2)
idw.out2df <- idw.out2 %>% filter(!is.na(var1.pred))
glimpse(idw.out2df)
##--Potassium
K.idw <- raster(K.idw) 
K.idw_crop <- raster::mask(K.idw,andra)
K.idw_rm <- rasterToPoints(K.idw_crop, spatial = TRUE)
gridded(K.idw_rm) <- TRUE

K.idw.out2 <- as.data.frame(K.idw_rm)
glimpse(K.idw.out2)
K.idw.out2df <- K.idw.out2 %>% filter(!is.na(var1.pred))
glimpse(K.idw.out2df)
##plotting IDW
palette <- c('darkblue','blue','green','yellow','orange','red')
#Sodium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = idw.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Inverse Distance Weight', subtitle = 'IDW interpolation for Sodium in Andra Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#Potassium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = K.idw.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Inverse Distance Weight', subtitle = 'IDW interpolation for Potassium in Andra Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

##############Nearest Neighbor
#Sodium
Na.NN <- idw(formula = sodium ~ 1, andra.sodium,grd,nmax = 1, nmin = 1, debug.level = -1)
plot(Na.NN)
NN.Na <- raster(Na.NN) #cropping
NN.Na_crop <- raster::mask(NN.Na,andra)
NN.Na_rm <- rasterToPoints(NN.Na_crop, spatial = TRUE)
gridded(NN.Na_rm) <- TRUE

NN.Na.out2 <- as.data.frame(NN.Na_rm)
glimpse(NN.Na.out2)
NN.Na.out2df <- NN.Na.out2 %>% filter(!is.na(var1.pred))

#potassium
K.NN <- idw(formula = potassium ~ 1, andra.potassium,grd,nmax = 1, nmin = 1, debug.level = -1)
plot(K.NN)
NN.K <- raster(K.NN) #cropping
NN.K_crop <- raster::mask(NN.K,andra)
NN.K_rm <- rasterToPoints(NN.K_crop, spatial = TRUE)
gridded(NN.K_rm) <- TRUE

NN.K.out2 <- as.data.frame(NN.K_rm)
glimpse(NN.K.out2)
NN.K.out2df <- NN.K.out2 %>% filter(!is.na(var1.pred))

##plotting Nearest Neighbor (NN)
#Sodium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = NN.Na.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Nearest Neighbor (NN)', subtitle = 'NN interpolation for Sodium in Andra Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#Potassium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = NN.K.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Nearest Neighbor (NN)', subtitle = 'NN interpolation for Potassium in Andra Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

######Thin Plate Spline
andra.sodium.tps<- data.frame(andra_wq$Amount.of.Sodium..Na.,andra_wq$Surface.Water.Station.Latitude,andra_wq$Surface.Water.Station.Longitude)
andra.potassium.tps<- data.frame(andra_wq$Amount.of.Potassium..K.,andra_wq$Surface.Water.Station.Latitude,andra_wq$Surface.Water.Station.Longitude)

names(andra.potassium.tps) <- c('potassium','latitude','longitude')
#sodium
Na_TPS <- Tps( x = as.matrix(andra.sodium.tps[, c('andra_wq.Surface.Water.Station.Longitude','andra_wq.Surface.Water.Station.Latitude')]),Y = andra.sodium.tps$andra_wq.Amount.of.Sodium..Na.)
Nainterp_TPS <- interpolate(raster(grd), Na_TPS)
Na.TPS <- raster::mask(Nainterp_TPS, andra)

NaTPS <- as(Na.TPS, "SpatialPixelsDataFrame")
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = as.data.frame(NaTPS), aes(x,y, fill = layer))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Thin Plate Spline (TPS)', subtitle = 'TPS interpolation for Sodium in Andra Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#potassium
K_TPS <- Tps( x = as.matrix(andra.potassium.tps[, c('longitude','latitude')]),
              Y = andra.potassium.tps$potassium)
Kinterp_TPS <- interpolate(raster(grd), K_TPS)
K.TPS <- raster::mask(Kinterp_TPS, andra)

KTPS <- as(K.TPS, "SpatialPixelsDataFrame")
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = as.data.frame(KTPS), aes(x,y, fill = layer))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Thin Plate Spline (TPS)', subtitle = 'TPS interpolation for Potassium in Andra Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

###KRIGING
Na.vgm <- variogram(log(sodium)~1, andra.sodium, cutoff=2) #creating an Empirical variogram
K.vgm <- variogram(potassium~1, andra.potassium, cutoff=2)
plot(Na.vgm)

Na.fit <- fit.variogram(Na.vgm,vgm(.1,'Sph',2), fit.method = 1) #fit vgm
K.fit <- fit.variogram(K.vgm,vgm(1,'Sph',2), fit.method = 6)

plot(Na.vgm,model=Na.fit, main = 'Sodium Spherical Variogram')
plot(K.vgm,model=K.fit, main = 'Potassium Spherical Variogram')

#############Ordinary Kriging

Na.clean <- andra.sodium[-zerodist(andra.sodium)[,1],] #removing duplicates
Na.OK <- krige(sodium~1,Na.clean,grd,Na.fit)

K.clean <- andra.potassium[-zerodist(andra.potassium)[,1],]
K.OK <- krige(potassium~1,K.clean,grd,K.fit)

spplot(Na.OK)
spplot(K.OK)

ok <- brick(Na.OK)
ok <- mask(ok, andra)
names(ok) <- c('prediction', 'variance')
plot(ok)
#cropping OK files to SHP
OK.Na <- raster(Na.OK) #sodium
OK.Na_crop <- raster::mask(OK.Na,andra)
OK.Na_rm <- rasterToPoints(OK.Na_crop, spatial = TRUE)
gridded(OK.Na_rm) <- TRUE
OK.Na <- as.data.frame(OK.Na_rm)
glimpse(OK.Na)

OK.K <- raster(K.OK) #potassium
OK.K_crop <- raster::mask(OK.K,andra)
OK.K_rm <- rasterToPoints(OK.K_crop, spatial = TRUE)
gridded(OK.K_rm) <- TRUE
OK.K <- as.data.frame(OK.K_rm)
glimpse(OK.K)

#Plotting OK 
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = OK.Na, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Ordinary Kriging', subtitle = 'Ordinary Kriging map for Sodium in Andra Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = OK.K, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = andra,fill=NA)+
  labs(title = 'Ordinary Kriging', subtitle = 'Ordinary Kriging map for Potassium in Andra Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

##CROSS VALIDATION
Na.gst <- gstat(formula = log(sodium)~1, data = Na.clean,model =  Na.fit)
cv <- gstat.cv(Na.gst)
spplot(cv)

K.gst <- gstat(formula = potassium~1, data = K.clean,model =  K.fit)
K.cv <- gstat.cv(K.gst)
spplot(K.cv)

##RMSE
Na.rmse <- sqrt(sum((cv$var1.pred - cv$observed)^2) / nrow(cv))
Na.rmse
K.rmse <- sqrt(sum((K.cv$var1.pred - K.cv$observed)^2) / nrow(K.cv))
K.rmse
cat('The RMSE for Sodium is: ',Na.rmse,'and for Potassium is: ',K.rmse)

##ENSEMBLE
#~CANT RELATE TO VARIABLE BEING USED IN THE REFERENCE
