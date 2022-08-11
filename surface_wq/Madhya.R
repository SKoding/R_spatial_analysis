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
madhya <- st_read('shps/madhya_pradesh.shp')
madhya_wq <- read.csv('Madhya_surfacewq.csv') 
#dataframes for Sodium and potassium
madhya.sodium<- data.frame(madhya_wq$Amount.of.Sodium..Na.,madhya_wq$Surface.Water.Station.Latitude,madhya_wq$Surface.Water.Station.Longitude)
madhya.potassium<- data.frame(madhya_wq$Amount.of.Potassium..K.,madhya_wq$Surface.Water.Station.Latitude,madhya_wq$Surface.Water.Station.Longitude)
#descriptive analysis
summary(madhya.sodium)
summary(madhya.potassium)
hist(madhya.sodium$madhya_wq.Amount.of.Sodium..Na., breaks=30, col = 'green', main = 'Sodium Histogram', xlab = 'Sodium')
hist(madhya.potassium$madhya_wq.Amount.of.Potassium..K., breaks=30, col = 'blue', main = 'Potassium Histogram', xlab = 'Potassium')
##from the histogram, data for both parameters is skewed to the right, there is therefore need to normalize > Gaussian distribution
madhya.sodium['log_Na']<-log(madhya.sodium$madhya_wq.Amount.of.Sodium..Na.)
madhya.potassium['log_K']<-log(madhya.potassium$madhya_wq.Amount.of.Potassium..K.)

hist(madhya.sodium$log_Na, breaks=30, col = 'green', main = 'LogNa Histogram', xlab = 'Sodium')
hist(madhya.potassium$log_K, breaks=30, col = 'green', main = 'LogK Histogram', xlab = 'Potassium')

#renaming dataframes for easier references
names(madhya.sodium) <- c('sodium','latitude','longitude','logNa')
names(madhya.potassium) <- c('potassium','latitude','longitude','logK')
#removing NAs
madhya.sodium <- na.omit(madhya.sodium)
madhya.potassium <- na.omit(madhya.potassium)
#making the dataframe spatial
coordinates(madhya.sodium) <- c('longitude','latitude')  #sodium
coordinates(madhya.potassium) <- c('longitude','latitude') #potassium
class(madhya.potassium) #checking class type of DF

##plotting DF

madhya.sodium %>% as.data.frame %>% 
  ggplot() + 
  geom_point(aes(longitude,latitude,size=sodium,color = sodium)) +
  scale_color_gradientn(colors = c("blue","green", "yellow", "red"))+
  ggtitle("madhya Pradesh Sodium Concentrations", subtitle = 'Surface water amount of Sodium in madhya Pradesh state') + 
  theme_bw() + geom_sf(data = madhya,fill=NA) + 
  xlab('Longitude (Degrees)') + ylab('Latitude (Degrees)') + 
  labs(size = 'Sodium conc (µg/l)')

madhya.potassium %>% as.data.frame %>% 
  ggplot() + 
  geom_point(aes(longitude,latitude,size=potassium)) + 
  ggtitle("madhya Pradesh Potassium Concentrations", subtitle = 'Surface water amount of Potassium in madhya Pradesh state') + 
  theme_bw() + geom_sf(data = madhya,fill=NA) + 
  xlab('Longitude (Degrees)') + ylab('Latitude (Degrees)') + 
  labs(size = 'Potassium conc (µg/l)')

spplot(madhya.sodium, "sodium", colorkey = TRUE, main='Sodium Plot')
bubble(madhya.sodium, "sodium", key.space = "bottom")
#creating GRID for Interpolation operations
x.range <- c(74,83.0)
y.range <- c(21,27.0)

grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=0.02), y=seq(from=y.range[1], to=y.range[2], by=0.02))
coordinates(grd) <- ~ x+y #spatial
gridded(grd) <- TRUE
plot(grd, cex=1.5)
points(madhya.sodium)

####INTERPOLATION####
######IDW#
Na.idw <- idw(sodium~1,madhya.sodium,grd) #sodium IDW
plot(Na.idw)

K.idw <- idw(potassium~1,madhya.potassium,grd) #Potassium IDW
plot(K.idw)
##cropping to AOI extents
##--Sodium
idw <- raster(Na.idw)
idw_crop <- raster::mask(idw,madhya)
idw_rm <- rasterToPoints(idw_crop, spatial = TRUE)#convert back
gridded(idw_rm) <- TRUE

idw.out2 <- as.data.frame(idw_rm)
glimpse(idw.out2)
idw.out2df <- idw.out2 %>% filter(!is.na(var1.pred))
glimpse(idw.out2df)
##--Potassium
K.idw <- raster(K.idw) 
K.idw_crop <- raster::mask(K.idw,madhya)
K.idw_rm <- rasterToPoints(K.idw_crop, spatial = TRUE)
gridded(K.idw_rm) <- TRUE

K.idw.out2 <- as.data.frame(K.idw_rm)
glimpse(K.idw.out2)
K.idw.out2df <- K.idw.out2 %>% filter(!is.na(var1.pred))
glimpse(K.idw.out2df)
##plotting IDW
palette <- c('black','blue','green','grey','yellow','orange','red')
#Sodium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = idw.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Inverse Distance Weight', subtitle = 'IDW interpolation for Sodium in madhya Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#Potassium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = K.idw.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Inverse Distance Weight', subtitle = 'IDW interpolation for Potassium in madhya Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

##############Nearest Neighbor
#Sodium
Na.NN <- idw(formula = sodium ~ 1, madhya.sodium,grd,nmax = 4, nmin = 1, debug.level = -1)
plot(Na.NN)
NN.Na <- raster(Na.NN) #cropping
NN.Na_crop <- raster::mask(NN.Na,madhya)
NN.Na_rm <- rasterToPoints(NN.Na_crop, spatial = TRUE)
gridded(NN.Na_rm) <- TRUE

NN.Na.out2 <- as.data.frame(NN.Na_rm)
glimpse(NN.Na.out2)
NN.Na.out2df <- NN.Na.out2 %>% filter(!is.na(var1.pred))

#potassium
K.NN <- idw(formula = potassium ~ 1, madhya.potassium,grd,nmax = 1, nmin = 1, debug.level = -1)
plot(K.NN)
NN.K <- raster(K.NN) #cropping
NN.K_crop <- raster::mask(NN.K,madhya)
NN.K_rm <- rasterToPoints(NN.K_crop, spatial = TRUE)
gridded(NN.K_rm) <- TRUE

NN.K.out2 <- as.data.frame(NN.K_rm)
glimpse(NN.K.out2)
NN.K.out2df <- NN.K.out2 %>% filter(!is.na(var1.pred))

##plotting Nearest Neighbor (NN)
pal <- c('blue','darkgreen','green','grey','yellow','red')
#Sodium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = NN.Na.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = pal)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Nearest Neighbor (NN)', subtitle = 'NN interpolation for Sodium in madhya Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#Potassium
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = NN.K.out2df, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = pal)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Nearest Neighbor (NN)', subtitle = 'NN interpolation for Potassium in madhya Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

######Thin Plate Spline
#subsetting
madhya.sodium.tps<- data.frame(madhya_wq$Amount.of.Sodium..Na.,madhya_wq$Surface.Water.Station.Latitude,madhya_wq$Surface.Water.Station.Longitude)
madhya.potassium.tps<- data.frame(madhya_wq$Amount.of.Potassium..K.,madhya_wq$Surface.Water.Station.Latitude,madhya_wq$Surface.Water.Station.Longitude)

names(madhya.potassium.tps) <- c('potassium','latitude','longitude')
#sodium
Na_TPS <- Tps( x = as.matrix(madhya.sodium.tps[, c('madhya_wq.Surface.Water.Station.Longitude','madhya_wq.Surface.Water.Station.Latitude')]),Y = madhya.sodium.tps$madhya_wq.Amount.of.Sodium..Na.)
Nainterp_TPS <- interpolate(raster(grd), Na_TPS)
Na.TPS <- raster::mask(Nainterp_TPS, madhya)

NaTPS <- as(Na.TPS, "SpatialPixelsDataFrame")
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = as.data.frame(NaTPS), aes(x,y, fill = layer))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Thin Plate Spline (TPS)', subtitle = 'TPS interpolation for Sodium in madhya Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')
#potassium
K_TPS <- Tps( x = as.matrix(madhya.potassium.tps[, c('longitude','latitude')]),
              Y = madhya.potassium.tps$potassium)
Kinterp_TPS <- interpolate(raster(grd), K_TPS)
K.TPS <- raster::mask(Kinterp_TPS, madhya)

KTPS <- as(K.TPS, "SpatialPixelsDataFrame")
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = as.data.frame(KTPS), aes(x,y, fill = layer))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Thin Plate Spline (TPS)', subtitle = 'TPS interpolation for Potassium in madhya Pradesh', fill = 'Potassium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

###KRIGING
Na.vgm <- variogram(sodium~1, madhya.sodium,cutoff=3) #creating an Empirical variogram
K.vgm <- variogram(potassium~1, madhya.potassium, cutoff=2)
plot(K.vgm)

Na.fit <- fit.variogram(Na.vgm,vgm('Exp'), fit.method =1 ) #fit vgm
K.fit <- fit.variogram(K.vgm,vgm('Exp'), fit.method = 1)
#plot(Na.fit,cutoff=2)

plot(Na.vgm,model=Na.fit, main = 'Sodium Spherical Variogram')
plot(K.vgm,model=K.fit, main = 'Potassium Spherical Variogram')

#############Ordinary Kriging

Na.clean <- madhya.sodium[-zerodist(madhya.sodium)[,1],] #removing duplicates
Na.OK <- krige(sodium~1,Na.clean,grd,Na.fit)

K.clean <- madhya.potassium[-zerodist(madhya.potassium)[,1],]
K.OK <- krige(potassium~1,K.clean,grd,K.fit)

spplot(Na.OK)
spplot(K.OK)

ok <- brick(Na.OK)
ok <- mask(ok, madhya)
names(ok) <- c('prediction', 'variance')
plot(ok)
#cropping OK files to SHP
OK.Na <- raster(Na.OK) #sodium
OK.Na_crop <- raster::mask(OK.Na,madhya)
OK.Na_rm <- rasterToPoints(OK.Na_crop, spatial = TRUE)
gridded(OK.Na_rm) <- TRUE
OK.Na <- as.data.frame(OK.Na_rm)
glimpse(OK.Na)

OK.K <- raster(K.OK) #potassium
OK.K_crop <- raster::mask(OK.K,madhya)
OK.K_rm <- rasterToPoints(OK.K_crop, spatial = TRUE)
gridded(OK.K_rm) <- TRUE
OK.K <- as.data.frame(OK.K_rm)
glimpse(OK.K)

#Plotting OK 
ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = OK.Na, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Ordinary Kriging', subtitle = 'Ordinary Kriging map for Sodium in madhya Pradesh', fill = 'Sodium (µg/l)')+
  xlab('Longitude (Degrees)')+
  ylab('Latitude (Degrees)')

ggplot()+
  theme(panel.grid = element_line(color = 'grey70',size = 0.2), legend.position = 'right')+
  geom_tile(data = OK.K, aes(x,y, fill = var1.pred))+
  scale_fill_gradientn(colours = palette)+
  geom_sf(data = madhya,fill=NA)+
  labs(title = 'Ordinary Kriging', subtitle = 'Ordinary Kriging map for Potassium in madhya Pradesh', fill = 'Potassium (µg/l)')+
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
