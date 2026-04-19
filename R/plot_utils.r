
dim.x=710
dim.y=640

swiss.corners=c(485,75,835,295)
radar.corners <- c(255,-160,965,480)
h=75
INCA.corners=c(swiss.corners[1]-h,swiss.corners[2]-h,swiss.corners[3]+h,swiss.corners[4]+h)
x.INCA.1 = INCA.corners[1] - radar.corners[1]                                           
x.INCA.2 = x.INCA.1 + (INCA.corners[3]-INCA.corners[1])                                               
y.INCA.1 = INCA.corners[2] - radar.corners[2]                                           
y.INCA.2 = y.INCA.1 + (INCA.corners[4]-INCA.corners[2])                                               
n = (x.INCA.2-x.INCA.1+1) * (y.INCA.2-y.INCA.1+1)

INCA.corners <- c(410,0,910,370)

# Set dimensions
dim.x <- 710
dim.y <- 640
input.shape <- "/users/antoumos/SHAPE/"
input.shape.file <- "CHE_adm0"
input.dem.file <- "/users/antoumos/ccs4.png"

# Load the Switzerland shapefile
Switzerland <- readOGR(input.shape, input.shape.file, verbose = FALSE)

t <- geocors.trafo(x = Switzerland@polygons[[1]]@Polygons[[1]]@coords[,1],
                  y = Switzerland@polygons[[1]]@Polygons[[1]]@coords[,2],
                  from.type = "lonlat", to.type = "swisscors")

# Transform coordinates
t[[1]] <- t[[1]] / 1000
t[[2]] <- t[[2]] / 1000
tt <- matrix(c(t[[1]], t[[2]]), ncol = 2)
Switzerland@polygons[[1]]@Polygons[[1]]@coords <- tt
Switzerland@bbox[[1]] <-  255
Switzerland@bbox[[2]] <- -165
Switzerland@bbox[[3]] <-  255 + dim.x
Switzerland@bbox[[4]] <- -160 + dim.y
proj4string(Switzerland) <- ""

# Load DEM data
dem.m <- readPNG(input.dem.file)

# Adjust dimensions if necessary
dem.m <- t(dem.m)[,dim(dem.m)[[1]]:1]		
dem.m <- assignCoordsToPrecip(dem.m)


		
color.ramp=c( "#FFFFFF00"
				 ,"#640064","#AF00AF","#DC00DC","#3232C8","#0064FF"
				 ,"#009696","#00C832","#64FF00","#96FF00","#C8FF00"
				 ,"#FFFF00","#FFC800","#FFA000","#FF7D00","#E11900","#000000" )
				 				
color.levels=c( 0.00, 0.16, 0.25, 0.40, 0.63								
				   ,1.00, 1.60, 2.50, 4.00, 6.30								
				   ,10.0, 16.0, 25.0, 40.0, 63.0								
				   ,100., 160., 250.)	



## rotate ##
dem.fixed <- t(dem.m[nrow(dem.m):1, ])
dem.fixed <- dem.m[nrow(dem.m):1, ncol(dem.m):1]

### do this 3 times?
dem.fixed <- t(dem.fixed[nrow(dem.fixed):1, ])
dem.fixed <- dem.fixed[nrow(dem.fixed):1, ncol(dem.fixed):1]
			

# ── background-only function ───────────────────────────────────────────────
topography_background <- function() {
  plot.image(
    x = dem.m,
    time.stamp = NULL,
    color.ramp = color.ramp,
    color.levels = color.levels,
    dem.m = dem.m,
    Switzerland = Switzerland,
    scale.title = "Topography",
    mlayout = TRUE,
    plot.dem = TRUE,
    plot.color.scale = TRUE,
    border.col = "black"
  )
}