# climate velocity project 6/1/25
#  this loads and preps data, then calls dVoCC (modified from VoCC::dVoCC) to find climate velocity (speed and angle)
 

library(terra)
library(climetrics)
library(data.table)


# Load the TIFF file as a raster object 
# unc-tavg data uses projection not great circle 
#class       : SpatRaster 
#dimensions  : 1476, 1003, 1  (nrow, ncol, nlyr)
#resolution  : 1000, 1000  (x, y)
#extent      : 1089000, 2092000, 4748000, 6224000  (xmin, xmax, ymin, ymax)
#coord. ref. : NZGD2000 / New Zealand Transverse Mercator 2000 (EPSG:2193)

# source paths for MDD laptop and PC
# abribitrarily chose 2 years (P=present/source year, F = future/target year) 30 years apart for testing
## MDD PC path
tiff_fileP <- rast("C:\\Users\\OEM\\Downloads\\unc-tavg-1989-02.tif") 
tiff_fileF <- rast("C:\\Users\\OEM\\Downloads\\unc-tavg-2019-02.tif") 

## MDD laptop path
## tiff_fileP <- rast("C:\\Users\\User\\Downloads\\unc-tavg-1989-02.tif") 
## tiff_fileF <- rast("C:\\Users\\User\\Downloads\\unc-tavg-2019-02.tif") 


# convert present and future tiff rasters to present and future matrices, extract cid
matP <- as.matrix(tiff_fileP, dim(tiff_fileP)[1], dim(tiff_fileP)[2])
matPcid <- matrix(1:ncell(tiff_fileP), dim(tiff_fileP)[1], dim(tiff_fileP)[2], byrow=TRUE)

matF <- as.matrix(tiff_fileF, dim(tiff_fileF)[1], dim(tiff_fileF)[2])
matFcid <- matrix(1:ncell(tiff_fileF), dim(tiff_fileF)[1], dim(tiff_fileF)[2], byrow=TRUE)
cidtiffv = 1:ncell(tiff_fileF)

## sample
## 100 x 100 sample
samplematP = matP[1200:1300, 100:200]
samplematF = matF[1200:1300, 100:200]
samplematcid = matrix(matPcid[1200:1300,100:200], 101, 101)

## 500 x 500 sample
##samplematP = matP[900:1400, 100:600]
##samplematF = matF[900:1400, 100:600]
##samplematcid = matrix(matPcid[900:1400,100:600], 501, 501)

## full raster
## samplematP = matP
## samplematF = matF
## samplematcid = matrix(matPcid)

#plot sample as heatmap
image(samplematP)

## inspect

if (1==2) {
   print(tiff_file)        # Displays basic information about the raster
   summary(tiff_file)      # Provides a summary of values and structure
   crs(tiff_file)          # Displays the Coordinate Reference System (if applicable)
   ext(tiff_file)          # Shows the spatial extent (bounding box)
   res(tiff_file)
   dim(tiff_file)
   
   nrow(tiff_file)
   ncol(tiff_file)
   nlyr(tiff_file)         # Number of layers (bands)
   
   # View detailed metadata
   tiff_info <- terra::describe(tiff_file)
   print(tiff_info)
   
   # Inspect individual values
   values <- values(tiff_file, row=1, nrows=2)  # Get values for the first 5 rows
   print(values)
}


# build climdf
# dVoCC
V1p = as.vector(samplematP)         ## present vector for dVoCC
V1f = as.vector(samplematF)         ## future vector for dVoCC
cid = 1:length(V1p)

# get cell reference data for usage and cross-checking
# NB: RC and xy are actual coords from original tiff, cidtiffvs is cid of original tiff, but cidV1p is cid restarted from 1 for V1psample

R = rowFromCell(tiff_fileF, as.vector(samplematcid)) # R of samplemat in original tiff
C = colFromCell(tiff_fileF, as.vector(samplematcid)) # C of samplemat in original tiff
X = xFromCell(tiff_fileF, as.vector(samplematcid))   # X of samplemat in original tiff
Y = yFromCell(tiff_fileF, as.vector(samplematcid))   # Y of samplemat in original tiff
dif = V1f-V1p

## dataframe passed to dVoCC
clim4VoCCdf = na.omit(data.frame(V1p = V1p, V1f = V1f, cidV1p = cid, x= X, y= Y)) #y=row, x= col      

## dataframe for cross-checking
climdf = data.frame( V1p = V1p,
                     V1f = V1f,
                     cidV1p = cid, #cidV1p is cid restarted from 1 for V1psample
                     cidOfTiff = as.vector(samplematcid), #cidtiffvs is cid of orginal tiff,
                     
                     # RC and xy are actual coords from original tiff, note R = y and C=x
                     R = R,
                     C = C,
                     x = X,
                     y = Y,
                     dif = dif
                     )

# define remaining args for dVoCC
N = 1             # one climatic variable, ie tavg
Tdiff = 30        # time difference between first/present and second/future in years
Method = "Single" # use one threshold for all observations
ClimTol = 1       # temp tolerance (threshold) to determine climate (ie temp) analogue . this is less than, so <1 is actually zero difference for this data as it is all integer numbers
                  # Just guessed this (NB median difference V1p-V1f ~ 4.3, so less than this otherwise cell matches itself in future on average)
                  # go lower to get fewer V1p matching V1f.
GeoTol = 10000   



# =1 km (in metres).  NZGD data is projected

Distfun = "Euclidean"   # projected data, so use Euclidean not great circle
#trans = NA             # not using least cost for distfun
Lonlat = FALSE          # Again, assuming projected data. 


# use (modified) dVoCC to get climate velocities

## parallel execution
## Vel = dVoCC(clim = clim4VoCCdf, n=N, tdiff=Tdiff, method=Method, climTol=ClimTol, geoTol=GeoTol, distfun= Distfun, lonlat = Lonlat )

## non-parallel execution
Vel = dVoCCNP(clim = clim4VoCCdf, n=N, tdiff=Tdiff, method=Method, climTol=ClimTol, geoTol=GeoTol, distfun= Distfun, lonlat = Lonlat )

# use these to inspect if data if putting a browse stop on a non-parallel execution
#Browse[1]> resu1 = resu[!is.na(resu$vel)]
#Browse[1]> resu3=resu[resu$focal!=resu$target & resu$vel==0]
#Browse[1]> resu4=resu[resu$focal==resu$target]
#Browse[1]> resu2=resu[resu$vel==0]
#Browse[1]> head(resu1)

