library(purrr)
library(rgdal)
library(raster)
library(SpatialPack) # compute crh stats
library(yaml)
library(zonator)

source("src/00_lib/utils.R")

# Setup -------------------------------------------------------------------

# Provide data
provide_file_paths <- list.files("data/processed/features/provide", pattern = "\\.tif$",
                                 full.names = TRUE, recursive = TRUE)
# Datadryad data
datadryad_file_paths <- list.files("data/processed/features/datadryad",
                                   pattern = "\\.tif$", full.names = TRUE,
                                   recursive = TRUE)

# UDR
udr_file_paths <- list.files("data/processed/features/udr", pattern = "\\.tif$",
                             full.names = TRUE, recursive = TRUE)

# ES and BD files
es_file_paths <- c(provide_file_paths, datadryad_file_paths)
bd_file_paths <- udr_file_paths
all_file_paths <- c(es_file_paths)

# Read into raster --------------------------------------------------------

morans <- {}
n_files <- length(all_file_paths)
for (i in 1:n_files) {
  infile <- all_file_paths[i]
  message(" [", i, "/", n_files, "] processing ", infile)
  morans[[basename(infile)]] <- raster::Moran(raster::raster(infile))
}

## FIGURE 1 - Plot environmental services and zones of Cornwall

import12=vector("list",12)
for ( x in c(1:12) ){
  import12[[x]]=raster(readGDAL(paste("~/harmonized/",map_names[x],sep="")))
}
MAPS=brick(import12[[1]],import12[[2]],import12[[3]],import12[[4]],import12[[5]],import12[[6]]
           ,import12[[7]],import12[[8]],import12[[9]],import12[[10]],import12[[11]],import12[[12]])
import12[[5]]@data@values[which(complete.cases(import12[[1]]@data@values)==FALSE)]=NA
mycolor=topo.colors(100)
mycolor[1]="#808080" # grey
TITLE=c("Carbon in soil: 0.72","Agriculture: 0.58","Urban development: 0.47","Carbon above
ground: 0.47","Plant production: 0.41","Aesthetic: 0.38","Recreation: 0.36","Flood mitigation:
0.26","Tourism: 0.18","Wind energy: 0.16","Solar energy: 0.11","Zones of Cornwall")
nn=0
png(paste("~/figures/Figure1_maps.png",sep=""), width = 1600, height = 1200)
nf <- layout(matrix(c(1:12),3,4,byrow=TRUE), c(1,1), c(1,1), TRUE) ; layout.show(nf) ;
par(mar=c(5,5,3,5))
for ( x in c(3,1,10,2,5,7,8,4,9,11,6) ){
  nn=1+nn
  r=import12[[x]]
  plot(r, col=mycolor, legend=F, axes=F, main=TITLE[nn],cex.main=3, cex=2, cex.legend=3, box=F)
  if (x==6){plot(r, legend.only=TRUE, col=mycolor, legend.width=6, legend.shrink=0.75,
                 axis.args=list(at=c(0,50,89), labels=c(0,50,100), cex.axis=2.2))}
}
r=import12[[12]]
plot(r, col=c("blue","red","green3","orange"), legend=F, axes=F, main="Zones of
Cornwall",cex.main=3, cex=1.8, box=F)
legend(131000,110000,c("Coastal","West","Centre","East"), bty = "n", pch=c(15,15,15,15),
       pt.cex=c(2.5,2.5,2.5), cex=3, col=c("blue","red","green3","orange"))
dev.off()
################################################################################
##########
##### Compute Moran I index
MI=matrix(ncol=1,nrow=11)
rownames(MI)=c("Agriculture value","Carbon above ground","Carbon in soil","Flood
mitigation","Plant production","Solar energy","Aesthetic","Recreation","Tourism","Urban
development","Wind energy")
for ( x in c(1:12) ){
  MI[x,1]=Moran(import12[[x]])
}
MI[order(MI[,1]),]
