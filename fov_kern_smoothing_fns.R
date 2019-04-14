###-----------------------###
#### Kernel & Focus Cone ####
#### Calculation         ####
###-----------------------###

# Size of object in focus can be calculated from FOV equations. 
# Distance is assumed to be constant, though it varies as participants 
# Area of base of a triangle also varies - though the effect isn't large enough
# for it to affect the current data (~5% of the total length of the base).

# helper f-ns for kernel construction:
# rep.row & rep.col duplicate vectors (imitating meshgrid in MATLAB)

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# to.odd rounds to the nearest odd number

to.odd<-function(x){
  y = 2 * floor(x/2)+1
  return(y)
}


# Returns Radius for images that have been scaled down (size = 432 x 243)
# change default parameters if used for another experiment
# d = participants distance from the screen

focus_radius <- function(V){
  d = 0.5                     # default distance
  cons = 3779.5275590551      # m to pixels
  V <- V*pi/180               # convert to rad
  S = d*(tan(V))              # radius of area in focus
  r = S*cons                  # radius in pixels
  r = r/4                     # scale down for small size image 
  return(r)
}

# gKern(V) creates a Gaussian kernel for a chosen vision span V, given in degrees
# more information in the papers below.

gKern<-function(V){
  r <- focus_radius(V)
  k_sigma = r/3
  gRad  = ceiling(k_sigma)
  # make a gausian kernel of this size
  vGaussianKernel = exp(-((-gRad:gRad) ^ 2) / (2 * k_sigma * k_sigma));
  vGaussianKernel = vGaussianKernel / sum(vGaussianKernel);
  
  X=rep.col(-gRad:gRad, (length(-gRad:gRad)))
  Y=rep.row(-gRad:gRad, (length(-gRad:gRad)))
  
  sigSq = k_sigma^2
  G=exp(-((X^2)+(Y^2) ) /(2* sigSq) ) /(2* pi * sigSq)
}

# angular size in degrees
# foveal<-1.5; central<-2.5; near_peripheral<-10; peripheral<-20
# https://www.jstor.org/stable/40575089?seq=1#page_scan_tab_contents
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5127899/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2851849/
# ~2.5 deg around the fovea = foveal cone
# ~2-5 deg  = central/parafoveal cone
# ~10 deg = near peripheral
# ~20 deg = peripheral (possible to read ~128 words per minute)

