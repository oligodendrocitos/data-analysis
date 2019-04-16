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


###---------------------------------###
#### Pool fixations into matrices: ####
###---------------------------------###

# increment function
inc <- function(x)
{
  eval.parent(substitute(x <- x + 1))
}

# df_subset - subset of the main df to pool fixations for.
# pools fixations from a trial into a single matrix.
# assumed to contain columns of x & y coorinates, 
# boolean coordinates indicating whether the sample is reliable (trackloss) or a saccade 
# rounds coordinates to a pixel, picks reliable samples & fixations 
# variables used: img_width = 1728 (ncol), img_height = 972 (nrow)

poolfixations <- function(df_subset, nrows, ncols, scale_coords=1){
  # initialise empty matrix:
  # variables used: img_width = 1728 (ncol), img_height = 972 (nrow)
  mat <- matrix(0, nrow = nrows, ncol = ncols)
  
  # round x and y coordinates:
  df_subset$x <- round(df_subset$x/scale_coords)
  df_subset$y <- round(df_subset$y/scale_coords)
  
  # choose aoi
  df_subset <- df_subset[df_subset$image_aoi==TRUE, ]
  #df_subset <- df_subset[df_subset$stack==TRUE, ]
  
  # add all coordinates to the empty matrix
  for (row in 1:(length(df_subset$x)) ) {
    temp <- df_subset[row, ]
    inc(mat[temp$y,temp$x])
  }
  return(mat)
}


# Crops a matrix dealing with rows/cols the same way 
# centercrop does in pytorch.
# th, tw - target height and width
# mat - matrix to crop
crop_torch <- function(mat, th, tw){
  # lower range
  i = to.odd((nrow(mat)-th)/2)
  j = to.odd((ncol(mat)-tw)/2)
  # arrays here start from 1 so need to subtract
  th=th-1; tw=tw-1;
  # upper range
  th = th+i; tw = tw+j;
  mat <- mat[(i:th), (j:tw)]
}

###---------------------###
#### Saccade detection ####
###---------------------###

# functions taken from Titus von der Malsburg
# package 'saccades' and altered slightly to 
# make it suitable for data structures used here
# https://cran.r-project.org/package=saccades


detect.saccades <- function(samples, lambda, smooth.saccades) {
  
  # Calculate horizontal and vertical velocities:
  #print(length(samples$x))
  #print(samples$trial, samples$sub_id)
  vx <- stats::filter(samples$x, -1:1/2)
  vy <- stats::filter(samples$y, -1:1/2)
  
  # We don't want NAs, as they make our life difficult later
  # on.  Therefore, fill in missing values:
  vx[1] <- vx[2]
  vy[1] <- vy[2]
  vx[length(vx)] <- vx[length(vx)-1]
  vy[length(vy)] <- vy[length(vy)-1]
  
  msdx <- sqrt(median(vx**2, na.rm=T) - median(vx, na.rm=T)**2)
  msdy <- sqrt(median(vy**2, na.rm=T) - median(vy, na.rm=T)**2)
  
  radiusx <- msdx * lambda
  radiusy <- msdy * lambda
  
  sacc <- ((vx/radiusx)**2 + (vy/radiusy)**2) > 1
  if (smooth.saccades) {
    sacc <- stats::filter(sacc, rep(1/3, 3))
    sacc <- as.logical(round(sacc))
  }
  samples$saccade <- ifelse(is.na(sacc), F, sacc)
  #samples$vx <- vx
  #samples$vy <- vy
  # drop timeseries or convert to vector
  #samples <- samples[c()]
  
  as.data.frame(samples)
  
}


detect.fixations <- function(samples, lambda=15, smooth.coordinates=T, smooth.saccades=T) {
  
  # Discard unnecessary columns:
  samples <- samples[c("x", "y", "trial", "time")]
  
  if (smooth.coordinates) {
    # Keep and reuse original first and last coordinates as they can't
    # be smoothed:
    x <- samples$x[c(1,nrow(samples))]
    y <- samples$y[c(1,nrow(samples))]
    kernel <- rep(1/3, 3)
    samples$x <- stats::filter(samples$x, kernel)
    samples$y <- stats::filter(samples$y, kernel)
    # Plug in the original values:
    samples$x[c(1,nrow(samples))] <- x
    samples$y[c(1,nrow(samples))] <- y
  }
  
  samples <- detect.saccades(samples, lambda, smooth.saccades)
  
  # not needed as this check is performed later and raw df may contain empty trials
  #if (all(!samples$saccade))
  #  stop("No saccades were detected.  Something went wrong.")

  samples
  
}

