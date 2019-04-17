#### Libraries ####

libs <- c("plyr", "readr", "psych", "eyetrackingR", "dplyr",
          "grid", "lattice", "tidyr") # not all are necessary

lapply(libs, library, character.only = TRUE)

# path to dir containing data
path = "/Analysis/R/"
path = "/home/maija/maija.fil@gmail.com/Research Placement/Analysis/R/"
setwd(path)

et_path = paste(path, "eye_data.txt", sep="")
et_data <- read_tsv(et_path, col_names = c("GazeX_px","GazeY_px","GazeTimestamp","L","R","trial_n",
                                           "participant","trial"),
                    col_types=cols("GazeX_px"=col_double(),
                                   "GazeY_px"=col_double(),
                                   "GazeTimestamp"=col_double(),
                                   "L"=col_double(),
                                   "R"=col_double(),
                                   "trial_n"=col_double(),
                                   "participant"=col_double(),
                                   "trial"=col_double()))#nrows=2332973)

# IDs need to match the e.t. data. Currently 1 to 30. 
exp_data <- read.csv("resp_data.txt", header = FALSE, fileEncoding="UTF-8", sep="\t")
aoi_table <- read.csv("aoi_table.csv", header=TRUE, fileEncoding="UTF-8", sep=",")
img_table <- read.table("img_dict.txt", header = FALSE, fileEncoding="UTF-8", sep="\t")
# has no headers to avoid having to open file and play around with it. should be easier to do it here
colnames(exp_data) <- c("trial_n", "sub_id", "trial", "scene", "resp", "rt")
colnames(img_table) <- c("idx", "file")
colnames(aoi_table) <- c("trial", "img_name", "left", "bottom", "right", "top")
# aoi table contains image filename for each trial and the aoi coordinates (cubestack)
# trial == image idx in global order (filenames from both folders ordered in MATLAB)


###----------------------------------------------------------------------###
############################    Preprocessing    ###########################
###----------------------------------------------------------------------###
# take first time point for each participant.
# subtract from all timepoints the first one, per participant
# Rezero each trial

# clean from zero data rows
# these mostly include erroneous data rows e.g. all zeros
# due to some issues with the eyetracker and the university machine
# missing trials are checked for later so this doesn't reflect
# on analysis
et_data <- et_data[et_data$GazeX_px!=0 & et_data$GazeY_px!=0, ]

# partitipant vector:
x <- c(et_data$participant)
(ux <- unique(x))

# participant vector:
a <- c(exp_data$sub_id)
(ax <- unique(a))

# et_data has Id values from a set containing test runs and breaks, exp_data has values 1:30 
# map participant id to eyetracking dataframe (for RT extraction)
et_data$sub_id <- mapvalues(et_data$participant, 
                            from=ux, 
                            to=ax)

for (i in ax){
  # index for condition, based on ID
  index <- et_data$sub_id == i
  # extract participant as subset
  s <- subset(et_data, et_data[,9] == i)
  # extract their trials
  n <- c(s$trial_n)
  (nx <- unique(n))
  for (j in nx){
    # index of trial in exp_data
    exp_data_idx <- exp_data$sub_id==i & exp_data$trial_n==j
    # index of 
    index <- et_data$sub_id == i & et_data$trial_n == j
    # extract a single trial
    subtrial <- subset(s, s[,6] == j)
    # get time of first gaze
    t <- subtrial$GazeTimestamp[1]
    l <- tail(subtrial, 1)
    # get RT for this trial
    rtvar <- exp_data$rt[exp_data_idx]*1000
    exp_data$rt[exp_data_idx] <- rtvar-l$GazeTimestamp
    # subtract time of first trial from all other datapoints for this trial
    et_data$GazeTimestamp[index] <- (et_data$GazeTimestamp[index] - t)
    remove(index)
  }
  # IDs 13, 16, 27, 31,52 have different values due to experiment crashing so this needs to be done for each trial
}

# check unique image values for each subject
for (i in ax){
  # extract participant as subset
  s <- subset(et_data, et_data[,9] == i)
  # extract their trials
  n <- c(s$trial_n)
  (nx <- unique(n))
  print(i)
  print(length(nx))
}

# set all unusable RT values as 9999 (7 in total)
index <- exp_data$rt < 0
exp_data$rt[index] <- 9999
index <- exp_data$rt > 2000
exp_data$rt[index] <- 9999


###----------------###
#### ADD RESPONSE ####
###----------------###

for (i in ax){
  # extract participant as subset
  s <- subset(exp_data, exp_data[,2] == i)
  # extract their trials 1:200
  n <- c(s$trial)
  (nx <- unique(n))
  print(length(nx))
  for (j in nx){
    # index for condition, based on ID
    index <- et_data$sub_id == i & et_data$trial == j
    # extract a single trial, response and scene
    subtrial <- subset(s, s[,1] == j)
    var1 <- subtrial$resp
    var2 <- subtrial$scene
    # add scene type and response to et_data
    et_data$resp[index] <- var1
    et_data$scene[index] <- var2
    
    remove(index)
  }
}

###------------------------------------###
#### DATASET WITHIN IMAGE COORDINATES ####
###------------------------------------###
# et_data with bool column corresponding to
# whether gaze is within the image, and 
# gaze coordinates corrected to image origin
# x&y offset
x_off <- 96
y_off <- 54

# monitor origin in screen coordinates (middle of screen)
scr_origin_x = 960;
scr_origin_y = 540;

img_width = 1728;
img_height = 972;

# image span in screen coordinates
x_min = scr_origin_x - img_width/2;
x_max = scr_origin_x + img_width/2;
y_min = scr_origin_y - img_height/2;
y_max = scr_origin_y + img_height/2;

et_data$image_aoi <-(et_data$GazeX_px>x_min & et_data$GazeX_px<x_max) & (et_data$GazeY_px>y_min & et_data$GazeX_px<y_max)

# set all negative values to 0 (participants looking off-screen)
index <- et_data$GazeX_px < 0
et_data$GazeX_px[index] <- 0
index <- et_data$GazeY_px < 0
et_data$GazeY_px[index] <- 0

# set all off screen coordinates to min and max values (min already done above)
index <- et_data$GazeX_px>1920
et_data$GazeX_px[index] <- 1920
index <- et_data$GazeY_px>1080
et_data$GazeY_px[index] <- 1080


# Note: Y COORDINATES NEED TO BE REVERSED TO PLOT OVER AN IMAGE
# (screen & image origin is upper left, but images are plotted knowing this)

###-------------------------------------------###
#### DATASET WITHIN BOUNDING BOX COORDINATES ####
###-------------------------------------------###

et_data = readRDS("pre-eyetracking-et-data.rds") # comment on first run

colnames(aoi_table)[1] <- "trial"
# alter the values as et_data is in screen coordinates
aoi_table$left = aoi_table$left + x_off
aoi_table$right = aoi_table$right + x_off
aoi_table$top = aoi_table$top + y_off
aoi_table$bottom = aoi_table$bottom + y_off

# Add AOI
et_data <- add_aoi(et_data, aoi_table, x_col="GazeX_px", y_col="GazeY_px", aoi_name="stack", x_min_col = "left",
                   x_max_col = "right", y_min_col = "top", y_max_col = "bottom")

#et_data <- add_aoi(et_data, aoi_table, x_col="GazeX_px", y_col="GazeY_px", aoi_name="stack", x_min_col = "left",
#                   x_max_col = "right", y_min_col = "top", y_max_col = "bottom")


###-------------------###
#### Other additions #### LR only on raw Eyetracking data, not RDS
###-------------------###
# Trackloss columns for losing one or both eyes
#et_data$trackloss <- (et_data$L == 0 & et_data$R==0)
#et_data$trackloss <- (et_data$L == 0 | et_data$R==0)

et_data<-et_data[-c(135903),] # remove duplicate row

# get rid of unnecessary columns
drops <- c("L","R","participant")
et_data<-et_data[ , !(names(et_data) %in% drops)]
#saveRDS(et_data, "pre-eyetracking-et-data-new.rds")

# make eyetracking data & treat non aoi looks as trackloss
eyetracking <- make_eyetrackingr_data(et_data, 
                                      participant_column = "sub_id",
                                      trial_column = "trial",
                                      time_column = "GazeTimestamp",
                                      trackloss_column = "trackloss",
                                      aoi_columns = c('image_aoi','stack'),
                                      treat_non_aoi_looks_as_missing = TRUE)


#saveRDS(exp_data, "exp_data-rt.rds")
#saveRDS(eyetracking, "eyetracking-new.rds")

###----------------------###
#### Trackloss Analysis ####
###----------------------###

# Check how much data is lost after cleaning
# & whether the amount left is reasonable 
# Trackloss not significant, though some trials will
# need to be removed due to lack of samples still
# threshold <25% of the trial (<1.75 seconds)

#eyetracking = readRDS("eyetracking.rds") 

(trackloss <- trackloss_analysis(data = eyetracking))
eyetracking_clean <- clean_by_trackloss(data = eyetracking, trial_prop_thresh = .25)
# removes 30 Trials
eyetracking_clean$Target <- (eyetracking_clean$resp==eyetracking_clean$scene)

trackloss_clean <- trackloss_analysis(data = eyetracking_clean)
(trackloss_clean_subjects <- unique(trackloss_clean[, c('sub_id','TracklossForParticipant')]))

# get mean samples contributed per trials, with SD
mean(1 - trackloss_clean_subjects$TracklossForParticipant)
sd(1- trackloss_clean_subjects$TracklossForParticipant)

# contribution per participant
(final_summary <- describe_data(eyetracking_clean, 'stack', 'sub_id'))

# save new eyetracking RDS
#saveRDS(eyetracking_clean, "eyetracking-new.rds") 
# out of bounds coordinates set to max/min

###----------------------------------###
#### Extract time window if desired ####
###----------------------------------###

# subset by window, 7 seconds
response_window <- subset_by_window(eyetracking_clean,
                                    window_start_time = 0,
                                    window_end_time = 7000, rezero = FALSE, remove = TRUE)

# subset by window: take last 7000 seconds from each trial
for (i in ax){
  # index for condition, based on ID
  index <- eyetracking$sub_id == i
  # extract participant as subset
  s <- eyetracking[eyetracking$sub_id==i, ]
  # extract their trials
  n <- c(s$trial_n)
  (nx <- unique(n))
  for (j in nx){
    # extract a single trial
    subtrial <- eyetracking[eyetracking$trial_n==j, ]
    # get time of first gaze
    l <- tail(subtrial, 1)
    # start of trial (if trial longer than 7k ms)
    timevar <- l$GazeTimestamp-7000
    # true if timestamp larger than start time
    idx_T <- eyetracking$sub_id == i & eyetracking$trial_n == j & eyetracking$GazeTimestamp>=timevar
    idx_F <- eyetracking$sub_id == i & eyetracking$trial_n == j & eyetracking$GazeTimestamp<timevar
    eyetracking$window[idx_T] <- TRUE
    eyetracking$window[idx_F] <- FALSE
    
    remove(index)
  }
}


###-------------------------###
#### Translate Coordinates ####
###-------------------------###

# translate coordinates to image frame
# this will create negative values, so save image only coordiantes
scr_only <- eyetracking # for saving - this has positive coords
eyetracking$GazeX_px <- (eyetracking$GazeX_px- x_off)
eyetracking$GazeY_px <- (eyetracking$GazeY_px- y_off)

# subset by constraining to image
img_only <- eyetracking[eyetracking$GazeX_px>=0 & eyetracking$GazeY_px>=0, ]
img_only <- eyetracking[eyetracking$GazeX_px<=1728 & eyetracking$GazeY_px<=972, ]

# save RDS
#saveRDS(scr_only, "eyetracking_scr_only.rds")
#saveRDS(img_only, "eyetracking_img_only_new.rds")
# use scr only for extracting saccades as img only loses information about movements

###-----------------------------###
#### Scale to small image size ####
###-----------------------------###

et_data$GazeX_px <- (et_data$GazeX_px/4)
et_data$GazeY_px <- (et_data$GazeY_px/4)
