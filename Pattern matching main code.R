# Pattern matching main code
# Author: Mark Turner 2020-2023
# 
# Seeks template-like patterns in reconstructions (DO like points, DOLPS)
# Identifies and filters warm intervals following DOLPs
# Templates are extracted from GIs in Greenland
# Uses reconstructions made by fxTWA-PLS

# Setup --------------------------------------------------
library(ggplot2) 
library(gridExtra)
library(ggrepel)
library(scales)
library(dplyr)
library(reshape2)
library(splus2R)

#****** ~~~~~~~~~~ Setup ~~~~~~~~~~~~ *****

setwd("C:/Users/markg/Dropbox/UoR/Taxon metrics/MPDB3")
# Load Kindler Greenland temp NB age (not age1) is GICC05
ref<-read.csv("C:/Users/markg/Dropbox/UoR/Taxon metrics/MPDB2/Kindler with mavg.csv")
Kind <- ref[,c("age","temp")]
# Load GICC05 GI dates etc
GI_dates<-read.csv("C:/Users/markg/Dropbox/UoR/Greenland/GI_dates_extended_v2.csv")
# Sites in longitude order
sites<-c("Villarquemado","Navarres", "Lake_Banyoles","Lac_du_Bouchet","Les_Echets_redone_2",
         "Lagaccione","Stracciacappa","Castiglione","Lago_Grande_di_Monticchio", 
         "Ioannina","Lake_Xinias","Megali_Limni", "Dead_Sea","Ghab", "Zeribar") 

#****** ~~~~~~~~~~ Template extraction ~~~~~~~~~~~~ *****
# Extracts templates for GIs 1 to 20 from Greenland Kindler temperature series
# Required only once per session; same templates work for all base reconstructions
#
# Template extraction parameters and defaults
# cki                        # Reference series e.g. Kindler 
# GI_dates                   # Dates of GI initiations GICC05modelext
lspan = 0.1                  # loess smoothing parameter 
before = 10                  # Window includes...Kindler data points older than GI date  
after  = 15                  #...younger than GI date 
safe  = 200                  # max datapoints between DOLPs to avoid multiple GIs

# Execute template extraction 
textr<-template_extract(cki = Kind,GI_dates = GI_dates)
# Output is two dataframes used by pattern_match function
tlist<-textr$tlist
nset<-textr$nset

#****** ~~~~~~~~~~ Pattern matching  ~~~~~~~~~~~~ *****
# Finds template-like points in reconstructions (DO like points, DOLPS)
# in a reconstruction of a given bioclimatic variable for a given core.
#
# Choose site, bioclimatic variable, base reconstruction version
p = 1
site = sites[[p]]
biovar = "gdd"
vtemp ="g15Jul433"                     # Base recontruction version to compare with templates;
# Load reconstruction for chosen component of this biovariable
core.sl<- read.csv(paste0("Recon fxtv1 age version for ", site, " for extended corr rtmi.csv"))
# Reconstruction contains age, depth, mean reconstructed value, and core ID
recon<-data.frame(mean = core.sl[,biovar], age = core.sl$age, core = core.sl$core, depth = core.sl$depth)

# Pattern matching parameters and defaults
vtemp = "g15Jul433"                    # specifies version of reconstruction to use
meas = "smmdist"                       # smmdist is mean of distances stdized by mean; mdist is unstandardised
lspan = 0.1                            # For loess detrending of series
before = 10; after = 15                # no of points younger and older than official date for template
tcent = TRUE                           # use y-centred templates? pref TRUE
sgn = -1                               # good Euclidean is low so -1 = find troughs not peaks
span = 101                             # How wide do you look to define troughs/peaks?
wshow = 20                             # +/- how many samples to highlight
simfilt1 = 1                           # 1 = mean stdized ED: >1 removed when finding troughs
simfilt2 = 0.9                         # for combined distances, >simfilt2 removed when finding troughs
tail = 100                             # max datapoints younger than DOLPS to consider

# Execute function
pattern_match(biovar = "gdd",vtemp = "g15Jul433")
# Output is files written out e.g. 
# paste0("testDistance based troughs ", corename," ", vtemp,tcent,".csv"))
# and paste0("testdbioa for ",corename," for ",biovar,vtemp,meas,tcent,".csv"))
# which are read in by filter_warm() function


#****** ~~~~~~~~~~ Filtering warm intervals  ~~~~~~~~~~~~ *****
# Filters the candidate DOLPS found in pattern_match() to select best
# Reads files created by pattern_match and applies filters;  
# filter settings can be varied in the function call. Result is
# shown as plot. Output is three dataframes.
#
# When the filters are deemed OK, run confirm_filtered() to
# remove DOLPS which did not make the cut.

# Filtering parameters and defaults
lspan = 0.75                                    # loess span to use if coldtest = "loess"
samfilt <- 2                                    # Removes this and lower sample counts
#                                               # Since min(sampleno) = 3, 2 removes filter 
N2filt <- 2                                     # Removes this and lower N2 intervals
tolage<-500                                     # widen window by this many years for broader measure
dfilt<-  0.75                                   # max smmdist (Euclidean distance)
rfilt<-  200                                    # min wrange (biovar range: depends on biovar)
afilt<- 1000                                    # min warea (biovar area-under-curve)
# Load Hill's N2 for all cores
coren2<-read.csv("C:/Users/markg/Dropbox/UoR/TWAPLS/Hills N2 for all cores.csv")
recondir<-"C:/Users/markg/Dropbox/UoR/Taxon metrics/MPDB3/" 

# Execute filtering function
fwarm<-filter_warm(GI_dates = GI_dates, HillN2 = HillN2, corename = sites[[p]],
                   biovar = "gdd",vtemp = "g15Jul433" )
#
N2left<-fwarm$N2left
Donly<-fwarm$Donly
dbiob<-fwarm$dbiob

#****** ~~~~~~~~~~ Confirm filtering of warm intervals  ~~~~~~~~~~~~ *****
# When the filters used in pattern_match() are deemed OK, run this to
# remove DOLPS which did not make the cut.
# This uses the three dfs made by pattern_match().

# Execute Confirm filtered result
confirm_filtered(N2left= N2left, Donly = Donly, dbiob = dbiob, vout = "g15Jul433")

# Output is three files written out
# Concatenation of vtemp and vfilt -> ensures base params AND filters used are recorded
# dbiob, paste0("testWarm patches detail for ",corename,vout,".csv"))
# dslim, paste0("testWarm patches summary for ",corename,vout,".csv"))
# filt, paste0("testWarm patches filters for ",corename,biovar,".csv"))

# Of these, dbiob is the key output used downstream.

# End
 

