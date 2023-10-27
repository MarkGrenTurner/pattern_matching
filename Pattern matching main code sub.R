# Pattern matching main code sub
# Author: Mark Turner 2020-2023
# 
# Extracts templates representing Dansgaard-Oeschger events from Greenland data
# Seeks template-like patterns in quantitative climate reconstructions (D-O like points, DOLPS)
# Identifies and filters warm intervals (possible interstadials) following DOLPs
# 
# Uses gdd (growing degree day) reconstructions made by fxTWA-PLS.
#
# -------- Local setup --------------------------------------------------
library(ggplot2) 
library(gridExtra)
library(ggrepel)
library(scales)
library(dplyr)
library(reshape2)
library(splus2R)

#****** ~~~~~~~~~~~~~~~~~ File/variable descriptions~~~~~~~~~~~~~~ *****

# Kind        Kindler temperature series with columns c("age","temp")
# GI_dates    start dates of GIs with columns c("type","event","ageb2K","uncertainty_1s")
# recon       quantitative reconstruction for one core, at physical sample level
#                 for one component for one bioclimatic variable
#                 with columns c("mean","age","core","depth") where mean = reconstructed value
# coren2      Hill's N2 for fossil core interpolated to Kindler resolution
#                 with columns c("core","depth","age","fossN2")
# vtemp       a string used to identify versions of input files
# vout        a string used to identify versions of output files

#****** ~~~~~~~~Local loading of files and naming as needed ~~~~~~~~ *****

setwd("C:/Users/markg/Dropbox/UoR/Taxon metrics/MPDB3")
# Load Kindler Greenland temp NB age (not age1) is GICC05
ref<-read.csv("C:/Users/markg/Dropbox/UoR/Taxon metrics/MPDB2/Kindler with mavg.csv")
Kind <- ref[,c("age","temp")]
# Load GICC05 GI dates etc
GI_dates<-read.csv("C:/Users/markg/Dropbox/UoR/Greenland/GI_dates_extended_v2.csv")
# Load Hill's N2 for all cores
coren2<-read.csv("C:/Users/markg/Dropbox/UoR/TWAPLS/Hills N2 for all cores.csv")
# Sites in longitude order, e.g.
sites<-c("Villarquemado","Navarres", "Lake_Banyoles","Lac_du_Bouchet","Les_Echets_redone_2",
         "Lagaccione","Stracciacappa","Castiglione","Lago_Grande_di_Monticchio", 
         "Ioannina","Lake_Xinias","Megali_Limni", "Dead_Sea","Ghab", "Zeribar") 
# Select site and biovar
p = 1
site = sites[[p]]
biovar = "gdd"            # bioclimatic variable nme
vtemp ="g15Jul433"        # Base reconstruction version to compare with templates
# Load reconstruction for chosen component of this biovariable
core.sl<- read.csv(paste0("Recon fxtv1 age version for ", site, " for extended corr rtmi.csv"))
# Reconstruction contains age, depth, mean reconstructed value, and core ID
recon<-data.frame(mean = core.sl[,biovar], age = core.sl$age, core = core.sl$core, depth = core.sl$depth)


#****** ~~~~~~~~~~ Template extraction ~~~~~~~~~~~~ *****
# Extracts templates for GIs 1 to 20 from Greenland Kindler temperature series
# Required only once per session; same templates work for reconstructions of any climatic value
#
# Template extraction parameters and defaults
# Requires
#   cki                      # Reference series e.g. Kindler or d18O
#   GI_dates                 # Dates of GI initiations on GICC05modelext chronology
lspan = 0.1                  # loess smoothing parameter for detrending
before = 10                  # Window includes...Kindler data points younger than GI date  
after  = 15                  # ...older than GI date 
safe  = 200                  # max data points between DOLPs to avoid multiple GIs

# Execute template extraction 
textr<-template_extract(cki = Kind, GI_dates = GI_dates)
# Output is two dataframes used by pattern_match function
tlist<-textr$tlist
nset<-textr$nset

#****** ~~~~~~~~~~ Pattern matching  ~~~~~~~~~~~~ *****
# Finds template-like points in reconstructions (DO like points, DOLPS)
# in a reconstruction of one bioclimatic variable for one core.
#
# Pattern matching parameters and defaults
vtemp = "g15Jul433"                    # specifies version of reconstruction to use
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
# Output is files written out 
#     paste0("testDistance based troughs ", corename," ", vtemp,tcent,".csv"))
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
samfilt <- 2                                    # Removes this and lower sample counts
#                                               # Since min(sampleno) = 3, 2 removes filter 
N2filt <- 2                                     # Removes this and lower N2 intervals
tolage<-500                                     # widen window by this many years for broader measure
dfilt<-  0.75                                   # max smmdist (Euclidean distance)
rfilt<-  200                                    # min wrange (biovar range: depends on biovar)
afilt<- 1000                                    # min warea (biovar area-under-curve)

# Load output from pattern_match
dbioa<-read.csv(paste0("testdbioa for ",corename," for ",biovar,vtemp,meas,tcent,".csv"))

# Execute filtering function
fwarm<-filter_warm(GI_dates = GI_dates, recon = recon, dbioa = dbioa, HillN2 = coren2, corename = sites[[p]],
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

# Of these, dbiob is the key output for use downstream and contains much 
# information collected en route.

#****** Function definitions ~~~~~~~~~~~~~~~~~~~~~~~~ *****

#****** Template extraction function ~~~~~~~~~~~~~~~~ *****
# Author: Mark Turner 2020-2023
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extracts templates for GIs 1 to 20 from Greenland Kindler temperature series
# Required only once per session; same templates work for all base reconstructions

template_extract <- function(cki, 
                             GI_dates, 
                             lspan = 0.1, 
                             before = 10, 
                             after = 15, 
                             safe = 200) {
  cvs <-GI_dates[GI_dates$type =="GI start",]              # GI starts only
  
  # Get GI ages and names excl GI 0 and > 20
  intval<-round(mean(diff(cki$age)),0)                    # age resolution of Kindler
  gnset<-cvs[cvs$event>0 & cvs$event<21,"ageb2K"]         # actual GI initiation ages (= cva): midpoint of rise
  nnset<-cvs[cvs$event>0 & cvs$event<21,"event"]          # names of GIs
  # Detrend and normalise Kindler before extraction ("pretreat")
  reference<- cki$temp                                     # reference is raw Kindler
  ages<-cki$age                                            # ages = age data in age-stdised core/biovar recon
  ref.norm<-scale(reference)                               # z score reference series
  sdf<-predict(loess(ref.norm~cki$age, span = lspan))      # smoothed "trend" to be extracted
  ref.nd<-ref.norm-sdf                                     # detrended normalised reference
  jsh<-as.data.frame(cbind(ages, ref.nd))                  # detrended normalised Kindler temp with ages
  colnames(jsh)<-c("age","temp")
  
  # Template extraction
  
  run<-seq(1:(before+after+1))                      # running "age" within template
  # Set up df for template using first GI
  this<-which.min(abs(ages - gnset[1]))             # nearest row to 'official' date of warming
  tlist<-jsh[(max(this-before,1)):(this+after),]    # get patch of detr, norm Kindler for this interval
  tlist$GI<-nnset[1]                                # name it
  tlist$warm<-NA
  tlist$rise<-NA
  # Find start/end of rise
  tlist[which(diff(tlist[,"temp"]) <0),"rise"]<-1   # label intervals where value rises
  tlist[c(1:(before-7), (before+8):(before+after+1)),"rise"]<-NA #.. within limits
  # Locate important point: Find range of rise
  top<-max(tlist[tlist$rise == 1, "temp"], na.rm = TRUE)
  foot<-min(tlist[tlist$rise == 1, "temp"], na.rm = TRUE)
  midpoint<-which.min(abs(tlist[tlist$rise == 1, "temp"]-((top-foot)/2+foot))) # finds row nearest to midpoint
  tlist[midpoint,"warm"]<-tlist[midpoint,"age"]
  # Running "age" within template: handle case where earliest available data younger than "before"
  if(before>this){tlist$run<-run[(before-this+2):(before+after+1)]
  }else{
    tlist$run<-run                                    
  }
  tlist$centred<-tlist$temp-mean(tlist$temp)        # Centre the data by subtracting mean
  pin<-cki[this,"age"]+safe                         # Set pin: "safe" distance older than this warming   
  # Extract remaining templates into tlist
  for (i in 2:length(gnset)){
    this<-which.min(abs(ages - gnset[i]))           # find nearest row to given warming
    tr<-jsh[(max(this-before,1)):(this+after),]     # get the defined age interval around this
    tr$GI<-nnset[i]                                 # name of GI
    tr$warm<-NA
    tr$rise<-NA
    # Find start/end of rise
    tr[which(diff(tr[,"temp"]) <0),"rise"]<-1       # all intervals where the value increases
    tr[c(1:(before-7), (before+8):(before+after+1)),"rise"]<-NA #.. only within limits (set for "centre")
    # Find range of rise
    top<-max(tr[tr$rise == 1, "temp"], na.rm = TRUE)
    foot<-min(tr[tr$rise == 1, "temp"], na.rm = TRUE)
    midpoint<-which.min(abs(tr[tr$rise == 1, "temp"]-((top-foot)/2+foot))) # finds row nearest to midpoint
    tr[midpoint,"warm"]<-tr[midpoint,"age"]
    if(before>this){tr$run<-run[(before-this+2):(before+after+1)]
    }else{
      tr$run<-run                                    
    }
    # Test for overlap with previous GI: do not get within "safe" distance of previous warming
    for (j in 1: nrow(tr)){                         # read along extract
      if(tr[j,"age"]<=pin){tr[j,"temp"]<-NA}        # temp = NA if age younger than previous pin
      # Special cases trimmed to avoid 2 GIs in 1 template
      if(tr[j,"GI"] == "14"   & tr[j,"run"]>125){tr[j,"temp"]<-NA}
      if(tr[j,"GI"] == "15.1" & tr[j,"run"]>115){tr[j,"temp"]<-NA}
      if(tr[j,"GI"] == "16"   & tr[j,"run"]>110){tr[j,"temp"]<-NA}
    }
    # Centre the data (subtract mean)
    tr$centred<-tr$temp-mean(tr$temp, na.rm = T)    # centred temperature
    tlist<-rbind(tlist, tr)                         # add this template (tr) to tlist
    pin<-cki[before,"age"]+safe                     # reset pin for this GI
  }
  tlist$rage<-NA                                    # Used for x axis for vline in plot
  for (k in 1:nrow(tlist)) {
    if(!is.na(tlist[k,"warm"])) tlist[k,"rage"]<-tlist[k,"run"]
  }
  #
  nset<-as.list(unique(tlist$GI))                    # Names of GI entries
  cvg<-cvs[which(cvs$event %in% unique(tlist$GI)),]  # Parallel list of cvg entries
  gset<-as.list(cvg$ageb2K)                          # Warming ages
  
  # return template details
  list <- list(tlist, nset)
  names(list)<-c("tlist", "nset")
  return(list)
}

#****** Pattern matching function ~~~~~~~~~~~~ *****
# Mark Turner 2020-2023
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Takes templates of D-O initiations provided by 'template_extract' function,
# finds Euclidean distance between templates and quantitative reconstruction,
# combines effect of multiple templates, and finds low points in the combined 
# ED curve which represent most D-O-like points (DOLP).

# Requires 'tlist' and 'nset' data frames from 'template_extract' function
#
# Outputs two files for use by 'filter_warm' function. Internally 
# these dfs are named 'tronly' and 'dbioa'.
#
# ~~~~~~~~~~~~~~~~~~ Meanings and typical values of variables ~~~~~~~~~~~~~
# vtemp = "g15Jul433"            # specifies version of reconstruction to use
# meas = "smmdist"               # smmdist is mean of distances stdized by mean; mdist is unstandardised
# lspan = 0.1                    # For loess detrending of series
# before = 10; after = 15        # no of points younger and older than official date for template
# tcent = TRUE                   # use y-centred templates? pref TRUE
# sgn = -1                       # good Euclidean is low so -1 = find troughs not peaks
# span = 101                     # How wide do you look to define troughs/peaks?
# wshow = 20                     # +/- how many samples to highlight
# simfilt1 = 1                   # 1 = mean stdized ED: >1 removed when finding troughs
# simfilt2 = 0.9                 # for combined distances, >simfilt2 removed when finding troughs
# tail = 100                     # max datapoints younger than DOLPS to consider
#
pattern_match <- function(site = sites[[p]],
                          core = recon,
                          cki = Kind,
                          biovar = biovar,
                          vtemp = vtemp,                       
                          meas = "smmdist",     
                          lspan = 0.1,                                              
                          before = 10,
                          after = 15,                                 
                          tcent = TRUE,                                            
                          sgn = -1,                                                
                          span = 101,                                               
                          wshow = 20,                                              
                          simfilt1 =1,                                             
                          simfilt2 = 0.9,                                           
                          tail = 100) {                   
  
  this<-site; corename<-this; print(corename)
  
  ### Make original core df for this biovar
  core$core<-as.character(core$core)
  corename<-core$core[1]; print(corename)
  coret<-as.data.frame(cbind(depth = core$depth,recm = core$mean,age = core$age))
  # coret is core raw data
  ### Find age limits and truncate
  print(paste("Core min age = ",min(coret$age)," and max =", max(coret$age)))
  minage<-min(coret$age); maxage<-max(coret$age)
  coret<-coret[coret$age>=minage&coret$age<=maxage,] # core data, age-truncated 
  coret<-cbind(coret,row=seq(1:nrow(coret)))       # with row numbers; coret is core raw data
  c3<-cki[cki$age>minage&cki$age<maxage,]          # Kindler temp data, age-truncated to match core
  kinb<-cbind(c3, row=seq(1:nrow(c3)))             # with row numbers: kinb is Kindler base (raw) data
  ### Standardise x axis for both
  stdx<-kinb$age                                                       # Standardise on Kindler axis
  corest<-as.data.frame(approx(coret$age, coret$recm, xout = stdx))    # interpolate core biovar on std axis
  # Interpolate depths
  if(all(!is.na(coret$depth))){                    # if there's depth given...
    corest1<-as.data.frame(approx(coret$age, coret$depth, xout = stdx))# interpolate core depth on std axis
    depint<-corest1$y
  }else{
    depint<-rep(NA,times = length(stdx))           # if no depth given set to NA
  }
  #
  colnames(corest)<-c("age","recm")
  corest<-cbind(corest, depth = depint, row=seq(1:nrow(corest)))       # corest is core data interpolated to std ages
  ###
  ctemp<-corest[,c(1:3)]                           # part of current age-standardised core/biovar recon
  colnames(ctemp)<-c("age","temp","depth")         # "temp" is climate variable, not necc temperature
  c2<-ctemp                                        # c2 is age-stdised core/biovar recon
  c2$yes<-0                                        # set up field
  reference<- c2$temp                              # reference = core recon, where you look for DOs
  n.reference = length(reference)
  ref.norm<-scale(reference)                       # z score whole ref series
  sdf<-predict(loess(ref.norm~c2$age, span = lspan))   # smoothed "trend" to be extracted (e.g. 0.1)
  # Validity of loess was tested with SSA: no significant difference
  ref.nd<-ref.norm-sdf                             # detrended normalised reference
  
  #*******************************************************************************
  # Inner loop comparing all templates with reconstruction for this core/biovar
  #*******************************************************************************
  for(k in 1:length(nset)){                            # No of GI templates
    if(tcent){                                         # use template centred on mean
      cset<-tlist[tlist$GI==nset[k],c("age","centred","rage")] # Current template centred: interval for this GI
    }else{
      cset<-tlist[tlist$GI==nset[k],c("age","temp","rage")]   # Current template: interval for this GI
    }
    colnames(cset)<-c("age","temp","rage")
    keypt<-unique(cset$rage) [2]                       # find row no of template regarded as 'key'
    query <- cset$temp                                 # query is template
    ages<-c2$age                                       # ages = age data in age-stdised core/biovar recon
    n.query = length(query)                            # how long is the template
    rpoint<-n.query-keypt                              # how far back from end of template do we record ED
    print(paste("Template  = ",nset[k],"template length =",n.query))          
    #********************************************************************************
    # Compute Euclidean distances by sliding template along reconstruction for this core/biovar
    #********************************************************************************
    dist = rep(NA, n.reference)                        # empty vector for EDs
    for(i in n.query:n.reference) {                    # start/end distances along which to slide window
      window.rnd = ref.nd[(i-n.query+1):i]             # window from normalised detrended full series
      diff<-rep(NA,n.query)                            # empty vector for difference in window
      diff<-query-window.rnd                           # diff between window and template
      diff<-diff/n.query                               # stdised for template length
      dist[i-rpoint]<-sqrt(sum(diff^2, na.rm = T))     # dist = Euclidean distance 
    }                   # dist is scalar stdised ED, recorded at *warming point* of sliding window 
    dmean<-mean(dist, na.rm = T)                       # mean of ED across all templates
    ### Preserve full set of distance data for this core/biovar/GI
    alldist<-cbind(c2,dist,dist/dmean)                 # "smdist" = ED/mean ED, so normalised
    colnames(alldist)<-c("age","temp","depth", "yes","dist","smdist")
    for(g in 1:nrow(alldist)) {                        # remove dist where there are NA e.g Navarres 
      if(is.na(alldist[g,"temp"])){
        alldist[g,"dist"]<-NA
      } 
    }
    write.csv(alldist,paste0("alldist",corename,vtemp,biovar,nset[[k]],tcent,".csv"))
    ### Find minima (most DO-like points) in ED for this core/biovar/GI
    filtdist<-tdist<-alldist$smdist                    # smdist standardises
    filtdist[filtdist>simfilt1]= NA                    # replace EDs above mean with NA
    try(pks<-splus2R::peaks(sgn*filtdist, span = span))# sign controls peaks v troughs; span = how wide you look to define trough
    min.index<-which(pks)                              # vector of indices in filtdist where troughs are found
    n.found = length(min.index)                        # how many troughs were found
    ### Flag the parts of the recon which match
    c2$yes<-0
    for(i in 1:n.found) {
      for(j in ((min.index[i]-wshow):(min.index[i]+wshow))) {  # fixed width to show
        c2[j,"yes"]<-1                                   # yes = 1 means this is highlightable
      }
    }
    
    ### Get matched points
    lkj<-as.data.frame(cbind(age = alldist$age[min.index], 
                             depth = alldist$depth[min.index], # new addition
                             pts =alldist$smdist[min.index]))  # get matched points
    lkj<-unique(lkj)
    hjk<-merge(c2,lkj,by = "age",all.x = T)
    hjk<-hjk[c(1:nrow(c2)),]                                   # removes repeating group at end of file
    ### Save peaks in core/biovar/GI
    write.csv(hjk,paste0("DOs like ",nset[[k]], vtemp,tcent," in ",corename,biovar,".csv")) 
    
    # *******************************************************
  } #  End of main inner all-template loop for single biovar
  # *******************************************************
  # Assemble core/biovar results --------------------------------------------
  
  # Assemble all template results for this core/biovar # rage is row of 'key' point in template
  for (s in 1:length(nset)) {                # no of templates
    pt<-read.csv(paste0("DOs like ",nset[[s]],vtemp,tcent," in ",corename,biovar,".csv"))
    pt$GI<-nset[[s]]; pt$rage<-pt$age
    alld<-read.csv(paste0("alldist",corename,vtemp,biovar,nset[[s]],tcent,".csv"))  # Continuous distance
    alld$GI <- nset[[s]]                     # Label with GI; record distance at warming age (centre of short)
    if(s == 1){
      ptsd<-pt
      alldists<-alld
    }else{
      ptsd<-rbind(ptsd,pt)
      alldists<-rbind(alldists, alld)
    }
  }
  #
  ptsdf<-ptsd[!is.na(ptsd$pts),]             # Strip out all but matched points
  # ptsdf is core/biovar matched points for all templates shown separately
  alldists$core<-corename
  # alldists is core/biovar distances for all templates shown separately
  # *************************************************************
  #  Summarising values for this core/biovar across all templates
  # *************************************************************
  # Calculate measures of combined distance across templates for this core/biovar
  alldists%>% group_by(age)%>% summarize(smmdist = mean(smdist, na.rm = T), # mean of mean-stdised distances across all templates
                                         mdist = mean(dist, na.rm = T),     # mean of raw ED across all templates
                                         temp = mean(temp, na.rm = T),
                                         depth = mean(depth, na.rm = T))->dbio 
  dbio<-as.data.frame(dbio)
  dbio<-dbio[!is.na(dbio$mdist),]                         # removes NaN, Inf, -Inf
  stage<-min(dbio$age);fage<-max(dbio$age)                # need to truncate corest similarly
  corest<-corest[ages>=stage & ages<=fage,]
  ### Find lowest distances
  filtdist2<-dbio$smmdist                                  
  ###
  filtdist2[filtdist2>simfilt2] = NA                      # replace distances above set point with NA
  
  # Locate warm parts of series --------------------------------------------
  simfilt3<-simfilt2
  filtdist2a<-filtdist2                     
  filtdist2a[filtdist2a>simfilt3]= NA                     # replace distances above set point with NA
  try(dpksa<-splus2R::peaks(sgn*filtdist2a, span = span)) # note -ve; span = how wide you look to define trough
  mina.index<-which(dpksa)                                # vector of indices where troughs are found
  m.found = length(mina.index)                            # how many troughs were found
  #
  dbioa<-dbio
  dbioa$DOLP<-FALSE
  
  # Flag the parts of dbioa (which has all ages) which form GI (based on interpolated samples)
  # DO warm part defined as recon >= that at DOLP, within 'tail' 
  dbioa$warm<-0                                           # flag for warm patch ('GI')
  dbioa$wrange<-NA                                        # value of DOLP-to-peak range
  DOno<-1
  for(z in mina.index) {                                  # for each DOLP
    dbioa[z,"DOLP"] = TRUE                                # flag this row as a DOLP
    dbioa[z, "DOno"]<-DOno                                # number them sequentially
    DOno<-DOno+1
  }
  #
  for(i in 1:m.found) {                                   # for each DOLP
    j = k = mina.index[i]                                 # find DOLPoint
    rg<-0                                                 # reset range for this DOLP
    while(j > 2                                           # don't go past beginning of file
          & dbioa[j,"temp"]>=dbioa[k,"temp"]              # for all points warmer than DOLP
          & dbioa[j-1,"DOLP"] == FALSE                    # but not running into next DOLP 
          & (k-j)<tail) {                                 # and only for a max of 'tail' samples
      dbioa[j,"warm"]<-1                                  # flag as warmer than at DOLP
      tdiff<- dbioa[j,"temp"]-dbioa[k,"temp"]             # t difference from DOLP
      if(tdiff > rg) rg<-tdiff                            # build up max difference
      j = j-1                                             # work towards younger
    }
    m <- k                                                # start from DOLP again
    while(dbioa[m,"warm"] == 1) {                         # if a warm patch, show range
      dbioa[m,"wrange"]<-rg
      m = m-1
    }
  }
  dbioa$core<-corename
  dbioa$biovar<-biovar
  # dbioa is warm intervals
  
  # Find just the rows with troughs in dbioa
  tronly<-dbioa[dbioa$DOLP == TRUE,]                             
  # Write file for use by "Find nearest matches..." and "Trial of of GIs.."
  write.csv(tronly,paste0("testDistance based troughs ", corename," ", vtemp,tcent,".csv"))
  write.csv(dbioa, paste0("testdbioa for ",corename," for ",biovar,vtemp,meas,tcent,".csv"))
  
  # end of finding warm patches
} 

#****** Filtering warm intervals functions *****
# Mark Turner 2020-2023
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Takes 'dbioa' Euclidean distance df from 'pattern_match' function,
# finds subsequent warm patches. Filters for ED, rise, area-under-curve, Hill's N2,
# and no of physical samples to propose acceptance/rejection of DOLP.
#
# This version is designed to use gdd as the base reconstruction. Other bases
# would need code modifications as well as different filter values. 
#
# When filter settings provide acceptable answers based on review of plot,
# run 'confirm_filtered' function. 
#
# Requires 'dbioa' dataframe from 'pattern_match' function
#
# Returns dataframes 'N2left', 'Donly',' and 'dbiob'

# ~~~~~~~~~~ Variable meanings: see 'pattern_match' ~~~~~~~~~~~

filter_warm<-function(corename,
                      recon,
                      dbioa,
                      GI_dates,              # for plotting only
                      HillN2, 
                      biovar,
                      vtemp,
                      meas = "smmdist",
                      tcent = TRUE,
                      limt =   20,           # fewest no of sequential cold points permitted
                      samfilt = 2,           # Removes this and lower sample counts
                      N2filt = 2,            # Removes this and lower N2 intervals
                      tolage = 500,          # widen window by this many years for broader measure
                      dfilt =  0.7,          # max smmdist (Euclidean distance) for DOLP
                      rfilt =  200,          # min wrange (biovar range: depends on biovar)
                      afilt =  1000) {       # min warea (biovar area-under-curve)
  
  # Find Hill's N2 for core
  cvs <-GI_dates[GI_dates$type =="GI start",] # GI starts only
  N2core<-HillN2[HillN2$core == corename,]
  # Rename sample-based reconstruction
  cgdd<-xcor_rtmi<-recon
  
  # Warm patch stuff ----------------------------------------------------
  
  dbiob<-dbioa[dbioa$age >13500,]             # exclude termination
  dbiob$OlDOLP<-dbiob$DOLP                    # preserve original DOLPs
  # Find start and end points chronologically in interpolated base variable of 'GI' patches
  dbiob$stend<-NA
  for(i in 1:(nrow(dbiob)-1)){
    if(dbiob[i+1,"warm"]== 1& dbiob[i,"warm"]== 0 ) dbiob[i+1,"stend"]<-"end"                                       # going backwards in time
    if(dbiob[i,"DOLP"] == TRUE) dbiob[i,"stend"]<-"start"
  }
  # Interpolate xcor_rtmi (sample based) onto interpolated axis as in dbioa
  stdx<-dbiob$age                                                       
  # Interpolate gdd (sample based) on std axis
  gdds<-as.data.frame(approx(cgdd$age, cgdd$mean, xout = stdx))
  colnames(gdds)<-c("age","gdd")
  dbiob<-inner_join(dbiob,gdds)             # tack on interpolated 
  
  # Work out warm interval values: area-under-curve, rise (range)
  dbiob$warea <- NA                         # base variable in warm area
  dbiob$gbump <- NA                         # the base variable bit above the DOLP
  
  for (k in 1:nrow(dbiob)){               
    if(dbiob[k,"DOLP"] == TRUE) {           # then k is a DOLP row
      # Work backwards from k to sum difference, etc
      inx = k
      star = 0   
      dtp = dbiob[k,"temp"]                          # base variable at DOLP
      while (dbiob[inx, "warm"] == 1 & inx>1) {      # stops at "end" of warm patch
        star<-star+dbiob[inx, "temp"]-dtp            # sum of base variable since DOLP: natural +ve
        dbiob[inx, "gbump"]<-dbiob[inx, "temp"]-dtp  # this is the excess over the DOLP
        inx <-inx-1
      }
      dbiob[k, "warea"]<-star                   # record base variable area-under-curve
    }
  }
  # Find number of original samples in warm intervals
  # Rationale: DOLP defining start is midpoint of rise, and end is point at which biovar declines below that 
  # at DOLP, neither of which are physical samples. So no of points which define warm interval shape include an
  # extra point both before and after the warm interval. Hence the +2.
  # Also, finds min(Hill's N2) among samples defining the warm interval
  dbiob$nsample1<-dbiob$nsample2<-dbiob$N2<-NA  # set up new cols
  warmbits<-dbiob[!is.na(dbiob$stend),]         # start and end of warm intervals already tagged
  warmbits$N2<-NA                               # original samples N2 minimum
  if(nrow(warmbits) == 1 ) {                    # | warmbits[1, "stend" ] == "start"
    warmbits[1,"nsample1"]<-warmbits[1,"nsample2"]<-NA # can't know - one line only
  }else{
    for (k in 2:nrow(warmbits)) {
      if(warmbits[k,"DOLP"] == TRUE) {          # start of warm patch
        ageold<-warmbits[k,"age"]               # age of start
        if(warmbits[k-1,"stend"] == "end") {    # end of warm patch
          agenew<-warmbits[k-1,"age"]           # age of end
        }
        # insert sample counts and N2 at start/DOLP row
        nsam1<-nrow(xcor_rtmi[xcor_rtmi$age >= agenew & xcor_rtmi$age <= ageold,])+2 # See above for why +2
        dbiob[dbiob$age == ageold,"nsample1"]<-nsam1
        nsam2<-nrow(xcor_rtmi[xcor_rtmi$age >= agenew-tolage & xcor_rtmi$age <= ageold+tolage,])
        dbiob[dbiob$age == ageold,"nsample2"]<-nsam2
        N2<-min(N2core[N2core$age >= agenew & N2core$age <= ageold, "fossN2"])
        dbiob[dbiob$age == ageold,"N2"]<-N2                              # original samples N2 minimum
      }
    }
  }
  dslim<-dbiob[dbiob$DOLP == TRUE , ]          # initial set of potential DOLPs
  
  
  # Filter warm patches -----------------------------------------------------
  
  # Use trial filters and find DOLPS which don't make the cut
  # No of official GIS in range, ignoring anything before 13.5 ka (avoid termination)
  GIs<-cvs[cvs$ageb2K >= max(min(dbioa$age),13500) & cvs$ageb2K <=max(dbioa$age),]
  ngiir<-paste("No of GIs in range",nrow(GIs))
  print(ngiir)
  # All candidates
  Donly<-dbiob[dbiob$DOLP == TRUE, ]
  ndcand<-paste("No of DOLP candidates", nrow(Donly))
  print(ndcand)
  # Ones to keep, based on filters 
  tryleft<-Donly[Donly$smmdist   < dfilt 
                 & Donly$wrange  > rfilt 
                 & Donly$warea   > afilt, ]
  ndlafter<-paste("No of DOLPs left after ED/area/range filtering",nrow(tryleft))
  print(ndlafter)
  # vector of DOLP numbers to keep
  Dkeep<-tryleft$DOno 
  # Further eliminating DOs for low sample count
  samfilt1<-samfilt
  #
  samleft<-tryleft[tryleft$nsample1 > samfilt1,]   
  Dkeep<-samleft$DOno 
  ndpostsam<-paste("No of DOLPs left after sample filter",nrow(samleft))
  print(ndpostsam)
  
  # Further eliminating DOs for min Hill's N2 below limit
  N2left<-samleft[samleft$N2 > N2filt,]
  Dkeep<-N2left$DOno 
  ndpostN2<-paste("No of DOLPs left after N2 filter",nrow(N2left))
  print(ndpostN2)
  
  # Reconstruction and DOLPs
  p1<-ggplot(dbiob)+
    geom_line(data = dbioa,aes(age, temp))+                                  # interpolated
    geom_line(aes(age, ifelse(warm == 1, temp, NA)), col = "red", size = 1)+ # ...with all warm patches
    geom_point(aes(age,ifelse(DOLP,temp,NA)), col = "#0072B2", size = 2)+    # all DOLPs                                               # All original DOLPS
    geom_point(data = N2left, aes(age,ifelse(DOLP,temp,NA)), col = "#D55E00", size = 2)+# surviving DOLPs 
    ylab(biovar)+                                      # DOLPS surviving filter
    annotate("text", x = 15000, y = 1.1*min(dbiob$temp,na.rm = T), label = corename)+
    xlim(c(10000,max(dbiob$age)))+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  
  # ED test
  g0<-ggplot(dslim)+
    geom_point(aes(age,smmdist), size = (1-dslim$smmdist)*5, col = "orange")+# orange will be binned on smmdist
    geom_line(data = dbiob, aes(age,smmdist), col = "orange")+# orange will be binned on smmdist
    #geom_point(aes(age,ifelse(dslim$smmdist<=dfilt,dslim$wrange,NA)), 
    #           size = (1-dslim$smmdist)*10)+  
    geom_hline(aes(yintercept = dfilt), col = "red")+
    xlim(c(10000,max(dbiob$age)))+
    ylim(c(0,2))+
    labs(y = "ED")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  # Range test
  g1<-ggplot(dslim)+
    geom_point(aes(age,wrange), size = (1-dslim$smmdist)*5, col = "orange")+# orange will be binned on smmdist
    geom_point(aes(age,ifelse(dslim$smmdist<=dfilt,dslim$wrange,NA)), 
               size = (1-dslim$smmdist)*10)+  
    geom_hline(aes(yintercept = rfilt), col = "red")+
    xlim(c(10000,max(dbiob$age)))+    
    labs(y = "Rise")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  # Area test
  g2<-ggplot(dslim)+
    geom_point(aes(age,warea), size = (1-dslim$smmdist)*5, col = "orange")+# orange will be binned on smmdist
    geom_point(aes(age,ifelse(dslim$smmdist<=dfilt,warea,NA)), 
               size = (1-dslim$smmdist)*10, col = "red")+      
    geom_hline(aes(yintercept = afilt), col = "blue")+
    xlim(c(10000,max(dbiob$age)))+
    labs(y = "Area")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  # Sample count test
  g3<-ggplot(dslim)+
    geom_point(aes(age,nsample1), col = "#56B4E9", size = 2)+ # orange will be binned on smmdist
    xlim(c(10000,max(dbiob$age)))+
    ylim(2,max(dslim$nsample1))+
    labs(y = "Sample count")+
    theme_bw()+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  if(samfilt>2) g3<-g3+geom_hline(aes(yintercept = samfilt), col = "#D55E00")
  
  # N2 test
  g4<-ggplot(dslim)+
    geom_point(aes(age,N2), col = "purple", size = 2)+ 
    geom_hline(aes(yintercept = N2filt), col = "red")+
    xlim(c(10000,max(dbiob$age)))+
    ylim(1,min(5,max(dslim$N2)))+
    labs(y = "N.2")+
    theme_bw()
  
  # Check GIs v DOLP ages. Assumes core age model is adequate
  mint<-round(min(dbiob$temp), -2)                                         # safe place to plot GIs
  maxt<-max(dbiob$temp)
  offs<-(maxt-mint)*0.2
  pg<-ggplot(GIs)+  
    geom_point(aes(ageb2K, y = mint), col = "blue", shape = 24, size = 4)+ # GIs in range                                    # 
    geom_text(aes(ageb2K, y = mint, label = number))+                      # GIs in range
    #
    geom_line(data = dbiob, aes(age,temp))  +                              # base recon
    #geom_point(data = cgdd, aes(age,X4), col = "#0072B2")  +               # base recon
    geom_vline(data = dbiob, aes(xintercept = ifelse(dbiob$OlDOLP == TRUE,dbiob$age,NA)), col = "grey50")+   # Not chosen DOLPs
    geom_vline(data = N2left, aes(xintercept = age), col = "red")+         # Surviving DOLPs
    annotate("text", x = 11000, y = (mint+maxt)/2,                         # How many found compared with GIs
             label = paste(nrow(N2left),"/",nrow(GIs)))+
    ylab(biovar)+
    xlim(c(10000,max(dbiob$age)))+
    #scale_x_continuous(breaks = seq(from = 10000, to = max(dbiob$age), by = 10000), 
    #                  limits =c(10000,max(dbiob$age)))+
    theme_bw()
  # Set parameters to find SD of bootstrapped recon; borrowed from 'Tel recon TWAPLS Jul12'
  set<-"12Jul23"          # defined in ...
  version1<-"1"
  if(file.exists(paste("C:/Users/markg/Dropbox/UoR/TWAPLS/Telford SD of reconstructions", 
                       biovar, sites[[p]], version1, set, ".csv"))) {           # no Telford bootstrapped for restricted
    # Get file
    this<-read.csv(paste("C:/Users/markg/Dropbox/UoR/TWAPLS/Telford SD of reconstructions", 
                         biovar, sites[[p]], version1, set, ".csv"))
    pg<-pg+geom_ribbon(data = this, aes(x = age, ymin = recm-2*y2, ymax = recm+2*y2), # Tel recon
                       col = "blue", fill = "skyblue", alpha = 0.3)
  }
  #
  sp<-egg::ggarrange(p1,g0, g1,g2,g3,g4,pg, nrow = 7,
                     heights = c(0.13,0.13,0.13,0.13,0.13,0.13,0.22) )
  grid::grid.draw(sp)
  list<-list(N2left, Donly, dbiob)
  names(list)<-c("N2left", "Donly", "dbiob")
  return(list)
} 
# end of filter_warm function

# ~~~~~~~~~~~~~~~ confirm_filtered description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Takes output of 'filter_warm' function, and removes DOLPs rejected by filtering..
#
# Writes out files containing 'dbiob', 'dslim', and 'filt', which lists the filter settings used.
# Of these, dbiob is the key output for use downstream, and contains much 
# information collected en route.

confirm_filtered<-function(N2left,
                           Donly,
                           dbiob,
                           vout) {
  print(corename)
  # Ones to bin
  trybin<- !(Donly$DOno %in% N2left$DOno)
  Dbin<-data.frame(Donly[trybin,"DOno"])
  
  # Blank out DOLPS which don't make the cut
  dbt<-dbiob
  for (k in 1:nrow(Dbin)){
    # Find binnable DOLP
    print(k)
    inx<-which(dbt$DOno == Dbin[k,1])     # index of the offending DOLP
    dbt[inx, "DOLP"]<-FALSE               
    dbt[inx, "stend"]<-NA                 # remove 'start'
    dbt[inx, "gbump"]<-NA                 # remove 'bumps'
    dbt[inx, "warea"]<-NA                 # remove gdd area-under-curve
    dbt[inx, "wrange"]<-NA                # remove gdd range
    dbt[inx, "mrange"]<-NA                # remove crtmi range
    
    if(nrow(Dbin) > 1){                   # v painful if only 1
      # Work backwards to eliminate its traces
      while (dbt[max(inx,1), "warm"] == 1 & inx > 0) {  # will stop at "end"   
        dbt[inx,"warm"]<-0
        dbt[inx,"DOno"]<-NA
        dbt[inx-1, c("stend","wrange","gbump","mbump","rbump","tbump")]<-NA  # remove 'bumps'
        inx <-inx-1
      }
    }
  }  
  #
  print(paste("Confirm: no of DOLPS left",sum(dbt$DOLP, na.rm = TRUE)))     # no of DOLPs after filtering
  # Accept results of final filters
  dbiob<-dbt                              # overwrite dbiob with filtered version
  dslim<-dbiob[dbiob$DOLP == TRUE, ]
  # Write the final choice of warm patch data
  # Concatenation of vtemp and vfilt -> ensures base params AND filters used are recorded
  write.csv(dbiob, paste0("testWarm patches detail for ",corename,vout,".csv"))
  write.csv(dslim, paste0("testWarm patches summary for ",corename,vout,".csv"))
  filt<-c(corename, dfilt, rfilt, afilt)
  write.csv(filt, paste0("testWarm patches filters for ",corename,biovar,".csv"))
}
# End




 

