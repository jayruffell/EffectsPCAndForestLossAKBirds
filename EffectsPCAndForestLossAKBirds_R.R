#=============================
#notes
#=============================

# ***NOTE in section 0 there is some code which is only relevant if I'm using raw counts as myresponse var (ie GLMMs) - eg survey conditions data. Don't read in 
#    ...otherwise! ***
# Copied and modified from ?ANALYSIS_CouncilY1to4? file in Auckland Council folder.

#=============================
# BOOKMARKS
#=============================

#0:   Read in data, and calculate PC per buffer
#0a:  Subset the data by tier, when I get the chance, so that whole rest of analyses doesnt need changing!
#0b:  No code to run here, but a note: need to decide prior to analyses whether I'm using the data summed within sites (for GLM's) or 
#     ...raw counts (for GLMMs). Depending on which one I use, the appropriate bird df needs to be read in, and the right set of code 
#     ...above needs to be hashed in/hashed out. It is also possible to use an offset variable to calculate summed data without 
#     ...removing sites <3 counts (see Zuur MM GEE chapters)
# 1. data preparation: transformations and PCA
# 1.a: create climate PC axes
# 1.b: Transform vars
# 2. graph correlations among variables:
# 3. Statistical modelling of community level vars - richness and total abundance          ***FINAL***
# 3.a poisson models, no overdispersion or autocorrelation accounted for    
# 3.b model validation ? overdispersion, autocorrelation, plots          *** Sp.nat ***     
# 3.c plotting fitted values          *** Sp.nat ***     
# 3.d repeating 3c, but with PatchPC cf PatchPC5yr          *** Sp.nat ***     
# 3.e model validation ? overdispersion, autocorrelation, plots          *** Ab.nat ***     
# 3.f plotting fitted values          *** Ab.nat ***     
# 3.g repeating 3c, but with PatchPC cf PatchPC5yr          *** Ab.nat ***  
# 3.h  stitching plots together          ***FINAL***                          [see also 3.l]
# 3.i re-running final models with "No control" as baseline, to get param estimates for tables,     ***FINAL***
#...and calculating sample size and NBT range for each PC category    
# 3.j re-plotting fitted values including 'vegClass' in models      ***Sp.nat***
# 3.k re-plotting fitted values including 'vegClass' in models      ***Ab.nat*** 
# 3.l  stitching vegClass plots together        
# 3.m re-plotting fitted values excluding re-intro Species      ***Sp.nat***
# 3.n re-plotting fitted values excluding re-intro Species      ***Ab.nat*** 
# 3.o  stitching ex-reintro plots together        
# 4. Statistical modelling of individual species                 ***FINAL***
# 4.a poisson models, no overdispersion or autocorrelation accounted for    
# 4.b model validation ? overdispersion, autocorrelation, plots          *** Sp.nat ***     
# 4.c plotting fitted values          *** Sp.nat ***     
# 4.d repeating 4c, but with PatchPC cf PatchPC5yr          *** Sp.nat ***     
# 4.e model validation ? overdispersion, autocorrelation, plots          *** Ab.nat ***     
# 4.f plotting fitted values          *** Ab.nat ***     
# 4.g repeating 4c, but with PatchPC cf PatchPC5yr          *** Ab.nat ***  
# 4.h  stitching plots together          ***FINAL***                          [see also 4.l]
# 4.i re-running final models with "No control" as baseline, to get param estimates for tables,     ***FINAL***
#...and calculating sample size and NBT range for each PC category    
# 4.j re-plotting fitted values including 'vegClass' in models      ***Sp.nat***
# 4.k re-plotting fitted values including 'vegClass' in models      ***Ab.nat*** 
# 4.l  stitching vegClass plots together        
# 4.m re-plotting fitted values excluding re-intro Species      ***Sp.nat***
# 4.n re-plotting fitted values excluding re-intro Species      ***Ab.nat*** 
# 4.o  stitching ex-reintro plots together        



#_________________________________________________________________________________________________________________________________________________________________________________________________________

#0. READ IN DATA, INCLUDING CALCULATING PC PER BUFFER      
#!*** NB THIS SECTION HAS DIFF CODE DEPENDING ON WHETHER USING "RAW" OR "SUMMED" BIRD DATA. NEED TO ENSURE UNUSED ONE IS "HASH OUT" ***
#!*** NB ALSO NEED TO DEFINE PEST CONTROL INTENSITY SCORES IN THIS SECTION, TO CALCULATE PC PER BUFFER
#_________________________________________________________________________________________________________________________________________________________________________________________________________                                                         


#=============================
#Read in Function to provide information on missing rows in a merge:
#=============================

#preliminary stuff:
mergeInfo <- function(x,y){
  
  x.name <-deparse(substitute(x))   #obtains data frame names for subsequent printouts
  y.name <-deparse(substitute(y))
  
  if(length(x)>2) {x <- x} else {x$tempVarX <- 1:nrow(x)}#This is just a dirty workaround for a prob I had below. The "rowSums" function won't work on
  if(length(y)>2) {y <- y} else {y$tempVarY <- 1:nrow(y)}#..DFs with<2 vars, so just this adds a fake var in these cases to make it work.
  
  #merge:
  merge1 <- merge(x, y, all.x=T)
  merge2 <- merge(x, y, all.y=T)
  
  #IDing rows within merge with ANY missing values (ie not just whole row missing):
  anyMissing1 <- merge1[!complete.cases(merge1),]  
  anyMissing2 <- merge2 [!complete.cases(merge2),]         
  
  #IDing rows COMPLETELY missing from merge (cf single N/A here and there)
  firstCol <- length(x)+1
  lastCol <- length(x)+length(y)-1     #this bit removes columns from other DF
  missing1 <- merge1[rowSums(is.na(merge1[,firstCol:lastCol])==1)==length(merge1[,firstCol:lastCol]),]#this bit extracts those columns with ALL NAs
  
  firstCol <- 2
  lastCol <- length(x)
  missing2 <- merge2[rowSums(is.na(merge2[,firstCol:lastCol])==1)==length(merge2[,firstCol:lastCol]),] 
  
  #printing output
  print(paste("All rows in",x.name,"COMPLETELY absent from", y.name,":"))
  print(missing1)
  print("-------------------------------------")
  print(paste("All rows in",y.name,"COMPLETELY absent from", x.name,":"))
  print(missing2)
  print("-------------------------------------")
  print(paste("Number of rows,",x.name,":",(nrow(x))))
  print(paste("Number of rows,",y.name,":",(nrow(y))))
  print("-------------------------------------")
  print(paste("Number of rows in",x.name,"COMPLETELY absent from", y.name,":",(nrow(missing1))))
  print(paste("Number of rows in",y.name,"COMPLETELY absent from", x.name,":",(nrow(missing2))))    
  print("-------------------------------------")
  print(paste("Number of rows in",x.name,"-",y.name,"merge with ANY missing values:",(nrow(anyMissing1))))
  print(paste("Number of rows in",y.name,"-",x.name,"merge with ANY missing values:",(nrow(anyMissing2))))    
  print("-------------------------------------")
  print("great function jay. you're coding like a boss.")
}



#=============================
# Read in DFs from "Preparing data for analysis" output, and remove missing rows/NAs if necessary (looking at you, veg):
#=============================

bird <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Birds.AllTiers.SUMMEDAcrossCounts.txt", header = TRUE)

#  bird <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Birds.AllTiers.SUMMEDAcrossCounts_Sub20m.txt", header = TRUE)
#  bird <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Birds.AllTiers.SUMMEDAcrossCounts_NoReintro.txt", header = TRUE)

nrow(bird); nrow(subset(bird, nCounts==3))    #lose ~30 sites if excluding sites with <3 counts - down to 249 (cf ~150 for CTC data)
bird <- subset(bird, nCounts==3)  #dont think offset will work with my data, cos not a linear increase with number of counts

#=============================
#!######!*** BIRD DATA IS ALL TIERS, + SUMMED ACROSS COUNTS. DECIDE ABOVE IF DOING <20M ONLY TOO, AND <3 COUNTS!!!***
#=============================

str(bird)                       #249 sites ?many sites missing because of not having three counts
bird <- subset(bird, select=-X)                                   #removing random "FID" variable
bird[!complete.cases(bird),]                       #no NAs
head(bird)
bird$Site
hist(bird$Sp.nat)    #maybe zero inflated if I use distance<20
graphics.off()  

#============================                                                            #===============================
### preparing site_count data for use with ***RAW*** bird count data [ie for glmms]:           *** DON'T READ IN OTHERWISE! ***
#============================                                                            #===============================

###--> Site var is actually Site_count for Bird.RAW data, e.g. "BAL_count1". Need to rename:
#    bird$Site_count <- bird$Site

### And need to create a "Site" variable that does not include count number:
#    truncRight <- function(x, n){substr(x, 1, nchar(x)-n)}  # user-defined function to remove n characters from end of STRING. (wont work for factors!)
#    bird$Site <- truncRight(as.character(bird$Site_count), 7)   #function works on characters, not factors
#    bird$Site <- as.factor(bird$Site)

###Move the variable "Site_count" to being the first, rather than the last variable:
#        Site_count <- subset(bird, select=Site_count) 
#        bird <- merge(Site_count, bird, by="Site_count")

#          str(bird)
#          bird[!complete.cases(bird),]                       #no NAs

###adding bird count conditions data:        
#    birdConds <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/BSum_CountConditions.txt", header = TRUE)
#    birdConds <- read.delim(file = "C:/Users/20901402/Documents/Auckland Council data and analyses/BSum_CountConditions.txt", header = TRUE)
#    birdConds = read.delim(file = "C:/Documents and Settings/Administrator/My Documents/Auckland Council data and analyses/BSum_CountConditions.txt", header = TRUE)
#    str(birdConds)
#         birdConds <- subset(birdConds, select=-c(X, Site))                                   #removing random "FID" variable

### i noticed in merge a couple of sites get lost. Find out why:        
#      tempDF <- merge(bird, birdConds, by="Site_count", all.x=T)
#     tempDF[!complete.cases(tempDF),]                                #missing all counts from CG40DD and CK41D, + 1 count from CO42. Are these sites important?
#        tempDF <- merge(bird, birdConds, by="Site_count", all.y=T)
#        tempDF[!complete.cases(tempDF),]

#    bird <- merge(bird, birdConds, by="Site_count")
#    nrow(bird)
#    length(unique(bird$Site))  
#    bird[1:10,]


#============================                                                            #===============================
### preparing site_count data for use with *SUMMED*** bird count data:                    *** DON'T READ IN OTHERWISE! ***
#============================                                                            #===============================

###--> this code adds 'ObserverID' and 'DaysSinceNov1' to Site-level data, since these vars are the same for all counts within a site.
#       ...BUT DON'T ADD THIS IN IF i'M USING ***RAW*** DATA, COS WILL HAVE ALREADY BEEN ADDED IN CODE ABOVE!!!

###adding bird count conditions data:        
birdConds <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/BSum_CountConditions.txt", header = TRUE)
str(birdConds)
###changing 'Site_count' to 'Site' for merge with other dfs:
truncRight <- function(x, n){substr(x, 1, nchar(x)-n)}  # user-defined function to remove n characters from end of STRING. (wont work for factors!)
birdConds$Site <- truncRight(as.character(birdConds$Site_count), 7)   #function works on characters, not factors
birdConds$Site <- as.factor(birdConds$Site)
birdConds <- subset(birdConds, select=-Site_count)
birdConds2 <- subset(birdConds, select=c(Site, Observer, DaysSinceNov1))      #these are the only 'site-level' vars
subset(birdConds2, Site=='TAG_23')
###since occasionally a single site had multiple observers, using aggregate() to pick most common observer per site: 
calcMode <- function(x){names(sort(-table(x)))[1]}  # need user-defined function to calc mode
birdConds3 <- aggregate(Observer ~ Site, data=birdConds2, FUN=calcMode)
birdConds3$Observer <- as.factor(birdConds3$Observer)                   #above function returns character, not factor
birdConds3
str(birdConds3)
#and repeating for daysSinceNov1 to the mix, then merging:
birdConds4 <- aggregate(DaysSinceNov1 ~ Site, data=birdConds2, FUN=calcMode)
birdConds4$DaysSinceNov1 <- as.numeric(birdConds4$DaysSinceNov1)                   #above function returns character
birdConds4[1:10,]
str(birdConds4)
#merging:
birdConds5 <- merge(birdConds3, birdConds4, by="Site")
#and adding a time since sunrise, rain, sun, noise, and wind value, averaged across counts:                                 
head(birdConds)
birdConds6a <- aggregate(Sun ~ Site, data=birdConds, FUN=mean)
birdConds6b <- aggregate(Rain ~ Site, data=birdConds, FUN=mean)
birdConds6c <- aggregate(Wind ~ Site, data=birdConds, FUN=mean)                              
birdConds6d <- aggregate(Noise ~ Site, data=birdConds, FUN=mean)
birdConds6e <- aggregate(MinSinceSunrise ~ Site, data=birdConds, FUN=mean)                    
birdConds6 <- merge(birdConds6a,birdConds6b, by="Site")
birdConds6 <- merge(birdConds6,birdConds6c, by="Site")
birdConds6 <- merge(birdConds6,birdConds6d, by="Site")
birdConds6 <- merge(birdConds6,birdConds6e, by="Site")
birdConds7 <- merge(birdConds5,birdConds6, by="Site")
head(birdConds7)

#=======================
#And back to sfuff to read in for all data, not just if I'm using raw [GLMM] Bird Count data
#=======================

veg <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/V.sum.txt", header = TRUE)
veg <- subset(veg, select=-X)
str(veg)
veg <- subset(veg, Site!="CK38C2")  #just found from NMDS exercise that this site was all macrocarpas - so removing
veg <- subset(veg, select=-AvDBH)  #this variable already accounted for by PCA
nrow(veg)  #282 sites

# CTC data (not pest control data ? see below for this):
pest <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Pest.txt", header = TRUE)
pest$nCTCcounts <- pest$nCounts; pest <- subset(pest, select=-nCounts) # renaming to avoid clash with variable name in Bird data
pest <- subset(pest, select=-c(X,Year))
pest[!complete.cases(pest),]                       #no NAs
pest$predRisk <- with(pest, ifelse(Rat>Poss, Rat, Poss))      #creating PredRisk var - whichever is greater out of rat and poss.
str(pest)                 #only 219 observations ? versus 286 for veg.
#renaming incorrectly named sites (easier to do here than in orig excel files & BSVP files)
levels(pest$Site)[match("TANGA1",levels(pest$Site))] <- "TANGA"
levels(pest$Site)[match("CR32BI7",levels(pest$Site))] <- "CR32B17"
levels(pest$Site)[match(" CP30AD4",levels(pest$Site))] <- "CP30AD4"

LandCov1000 = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/LandCover1000.txt", header = TRUE)
LandCov1000 <- subset(LandCov1000, select=-X)                                   #removing random "FID" variable
LandCov1000$NBT <- LandCov1000$NATBRO+LandCov1000$teatree        #now defining habitat as NBT.
LandCov1000 <- subset(LandCov1000, select=-c(NATBRO, teatree))   #these no longer needed
str(LandCov1000)
levels(LandCov1000$Site)
IslandMainland = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Island.mainland.txt", header = TRUE)
str(IslandMainland)

UrbanYN = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/BVSP data prepared in R/Urban.YN.txt", header = TRUE)
UrbanYN <- subset(UrbanYN, select=c(Site, Urban.YN))                                   #removing random "FID" variable
str(UrbanYN)


#=============================
# DFs From "Auckland council data and analyses" folder
#=============================

PAreaAndEdgeD = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/PatchIDPatchAreaEdgeDist_PhD+CouncilCOMBINED.txt", header = TRUE)
PAreaAndEdgeD$Site <- toupper(PAreaAndEdgeD$Site)
PAreaAndEdgeD <- subset(PAreaAndEdgeD, Site != "CE32BA")  #this site had no patch as far as I could tell! Middle of pine forest. Best to just remove - one 
str(PAreaAndEdgeD)                                        #...of many Tapora sites anyway.

lenz = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/LENZ2_PhD+CouncilCOMBINED.txt", header = TRUE)
str(lenz)

frag_AK = read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Fragstats_NATBRO_AKY1to4.txt", header = TRUE)       

frag_PhD = read.delim(file = "C:/Users/new user/Documents/PhD Study Sites_includes some council analyses/Fragstats_NATBRO_PhD.txt", header = TRUE)       
frag <- rbind(frag_PhD, frag_AK)
colnames(frag)[which(colnames(frag)=="TITLE")] <- "Site"   #! *** FRAG SITES ARE COUNCIL ONLY! NEED FOR PH.D. SITES AS WELL!***
str(frag)

coords <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/GPS coords_PhD+CouncilCOMBINED.txt", header = TRUE)
coords$XCOORD <- as.numeric(coords$XCOORD)
coords$YCOORD <- as.numeric(coords$YCOORD)   #were stored as text
str(coords)

tier <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Tiers_CouncilY1to4.txt", header = TRUE)
str(tier)
tier$Tier <- as.factor(tier$Tier)

#adding Ph.D. sites and assigning them a value of " tier one ":     
#! NEED TO MANUALLY ADD PHD SITES FOR NEXT YEAR!
SitePhD <- c("BAL", "LEN", "MAP", "RID", "BEN", "TANGA", "TANGB", "KAK2", "GC1", "DRE")
myTiers <- rep(1, times=length(SitePhD))
tierPhD <- as.data.frame(cbind(SitePhD, myTiers))
str(tierPhD)
colnames(tierPhD) <- c("Site", "Tier")
tier <- rbind(tier,tierPhD)

#now, 2dfs which were created in R but stored in AkCouncil folder:
vegClasses <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/vegClasses_PhD+CouncilCOMBINED.txt", header = TRUE)
vegClasses <- subset(vegClasses, select=-X)
str(vegClasses)

vegPCA <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/vegStructurePCAscores.txt", header = TRUE)
vegPCA <- subset(vegPCA, select=-X)
str(vegPCA)

#raw climate data (used to create PCA):
clim <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/ClimateSlopeAltitude_PhD+CouncilCOMBINED.txt", header = TRUE)

#climate PCA (differs from pest climate PCA in that all sites are included):
climPCA_bird <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/Analyses Y1-4 data+PhD sites/climatePCAscores_bird.txt", header = TRUE)
climPCA_bird <- subset(climPCA_bird, select=-X)
str(climPCA_bird)

#=============================
# read in PC data
#============================= 
pc <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/PestControlPerPatchAndBuffer_PhD+CouncilCOMBINED_1000m_AllPC_FINAL.txt", header = TRUE)
pc <- subset(pc, select=-NBT) #to avoid subsequent variable clash
str(pc)
pc[!complete.cases(pc),]                       #no NAs
pc$PatchPC <- factor(pc$PatchPC, levels=c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")) #makes pred free baseline

#and add 5yr PatchPC (for bird sensitivity analysis):
pc5yr <- read.delim(file = "C:/Users/new user/Documents/Auckland Council data and analyses/PestControlPerPatchAndBuffer_PhD+CouncilCOMBINED_1000m_5YrPC_FINAL.txt", header = TRUE)
pc5yr$PatchPC5yr <- pc5yr$PatchPC 
pc5yr <- subset(pc5yr, select=c(Site,PatchPC5yr)) 
pc[!complete.cases(pc5yr),]                       #no NAs
pc5yr$PatchPC5yr <- factor(pc5yr$PatchPC5yr, levels=c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")) #makes pred free baseline
str(pc5yr)
pc <- merge(pc, pc5yr, by="Site")    
str(pc)

#=================================
#! CHANGING TAPORA SITES **PATCHPC** TO INTERMEDIATE (WILL DO WITHIN GIS SOON, BUT FOR NOW THIS IS QUICKER)  
#=================================
pc[which(pc$Site=="CD32D"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CD32DA"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32AB"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32AD"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32BD"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32D"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE33A"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CF40DB16"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CD32AB"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CD32CA"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32BA"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE33AB"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CE32A"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CD32A"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"
pc[which(pc$Site=="CD33A"),which(colnames(pc)=="PatchPC")] <- "Intermediate RP"

##NB IntermedRP started 2007, with cyclic poss since 2004 - so manually ammending for _5yr PC (looking in ArcGIS surv year)
pc[which(pc$Site=="CD32D"),which(colnames(pc)=="PatchPC5yr")] <- "Cyclic possum"
pc[which(pc$Site=="CD32DA"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32AB"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32AD"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32BD"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32D"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"   #all surveyed in Y4 except 2y3 sites, 
pc[which(pc$Site=="CE33A"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"   #..1y2 and 1y1
pc[which(pc$Site=="CF40DB16"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CD32AB"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CD32CA"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32BA"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE33AB"),which(colnames(pc)=="PatchPC5yr")] <- "Intermediate RP"
pc[which(pc$Site=="CE32A"),which(colnames(pc)=="PatchPC5yr")] <- "Cyclic possum"
pc[which(pc$Site=="CD32A"),which(colnames(pc)=="PatchPC5yr")] <- "Cyclic possum"
pc[which(pc$Site=="CD33A"),which(colnames(pc)=="PatchPC5yr")] <- "Cyclic possum"

#=============================
# CALCULATE PEST CONTROL PER LANDSCAPE BASED ON INTENSITY SCORES:
#=============================

#==> actually, have deleted this since I'm only planning on using patch PC (highly correlated with landscape PC).
#... If I do want to calculate, see the Auckland CTC manuscript code

#=============================
# merge all data frames together:
#=============================

#GETTING INFO ON ANY MISSING SITES:

mergeInfo(coords,pc)    
mergeInfo(coords,LandCov1000)
mergeInfo(coords,IslandMainland) 
mergeInfo(coords,UrbanYN)
mergeInfo(coords,lenz)
mergeInfo(coords,frag)
mergeInfo(coords,tier)           #All rows match and no missing values for this paragraph. yeah boi!

mergeInfo(coords,PAreaAndEdgeD)   #single site is not represented in coords - ce32ba. This is a tapora site, and from memory it was all pine, so couldnt calc p size.

mergeInfo(coords,vegClasses)    #veg is only missing "CK38C2" - this is macrocarpa site I removed intentionally. coords is missing CD32AD, which i believe is a typo but i cant figure out how to fix
mergeInfo(coords,veg)           #identical to above. Yeah boi!
mergeInfo(coords,vegPCA)        #veg PCA is missing 5 sites: CK38C2 above, + 4 which had missing values for CanCov (had to remove for PCA)
mergeInfo(coords,climPCA_bird)    

mergeInfo(coords,pest)         #coords not missing any. as expected, pest has 65 sites missing.

mergeInfo(coords,birdConds7) #all birdConds sites match coords, but 7 coords sites are missing from birdConds. I already know 5 of these are sites 
#...where no birds were surveyed. Of other two, 1 is in Waitaks, 1 is urban AK

mergeInfo(coords,bird)       #33 bird sites missing 
mergeInfo(birdConds7,bird)   #...birdConds is missing 2 sites from bird

# MERGE ALL GIS DFS:
gisData <- merge(LandCov1000,IslandMainland, by="Site") 
gisData <- merge(gisData,UrbanYN, by="Site") 
gisData <- merge(gisData,PAreaAndEdgeD, by="Site")            
gisData <- merge(gisData,lenz, by="Site") 
gisData <- merge(gisData,frag, by="Site") 
gisData <- merge(gisData,coords, by="Site") 
gisData <- merge(gisData,tier, by="Site")   
gisData <- merge(gisData,pcPerBuffAll, by="Site")   
gisData <- merge(gisData,pcPerBuff5Yr, by="Site")   
gisData[!complete.cases(gisData),]; nrow(gisData)   #no NAs & 281 sites. Single site lost from PatchArea merge.
nrow(subset(gisData, Mainland.Island=="mainland"))  #225 mainland sites

veg2 <- merge(veg, vegClasses, by="Site")
veg2 <- merge(veg2, vegPCA, by="Site")

#MERGING BIRD DATA WITH OTHER DFs:
bird <- merge(bird, birdConds7, by="Site")
nrow(bird)                                       #275 sites (when all counts inc.)
bird <- merge(bird, pc, by="Site")
nrow(bird)                                       

#=============================
# including/excluding vegetation data:
#=============================
#  bird <- merge(bird, veg2, by="Site")          # hash out if not using veg data -- get another four sites this way

#...And back to other merges:
bird <- merge(bird,gisData, by="Site")
bird <- merge(bird,climPCA_bird, by="Site")
nrow(bird)                                       #269 sites ( a few sites were lost in veg PCA cos of missing data)
nrow(subset(bird, Mainland.Island=="mainland"))  # 215 sites (< 20 m counts)               
bird <- subset(bird, Site!="CK38C2")  #all macrocarpas - so removing [have already removed from veg data; only need to do again if not merging vegetation]
head(bird)

#=======================================
#removing/renaming a bunch of variables that have been superseded: 
#=======================================
#...? e.g. for the Shannon indices, I think species richness is just as good. And for vegetation variables,
#... I don't need all of those variables and some are redundant (e.g. "tree abundance" has been replaced by tree density in PCA; 
#... don't need tree diversity and vegetation diversity).

names(bird)        
# bird <- subset(bird, select=-c(Shan.tot, Shan.ex, Shan.nat, V.Sp.nat, V.Sp.ex, V.Sp.tot, V.Prop.nat, Tree.Ab.tot, Tree.Ab.nat, 
#         Tree.Ab.ex, Tree.Sp.nat, Tree.Sp.ex, Tree.Shan.tot,  Tree.Shan.nat, Tree.Shan.ex, Tree.Prop.Sp.nat,
#         vegStrPC2, vegStrPC3, climatePC2, grassland, urban))#, BEITAR, VITLUC, PRUFER, METROB))
###SUBSET WONT WORK WITH HASHED OUT VEG2 MERGE ABOVE
colnames(bird)[which(colnames(bird)=="grey.warbler")] <- "greywarbler"
colnames(bird)[which(colnames(bird)=="North.Island.kaka")] <- "kaka"
colnames(bird)[which(colnames(bird)=="North.Island.tomtit")] <- "tomtit"
colnames(bird)[which(colnames(bird)=="New.Zealand.fantail")] <- "fantail"
colnames(bird)[which(colnames(bird)=="grey.warbler")] <- "greywarbler"
colnames(bird)[which(colnames(bird)=="grey.warbler")] <- "greywarbler"    
names(bird)        

#=============================
# ADDING CTC DATA, BUT SAVING AS SEPARATE DF COS OF ALL THE LOST SITES:
#=============================
birdPest <- merge(bird,pest, by="Site")
nrow(birdPest)                                    ###==> 179 SITES
str(birdPest)



#=======================================
####### TRANSFORMATIONS COPIED FROM LOWER DOWN #############
#=======================================

#note variables that were not transformed either were beyond hope (i.e. mega zero inflated, as most land Variables were)or already normal
#...or are response variables

#bird df:
bird_t <- bird                                                
bird_t$EdgeD <- log(bird_t$EdgeD+0.01)
bird_t$maxDBH <- log(bird$maxDBH) #wont run if I don't merge veg2 data
bird_t$pine <- log(bird$pine+0.01)
bird_t$NBT <- sqrt(bird_t$NBT+0.01) 
bird_t$PatchArea_ha <- log(bird_t$PatchArea_ha+0.01) 
bird_t$Tree.Sp.tot <- sqrt(bird_t$Tree.Sp.tot+0.01) #as above
names(bird_t)

#birdPest df:
birdPest_t <- birdPest                                             
birdPest_t$EdgeD <- log(birdPest_t$EdgeD+0.01)
birdPest_t$maxDBH <- log(birdPest$maxDBH) #as above
birdPest_t$pine <- log(birdPest$pine+0.01) 
birdPest_t$NBT <- sqrt(birdPest_t$NBT+0.01) 
birdPest_t$PatchArea_ha <- log(birdPest_t$PatchArea_ha+0.01) 
birdPest_t$Prop.Ab.nat <- logit(birdPest_t$Prop.Ab.nat)
birdPest_t$Tree.Sp.tot <- sqrt(birdPest_t$Tree.Sp.tot+0.01) #as above
names(birdPest_t)

### Edge effects stuff From preliminary analyses ? may still be useful for the final stuff?:
#ED.patch <- ED[ED$Patch.ID=="LBI"|ED$Patch.ID=="GBI"|ED$Patch.ID=="Waitaks"|    # defining new DF where only multiple sites per patch (more powerful test of edge effects cos can control for patch characteristics). 
#            ED$Patch.ID=="Hunua"|ED$Patch.ID=="Hunua.North",]
#ED.patch$Patch.ID <- factor(ED.patch$Patch.ID)                                  # neeed to reset factor levels for coplot to work
#ED.patch$Patch.ID2 <- ifelse(ED.patch$PatchPC=="Intensive multispecies", paste(ED.patch$Patch.ID,"_intensePC"),paste(ED.patch$Patch.ID))    # when looking at edge effects within large patches need to treat KMA/AIP separately - 
#ED.patch$Patch.ID2 <- factor(ED.patch$Patch.ID2)                                                                                           # ... otherwise 'edge effects' may be PC effects (KMA is in middle of patch)
#ED$Patch.ID <- factor(ED$Patch.ID)

#=============================
# remove island sites:
#=============================

bird <- subset(bird, Mainland.Island=="mainland")  
bird_t <- subset(bird_t, Mainland.Island=="mainland")
birdPest <- subset(birdPest, Mainland.Island=="mainland")
nrow(bird)
nrow(bird_t)

#=============================
### prior to modelling, need to rescale continuous variables (so says lmer)
#=============================
bird_Rescaled <- bird
bird_Rescaled$NBT <- scale(bird_Rescaled$NBT)
bird_Rescaled$climatePC1 <- scale(bird_Rescaled$climatePC1)
bird_Rescaled$Year <- scale(bird_Rescaled$Year)
bird_Rescaled$XCOORD <- scale(bird_Rescaled$XCOORD)
bird_Rescaled$vegStrPC1 <- scale(bird_Rescaled$vegStrPC1) #wont run if i dont include veg2 merge above
bird_Rescaled$Tree.Sp.tot <- scale(bird_Rescaled$Tree.Sp.tot) #as above
bird_Rescaled$Tree.Prop.Ab.nat <- scale(bird_Rescaled$Tree.Prop.Ab.nat) #as above
bird_Rescaled$maxDBH <- scale(bird_Rescaled$maxDBH) #as above
bird_Rescaled$DaysSinceNov1 <- scale(bird_Rescaled$DaysSinceNov1)
bird_Rescaled$MinSinceSunrise <- scale(bird_Rescaled$MinSinceSunrise)
bird_Rescaled$Sun <- scale(bird_Rescaled$Sun)
bird_Rescaled$Rain <- scale(bird_Rescaled$Rain)
bird_Rescaled$Wind <- scale(bird_Rescaled$Wind)
bird_Rescaled$Noise <- scale(bird_Rescaled$Noise)
# subsequently adding transformed NBT, since untransformed version doesn't fit data well:
bird_Rescaled$sqrtNBT <- scale(sqrt(bird$NBT))
bird_Rescaled$logNBT <- scale(log(bird$NBT + 0.1))
bird_Rescaled$ID <- seq(from=1, to=nrow(bird_Rescaled))    


# reading in vegetation PCA variables pre-PCA - easier than calculating in situ from veg data if I dont want to use PCA:
vegVars <- read.csv(file = "C:/Users/new user/Documents/Auckland Council data and analyses/vegStructureVarsPrePCA.csv", header = TRUE)
vegVars <- subset(vegVars, select=-X)
str(vegVars)
####NB cancov untransformed; canopy height square root transformed; others log transformed.

library("latticeExtra")                                                           # for doing dot plots with transparent points
library("AICcmodavg")                                                   #AIC model averaging, also predictions w SE for mixed modelsthis
library("MuMIn")                                                   #AIC model averaging, also predictions w SE for mixed modelsthis
library("psych")                                                        #pairs panels
library("car")                                                          #variance inflation factor is
library("ncf")                                                          #spatial autocorrl correlograms)
library("ggplot2")
library("gridExtra")
library("lme4")
library("afex")                                                           # calculate P values for lmer objects 
library("multcomp")
library("zoo")                                                            #for merging by *nearest*, cf exact, value  
library("boot")
library("nlme")
#*** NB ABOVE SECTION HAS DIFF CODE DEPENDING ON WHETHER USING "RAW" OR "SUMMED" BIRD DATA, AND <20 M COUNTS OR ALL.
#.... NEED TO ENSURE UNUSED ONE IS "HASH OUT" ***

#######################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
########################################################################################################################################################################################################
##################################################################################()######################################################################################################################


#______________________________________________________________________________________________________________________________________________________

#1.a.  create PCA axis for climate and altitude variables
#______________________________________________________________________________________________________________________________________________________


### NB differs from PCA axis used for pests; this is only calculated scores for the pest subset of sites

#=============================
#run PCA
#=============================

# SEE http://www.youtube.com/watch?annotation_id=annotation_510466&src_vid=5zk93CpKYhg&v=oZ2nfIPdvjY for methods.

# WARNING: PCA ASSUMES (1) LINEAR RELNS (2) NORMALITY (3) CONTINUOUS VARS (4) MORE?    
# NB PCA DONE ON TRANSFORMED VARS - BUT EARLIER PC ON UNTRANSFORMED GAVE IDENTICAL RESULTS.

PCAvars  <- subset(clim, select=c(Slope, Elev, TempMin, TempAv))
PCA <- prcomp(PCAvars, center=T, scale=T)                    # read somewhere that 'prcomp' function is generally slightly superior to 'princomp' and should be used.
summary(PCA)                                                 # alsways need to centre the data, but scaling (ie making var=1) only needed when vars measured on very diff scales.


# HOW MANY COMPONENTES TO KEEP? 

# 1. Under "80% rule" (see Zuur Book1):  KEEP FIRST 2 (explain 90% variance)
# 2. Using 'Kaiser criterion', keep any with eigenvalue >1.0:
PCA$sdev^2                                                  # --> keep first 1, which explain 69% of variance

# 3. scree plot:
screeplot(PCA, main = "NB 'Kaiser criterion' suggests keeping 1st components", xlab="Components") 
# suggests keeping only first

# WHICH VARS ARE DRIVING TRENDS:

biplot(PCA, cex=0.5)   #Temp vars and elevn are strongly correlated, slope is almost orthogonal. Increasing axis values reflect high elev and low temps.
#...and to a much lesser extent, low slopes. (worth including is slope as its own variable?)
varimax(PCA$rotation)      # varimax rotation 'tweaks' orignal rotation- uses orig rotation variabless as input. Makes it easier to understand which vars are contributing to each component. 

PCA     #-->all variables loaded on to the first axis, but only slopes loads onto the second



#=============================
# extract PCA scores and write results to csv
#=============================

str(PCA)
climatePCA_bird <- PCA$x[,1:2]     #only need first PC
colnames(climatePCA_bird) <- c("climatePC1","climatePC2")
climatePCA_bird

#attaching site labels 
Site <- as.character(clim[,1])
climatePCA_bird <- cbind(Site,climatePCA_bird) 
climatePCA_bird <- as.data.frame(climatePCA_bird)
climatePCA_bird$climatePC1 <- as.numeric(as.character(climatePCA_bird$climatePC1)) #as.dF() converted scores to factors, so this is how they have to be converted back to numeric
climatePCA_bird$climatePC2 <- as.numeric(as.character(climatePCA_bird$climatePC2))
str(climatePCA_bird)

#Checking scores and sites match up with original data:
climatePCA_bird[1:2,]; PCA$x[1:2,]      #yip, scores match up with orignial scores.

subset(climatePCA_bird, climatePC1==max(climatePCA_bird$climatePC1))
subset(climatePCA_bird, climatePC1==min(climatePCA_bird$climatePC1))  #yip, the sites with the most extreme vals for PC1 also have extreme vals for temp/elevn
subset(clim, Site=="CM30C58")
subset(clim, Site=="CN42")


#and finally ...
mergeInfo(clim,climatePCA_bird)   #no missing rows
write.csv(climatePCA_bird, file="climatePCAscores_bird.csv")




#______________________________________________________________________________________________________________________________________________________

# 1.b. transform vars
#______________________________________________________________________________________________________________________________________________________

#================
#plot diff transformations:
#================

#UPDATE: redoing vars selected as arcsine transformation with a logit transform, since this is better:

library("car")                                                        #logit() function

op <- par(mfrow=c(3,3))
hist(ifelse(bird$Prop.Ab.nat>99, asin(sqrt((99.9)/100)), asin(sqrt((bird$Prop.Ab.nat+0.01)/100)))) # trans won't work for 100%
hist(logit(bird$Prop.Ab.nat))
hist(bird$Prop.Ab.nat)
hist(ifelse(bird$pcPerBuffAll>99, asin(sqrt((99.9)/100)), asin(sqrt((bird$pcPerBuffAll+0.01)/100)))) # trans won't work for 100%
hist(logit(bird$pcPerBuffAll))
hist(bird$pcPerBuffAll)
hist(ifelse(bird$pcPerBuff5Yr>99, asin(sqrt((99.9)/100)), asin(sqrt((bird$pcPerBuff5Yr+0.01)/100)))) # trans won't work for 100%
hist(logit(bird$pcPerBuff5Yr))
hist(bird$pcPerBuff5Yr)
par(op)                              #yip logit does a good job


#MAIN VARS:
names(bird)
plotsDF <- subset(bird, select=-c(Site, vegClass, kereru, tui, grey.warbler, New.Zealand.fantail, silvereye, North.Island.kaka, 
                                  North.Island.tomtit, Observer, Tier, BEITAR, VITLUC, PRUFER, METROB, Mainland.Island, Urban.YN, PatchID, 
                                  LENZ2, LENZ1, PLAND, NP, PD, TE, ED, PARA_MN, ENN_MN, CLUMPY, PatchPCAll, PatchPC5Yr))   
#--> only keeping vars I'm interested in transforming for subsequent data eyeballing


VAR.NAMES  = names(plotsDF)
for(i in 1:length(VAR.NAMES)){  
  X.VAR.i = plotsDF[,VAR.NAMES[i]]                                   
  YourFileName=paste("BoxplotsAndHists_",VAR.NAMES[i],".jpg",sep="")
  jpeg(file=YourFileName, res=300,width=4000,height=2000)                               
  par(mfrow=c(2,4))                              
  if(max(X.VAR.i)<100){                                 #ifelse is necessary because arcsine wont work - and will stop code - for values >100
    boxplot(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
    boxplot(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
    boxplot(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
    boxplot(asin(sqrt((X.VAR.i+0.01)/100)), main=paste(VAR.NAMES[i], " arc.sqrt+0.01"))
    hist(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
    hist(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
    hist(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
    hist(asin(sqrt((X.VAR.i+0.01)/100)), main=paste(VAR.NAMES[i], " arc.sqrt+0.01"))
  }else{ 
    boxplot(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
    boxplot(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
    boxplot(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
    boxplot(1:10, main="IGNORE: CANT PLOT ARCSINE, MAX>100")
    hist(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
    hist(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
    hist(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
    hist(1:10, main="IGNORE: CANT PLOT ARCSINE, MAX>100")             
  }      
  dev.off()
}
graphics.off()





#CTC VARS
plotsDF <- subset(birdPest, select=c(Rat, Poss, Mouse))

VAR.NAMES  = names(plotsDF)
for(i in 1:length(VAR.NAMES)){  
  X.VAR.i = plotsDF[,VAR.NAMES[i]]                                   
  YourFileName=paste("BoxplotsAndHists_",VAR.NAMES[i],".jpg",sep="")
  jpeg(file=YourFileName, res=300,width=4000,height=2000)                               
  par(mfrow=c(2,4))                              
  boxplot(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
  boxplot(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
  boxplot(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
  boxplot(asin(sqrt((X.VAR.i+0.01)/100)), main=paste(VAR.NAMES[i], " arc.sqrt+0.01"))         
  hist(X.VAR.i, main=paste(VAR.NAMES[i],"raw"))     
  hist(sqrt(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " sqrt+0.01"))
  hist(log(X.VAR.i+0.01), main=paste(VAR.NAMES[i], " log+0.01"))
  hist(asin(sqrt((X.VAR.i+0.01)/100)), main=paste(VAR.NAMES[i], " arc.sqrt+0.01"))         
  dev.off()
}
graphics.off()

#OTHER LATECOMER VARS:
op <- par(mfrow=c(2,2))
hist(bird$pine)
hist(log(bird$pine))   # easily best
hist(sqrt(bird$pine))
hist(logit(bird$pine))
par(op)

op <- par(mfrow=c(2,2))
hist(bird$maxDBH)
hist(log(bird$maxDBH))   # easily best
hist(sqrt(bird$maxDBH))
hist(logit(bird$maxDBH))
par(op)

#================
#transform and save to df:
#================

#note variables that were not transformed either were beyond hope (i.e. mega zero inflated, as most land Variables were)or already normal

bird_t <- bird                                                
bird_t$Ab.ex <- sqrt(bird_t$Ab.ex+0.01) 
bird_t$Ab.tot <- sqrt(bird_t$Ab.tot+0.01) 
bird_t$Ab.nat <- sqrt(bird_t$Ab.nat+0.01) 
bird_t$EdgeD <- log(bird_t$EdgeD+0.01) 
bird_t$NBT <- sqrt(bird_t$NBT+0.01) 
bird_t$PatchArea_ha <- log(bird_t$PatchArea_ha+0.01) 
bird_t$Prop.Ab.nat <- logit(bird_t$Prop.Ab.nat)
bird_t$Sp.ex <- sqrt(bird_t$Sp.ex+0.01)
bird_t$Sp.nat <- log(bird_t$Sp.nat+0.01)
bird_t$Tree.Sp.tot <- sqrt(bird_t$Tree.Sp.tot+0.01)
bird_t$pcPerBuffAll <- logit(bird_t$pcPerBuffAll)
bird_t$pcPerBuff5Yr <- logit(bird_t$pcPerBuff5Yr)
names(bird_t)

birdPest_t <- birdPest                                             
birdPest_t$Ab.ex <- sqrt(birdPest_t$Ab.ex+0.01) 
birdPest_t$Ab.tot <- sqrt(birdPest_t$Ab.tot+0.01) 
birdPest_t$Ab.nat <- sqrt(birdPest_t$Ab.nat+0.01) 
birdPest_t$EdgeD <- log(birdPest_t$EdgeD+0.01) 
birdPest_t$NBT <- sqrt(birdPest_t$NBT+0.01) 
birdPest_t$PatchArea_ha <- log(birdPest_t$PatchArea_ha+0.01) 
birdPest_t$Prop.Ab.nat <- logit(birdPest_t$Prop.Ab.nat)
birdPest_t$Sp.ex <- sqrt(birdPest_t$Sp.ex+0.01)
birdPest_t$Sp.nat <- log(birdPest_t$Sp.nat+0.01)
birdPest_t$Tree.Sp.tot <- sqrt(birdPest_t$Tree.Sp.tot+0.01)
birdPest_t$pcPerBuffAll <- logit(birdPest_t$pcPerBuffAll)
birdPest_t$pcPerBuff5Yr <- logit(birdPest_t$pcPerBuff5Yr)
birdPest_t$predRisk <- logit(birdPest_t$predRisk)
birdPest_t$Rat <- logit(birdPest_t$Rat)
birdPest_t$Poss <- logit(birdPest_t$Poss)
names(birdPest_t)




#_________________________________________________________________________________________________________________________________________________________________________________________________________

#2. graph correlations among variables:
#_________________________________________________________________________________________________________________________________________________________________________________________________________                                                         

#=============================
#first, look at correlations within each "variable group":
#=============================

names(bird_t)

pairs.panels(subset(bird_t, select=c(Sp.nat,NBT,pine,ocean,EdgeD, PatchArea_ha, XCOORD)))
pairs.panels(subset(bird_t, select=c(Ab.tot, Ab.nat, Ab.ex, Sp.tot, Sp.nat, Sp.ex, Prop.Sp.nat, Prop.Ab.nat)))
pairs.panels(subset(bird_t, select=c(Sp.nat,Tree.Sp.tot, Tree.Prop.Ab.nat, maxDBH, vegStrPC1)))
pairs.panels(subset(bird_t, select=c(NBT, Tree.Sp.tot, Tree.Prop.Ab.nat, maxDBH, vegStrPC1, climatePC1)))
pairs.panels(subset(bird_t, select=c(NBT, DaysSinceNov1, Sun,Rain,Wind,Noise,MinSinceSunrise))) # a few correlations

#and some graphs to look for confounding with categorical variables:
op <- par(mfrow=c(2,2))
plot(NBT ~PatchPC5yr, data=bird_t)
plot(NBT ~Urban.YN, data=bird_t)      # NBT verses other categorical variables
plot(NBT ~vegClass, data=bird_t)
plot(NBT ~Observer, data=bird_t)
par(op) 

op <- par(mfrow=c(3,2))
plot(Tree.Sp.tot ~PatchPC5yr, data=bird_t)
plot(Tree.Prop.Ab.nat ~PatchPC5yr, data=bird_t)
plot(Tree.Sp.tot ~PatchPC5yr, data=bird_t)
plot(vegStrPC1 ~PatchPC5yr, data=bird_t)      # pc versus other variables
plot(maxDBH ~PatchPC5yr, data=bird_t)
plot(climatePC1 ~PatchPC5yr, data=bird_t)
par(op) 

op <- par(mfrow=c(3,2))
plot(MinSinceSunrise ~PatchPC5yr, data=bird_t)
plot(Noise ~PatchPC5yr, data=bird_t)
plot(Wind ~PatchPC5yr, data=bird_t)
plot(Rain ~PatchPC5yr, data=bird_t)      # pc versus other variables, take 2
plot(Sun ~PatchPC5yr, data=bird_t)
plot(DaysSinceNov1 ~PatchPC5yr, data=bird_t)
par(op) 

###FOR BOX & WHISKER PLOTS, CHECKING WHETHER THESE RELATIONSHIPS ARE SIGNIFICANT WITH STATISTICAL MODELS:

# NBT verses other categorical variables
mixed(NBT ~PatchPC5yr + (1|PatchID), data=bird_t)     #'mixed()', from afex package, calculates P values for anova test of lmer models.
mixed(NBT ~Urban.YN + (1|PatchID), data=bird_t)
summary(glht(lmer(NBT ~Urban.YN + (1|PatchID), data=bird_t), mcp(Urban.YN="Tukey")))  #glht pvals near identical... excellent   
mixed(NBT ~vegClass + (1|PatchID), data=bird_t)
summary(glht(lmer(NBT ~vegClass + (1|PatchID), data=bird_t), mcp(vegClass="Tukey")))  
#==> urban and scrub vs natbro vary with NBT (PC does too, but unimportant since I'm including both in models
# NBT versus continuous variables:
mixed(XCOORD~ NBT + (1|PatchID), data=bird_t) 
mixed(pine ~NBT + (1|PatchID), data=bird_t) 
mixed(Sun ~NBT + (1|PatchID), data=bird_t) 
mixed(Rain ~NBT + (1|PatchID), data=bird_t)       #p=0.06
mixed(Wind ~NBT + (1|PatchID), data=bird_t) 
mixed(Noise ~NBT + (1|PatchID), data=bird_t)     #0.001
mixed(MinSinceSunrise ~NBT + (1|PatchID), data=bird_t)    #highly sig
mixed(DaysSinceNov1 ~NBT + (1|PatchID), data=bird_t)      #0.04
mixed(Tree.Sp.tot ~NBT + (1|PatchID), data=bird_t)        #0.0003
mixed(Tree.Prop.Ab.nat ~NBT + (1|PatchID), data=bird_t) 
mixed(maxDBH ~NBT + (1|PatchID), data=bird_t)           # 0.07
mixed(vegStrPC1 ~NBT + (1|PatchID), data=bird_t)           
mixed(climatePC1 ~NBT + (1|PatchID), data=bird_t)           # high sig
mixed(vegClass ~NBT + (1|PatchID), data=bird_t)           # 0.003 

# PC vs veg
mixed(Tree.Sp.tot ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(Tree.Prop.Ab.nat ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(vegStrPC1 ~PatchPC5yr + (1|PatchID), data=bird_t)    # 0.06       #should probably have PC as a function of vegetation, not the other way around?
mixed(maxDBH ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(vegClass~PatchPC5yr + (1|PatchID), data=bird_t)     #0.08
mixed(climatePC1 ~PatchPC5yr + (1|PatchID), data=bird_t) #climate highly sig effect, but shouldnt affect birds?
mixed(XCOORD~PatchPC5yr + (1|PatchID), data=bird_t)   
summary(glht(lmer(climatePC1 ~PatchPC5yr + (1|PatchID), data=bird_t), mcp(PatchPC5yr="Tukey")))  #PP diff from LRP 
#==> climate vary with PatchPC, but shouldnt affect birds?
# PC vs count conditions
mixed(MinSinceSunrise ~PatchPC5yr + (1|PatchID), data=bird_t)
summary(glht(lmer(MinSinceSunrise ~PatchPC5yr + (1|PatchID), data=bird_t), mcp(PatchPC5yr="Tukey")))  #PP diff from LRP 
mixed(Noise ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(Wind ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(Rain ~PatchPC5yr + (1|PatchID), data=bird_t)     
summary(glht(lmer(Rain ~PatchPC5yr + (1|PatchID), data=bird_t), mcp(PatchPC5yr="Tukey")))  #PP highly sig rainier than LRP 
mixed(Sun ~PatchPC5yr + (1|PatchID), data=bird_t)
mixed(DaysSinceNov1 ~PatchPC5yr + (1|PatchID), data=bird_t)
summary(glht(lmer(DaysSinceNov1 ~PatchPC5yr + (1|PatchID), data=bird_t), mcp(PatchPC5yr="Tukey")))  #few differences here, but none with NC 
#==> a few significant differences, but none between NC and other categories                                          

#=============================                                            
#now checking for confounding among potentially important variables:      ###OBSOLETE###
#=============================

### NBT versus everything else:

names(bird_t)

Y.VAR.NAMES  = c("NBT")   #needs to be yvar, cos cant have factors as yvar...
X.VAR.NAMES  = c("Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.nat","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","Year","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit","Observer","DaysSinceNov1","Tier","Tree.Sp.tot","Tree.Prop.Ab.nat","BEITAR","VITLUC",
                 "PRUFER","METROB","grassland","pine","urban","ocean","Mainland.Island","Urban.YN","PatchArea_ha","EdgeD","PatchID","LENZ2",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3")            

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- bird_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland"|PatchID=="LBI") #removing island sites except LBI
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    if(is.numeric(X.VAR.i)=="TRUE") {                                       # if-else required because the loess function wont work on factors, eg Mainland.island, and stops loop.
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i],
           main="blue triangles=LBI sites, red circles=my sites. \nAll other island sites removed", 
           col=ifelse(tempDF$PatchID=="LBI","blue",ifelse(tempDF$Observer=="Jay","red","black")),
           pch=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)),
           cex=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)))
      M.Loess=loess(Y.VAR.j~X.VAR.i)                                        # ...this says 'only do Loess if the XVAR is numeric'
      Fit=fitted(M.Loess)
      Ord1=order(X.VAR.i)
      lines(X.VAR.i[Ord1], Fit[Ord1], col="red")
      dev.off()
    } else {
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      dev.off()
    }
  }
}
dev.off()


### Observer versus everything else:  [NB confounding is only important if observer correlates with important variables - e.g. observer ID probably correlates closely with shoe size, but who cares!]

X.VAR.NAMES  = c("Observer")   
Y.VAR.NAMES  = c("Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.nat","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","Year","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit","Observer","DaysSinceNov1","Tier","Tree.Sp.tot","Tree.Prop.Ab.nat","BEITAR","VITLUC",
                 "PRUFER","METROB","grassland","NBT","pine","urban","ocean","Mainland.Island","Urban.YN","PatchArea_ha","EdgeD","PatchID","LENZ2",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3")            

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- bird_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland") #removing island sites
  
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = bird_t[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    plot(Y.VAR.j~X.VAR.i, las=3, ylab=Y.VAR.NAMES[j], main="mainland only")
    graphics.off()
  }
}
dev.off()


### birds versus everything else:

Y.VAR.NAMES  = c("Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.nat","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit")
X.VAR.NAMES  = c("DaysSinceNov1","Tree.Sp.tot","Tree.Prop.Ab.nat","grassland","NBT","pine","urban","ocean","Urban.YN","PatchArea_ha","EdgeD",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3")            

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- bird_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland"|PatchID=="LBI") #removing island sites except LBI
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    if(is.numeric(X.VAR.i)=="TRUE") {                                       # if-else required because the loess function wont work on factors, eg Mainland.island, and stops loop.
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i],
           main="blue triangles=LBI sites, red circles=my sites. \nNo other island sites", 
           col=ifelse(tempDF$PatchID=="LBI","blue",ifelse(tempDF$Observer=="Jay","red","black")),
           pch=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)),
           cex=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)))
      M.Loess=loess(Y.VAR.j~X.VAR.i)                                        # ...this says 'only do Loess if the XVAR is numeric'
      Fit=fitted(M.Loess)
      Ord1=order(X.VAR.i)
      lines(X.VAR.i[Ord1], Fit[Ord1], col="red")
      dev.off()
    } else {
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      dev.off()
    }
  }
}
dev.off()

### repeating for birds vs pest vars, based on birdPest DF:
X.VAR.NAMES  = c("Mouse", "Rat", "Poss","predRisk")
Y.VAR.NAMES  = c("Sp.nat","Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit")            

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- birdPest_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland"|PatchID=="LBI") #removing island sites except LBI
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    if(is.numeric(X.VAR.i)=="TRUE") {                                       # if-else required because the loess function wont work on factors, eg Mainland.island, and stops loop.
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i],
           main="blue triangles=LBI sites, red circles=my sites. \nNo other island sites", 
           col=ifelse(tempDF$PatchID=="LBI","blue",ifelse(tempDF$Observer=="Jay","red","black")),
           pch=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)),
           cex=ifelse(tempDF$PatchID=="LBI",2,ifelse(tempDF$Observer=="Jay",2,1)))
      M.Loess=loess(Y.VAR.j~X.VAR.i)                                        # ...this says 'only do Loess if the XVAR is numeric'
      Fit=fitted(M.Loess)
      Ord1=order(X.VAR.i)
      lines(X.VAR.i[Ord1], Fit[Ord1], col="red")
      dev.off()
    } else {
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      dev.off()
    }
  }
}
dev.off()


### LENZ versus veg & LS vars:

Y.VAR.NAMES  = c("vegClass","vegStrPC1","vegStrPC2","vegStrPC3","Tree.Sp.tot","Tree.Prop.Ab.nat","NBT","pine","urban","ocean",
                 "Urban.YN","PatchArea_ha","EdgeD")
X.VAR.NAMES  = c("LENZ1")           

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- bird_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland") #removing island sites
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], main="mainland only")
    graphics.off()
  }
}
dev.off()


### Pest vars vs everything else:
Y.VAR.NAMES  = c("Mouse", "Rat", "Poss")
X.VAR.NAMES  = c("DaysSinceNov1","Tree.Sp.tot","Tree.Prop.Ab.nat","grassland","NBT","pine","urban","ocean","Urban.YN","PatchArea_ha","EdgeD","PatchID",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3",
                 "pcPerBuffAll","PatchPCAll","Mouse","Rat","Poss")            

for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  tempDF <- birdPest_t
  tempDF <- subset(tempDF, Mainland.Island=="mainland") #removing island sites
  Y.VAR.j = tempDF[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    if(is.numeric(X.VAR.i)=="TRUE") {                                       # if-else required because the loess function wont work on factors, eg Mainland.island, and stops loop.
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i], main="mainland only")
      M.Loess=loess(Y.VAR.j~X.VAR.i)                                        # ...this says 'only do Loess if the XVAR is numeric'
      Fit=fitted(M.Loess)
      Ord1=order(X.VAR.i)
      lines(X.VAR.i[Ord1], Fit[Ord1], col="red")
      dev.off()
    } else {
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      dev.off()
    }
  }
}
dev.off()


### Repeating for pest control variables, which I added later. MAINLAND SITES ONLY :

#pcPerBuff:

Y.VAR.NAMES  = c("pcPerBuffAll","pcPerBuff5Yr")
X.VAR.NAMES  = c("Sp.nat","Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","Year","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit","Observer","DaysSinceNov1","Tier","Tree.Sp.tot","Tree.Prop.Ab.nat","BEITAR","VITLUC",
                 "PRUFER","METROB","grassland","NATBRO","pine","teatree","urban","ocean","Mainland.Island","Urban.YN","PatchArea_ha","EdgeD","PatchID","LENZ2",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3")            

tempDF_t <- subset(bird_t, Mainland.Island=="mainland")     #removing island sites
for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  Y.VAR.j = tempDF_t[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF_t[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    if(is.numeric(X.VAR.i)=="TRUE") {                                       # if-else required because the loess function wont work on factors, eg Mainland.island, and stops loop.
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      M.Loess=loess(Y.VAR.j~X.VAR.i)                                        # ...this says 'only do Loess if the XVAR is numeric'
      Fit=fitted(M.Loess)
      Ord1=order(X.VAR.i)
      lines(X.VAR.i[Ord1], Fit[Ord1], col="red")
      dev.off()
    } else {
      plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i])
      dev.off()
    }
  }
}
dev.off()

#PatchPC:

X.VAR.NAMES  = c("PatchPCAll","PatchPC5Yr")
Y.VAR.NAMES  = c("Sp.nat","Ab.tot","Ab.nat","Ab.ex","Sp.tot","Sp.ex","Prop.Sp.nat","Prop.Ab.nat","Year","kereru","tui","grey.warbler","New.Zealand.fantail",
                 "silvereye","North.Island.kaka","North.Island.tomtit","Observer","DaysSinceNov1","Tier","Tree.Sp.tot","Tree.Prop.Ab.nat","BEITAR","VITLUC",
                 "PRUFER","METROB","grassland","NATBRO","pine","teatree","urban","ocean","Mainland.Island","Urban.YN","PatchArea_ha","EdgeD","PatchID","LENZ2",
                 "LENZ1","PLAND","NP","PD","TE","ED","PARA_MN","ENN_MN","CLUMPY","XCOORD","YCOORD","vegClass","vegStrPC1","vegStrPC2","vegStrPC3")            

tempDF_t <- subset(bird_t, Mainland.Island=="mainland")     #removing island sites
for(j in 1:length(Y.VAR.NAMES)){                                                # 'outer loop': loops through each bird variable
  Y.VAR.j = tempDF_t[,Y.VAR.NAMES[j]]
  
  for(i in 1:length(X.VAR.NAMES)){                                            # 'inner loop': for each bird var, loop through each Xvar and graph it
    X.VAR.i = tempDF_t[,X.VAR.NAMES[i]]
    YourFileName=paste(Y.VAR.NAMES[j], "vs.",X.VAR.NAMES[i],".jpg")
    jpeg(file=YourFileName)
    plot(Y.VAR.j~X.VAR.i, ylab=Y.VAR.NAMES[j], xlab=X.VAR.NAMES[i], main="mainland sites only")
    dev.off()
  }
}
dev.off()





#______________________________________________________________________________________________________________________________________________________

#3. preliminary statistical modelling            ***FINAL***
#______________________________________________________________________________________________________________________________________________________

#=============================
#initially , stuffign round with variance inflation factors to see which vars I should remove   ***requires library "car"***
#=============================

#=============================
### NB ZUUR (2010) SAYS: "one approach is to calculate the VIFs, then sequentially drop the term with the highest VIF until 
#..all VIFs are below particular threshold. 
#..Montgomery and Pack used a [VIF] value of 10,but a more stringent approach would be to use values as low as 3, 
#..as we did here. [When ecological signals are weak], even a VIF of 2 may cause non-significance" 
#=============================

###read in funciton to calc VIFs for GLMMs, from https://github.com/aufrank/R-hacks/blob/master/mer-utils.R:
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }   
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

#make no control baseline level of PC
bird_Rescaled$PatchPC5yr <- relevel(bird_Rescaled$PatchPC5yr, ref="No control")

#=============================
### prior to modelling, need to rescale continuous variables (so says lmer)         [already done in bookmark zero]
#=============================
bird_Rescaled <- bird
bird_Rescaled$NBT <- scale(bird_Rescaled$NBT)
bird_Rescaled$climatePC1 <- scale(bird_Rescaled$climatePC1)
bird_Rescaled$Year <- scale(bird_Rescaled$Year)
bird_Rescaled$XCOORD <- scale(bird_Rescaled$XCOORD)
bird_Rescaled$vegStrPC1 <- scale(bird_Rescaled$vegStrPC1). #won't read if I haven't read in the veg data. same for other veg vars.
bird_Rescaled$Tree.Sp.tot <- scale(bird_Rescaled$Tree.Sp.tot)
bird_Rescaled$Tree.Prop.Ab.nat <- scale(bird_Rescaled$Tree.Prop.Ab.nat)
bird_Rescaled$maxDBH <- scale(bird_Rescaled$maxDBH)
bird_Rescaled$DaysSinceNov1 <- scale(bird_Rescaled$DaysSinceNov1)
bird_Rescaled$MinSinceSunrise <- scale(bird_Rescaled$MinSinceSunrise)
bird_Rescaled$Sun <- scale(bird_Rescaled$Sun)
bird_Rescaled$Rain <- scale(bird_Rescaled$Rain)
bird_Rescaled$Wind <- scale(bird_Rescaled$Wind)
bird_Rescaled$Noise <- scale(bird_Rescaled$Noise)
# subsequently adding transformed NBT, since untransformed version doesn't fit data well:
bird_Rescaled$sqrtNBT <- scale(sqrt(bird$NBT))
bird_Rescaled$logNBT <- scale(log(bird$NBT+0.1))    

#______________________________________________________________________________________________________________________________________________________

# 3.a. poisson models, no overdispersion or autocorrelation accounted for       ***FINAL***
#______________________________________________________________________________________________________________________________________________________

#=============================
# STUFFING AROUND WITH INCLUDING CONFOUNDING VARIABLES (IDENTIFIED IN BOOKMARK 2), CHECKING VARIANCE INFLATION, TRANSFORMED VS UNTFORMED NBT,
# ...AND NBT*PC INTERACTION VS MAIN EFFECTS ONLY:
# (note originally did not transformed NBT, but model predictions were way off)
#=============================

## Note that'mFull' models won't run unless I've merged in the 'veg2' dataframe at the top of the code.

#=====================
#1 Sp.nat
#=====================

# a) using raw NBT
mFull_raw <- glmer(Sp.nat ~ NBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)       #may need to remove veg class because of high variance inflation

mMain_raw <- glmer(Sp.nat ~ NBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)

mInt_raw <- glmer(Sp.nat ~ NBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID), data=bird_Rescaled, family=poisson)   

# b) using sqrt NBT - fits best in terms of normalising NBT distribution
mFull_sqrt <- glmer(Sp.nat ~ sqrtNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)       #may need to remove veg class because of high variance inflation

mMain_sqrt <- glmer(Sp.nat ~ sqrtNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)

mInt_sqrt <- glmer(Sp.nat ~ sqrtNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   

# c) using log NBT - as per SP-Area curve (ie linear when log-log transformed)
mFull_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull_log)       #may need to remove veg class because of high variance inflation

mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)

mInt_log <- glmer(Sp.nat ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID), data=bird_Rescaled, family=poisson)   

AIC(mFull_raw, mMain_raw, mInt_raw, mFull_sqrt, mMain_sqrt, mInt_sqrt, mFull_log, mMain_log, mInt_log)


###==> very similar AIC for all models that are main effects only; log NBT is best though
###==> logNBT main is better than logNBT interaction - dAIC= 3.6

summary(mMain_log)  
summary(glht(mMain_log, mcp(PatchPC5yr="Tukey")))  #No sig diffs!!!!  

#==============
# dAIC is <4, so model averaging:
#==============
CM <- list()                   
CM[[1]] <-  mMain_log
CM[[2]] <-  mInt_log

model.sel(CM)                                           #===> actually, when only two models dAIC is 5.1; weight =0.93 
r.squaredGLMM(mMain_log)          #low R2...12% for fixed effects.


#=====================
#2 Ab.nat                  
#=====================

#Adding "ID" term added to account for overdisp. detected in bookmark #3.e
bird_Rescaled$ID <- 1:nrow(bird_Rescaled)
head(bird_Rescaled)

# a) using raw NBT
mFull_raw <- glmer(Ab.nat ~ NBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)       #may need to remove veg class because of high variance inflation

mMain_raw <- glmer(Ab.nat ~ NBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)

mInt_raw <- glmer(Ab.nat ~ NBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   

# b) using sqrt NBT - fits best in terms of normalising NBT distribution
mFull_sqrt <- glmer(Ab.nat ~ sqrtNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)       #may need to remove veg class because of high variance inflation

mMain_sqrt <- glmer(Ab.nat ~ sqrtNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull)

mInt_sqrt <- glmer(Ab.nat ~ sqrtNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   

# c) using log NBT - as per SP-Area curve (ie linear when log-log transformed)
mFull_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mFull_log)       #may need to remove veg class because of high variance inflation

mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
vif.mer(mMain_log)

mInt_log <- glmer(Ab.nat ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   

AIC(mFull_raw, mMain_raw, mInt_raw, mFull_sqrt, mMain_sqrt, mInt_sqrt, mFull_log, mMain_log, mInt_log)

###==> very similar AIC for all models that are main effects only; log NBT is best though

summary(mMain_log)  
summary(glht(mMain_log, mcp(PatchPC5yr="Tukey")))  #sig or near sig diffs between pred free and all except HRP, and betw HRP and LRP.

#==============
# dAIC is <4, so model averaging:
#==============
CM <- list()                   
CM[[1]] <-  mMain_log
CM[[2]] <-  mInt_log

model.sel(CM)                                           #===> actually, when only two models dAIC is 9.2; weight =0.99 
r.squaredGLMM(mMain_log)          #low R2...19% for fixed effects.


#______________________________________________________________________________________________________________________________________________________

# 3.b. model validation ? overdispersion, autocorrelation, plots        *** Sp.nat ***     
#______________________________________________________________________________________________________________________________________________________

mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   

#=============================
###TESTING GLMM OVERDISPERSION
#=============================

# first, try quasi poisson glm:
mQ <- glm(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise, data=bird_Rescaled, family=quasipoisson)
summary(mQ)  #no overdispersion.

#now with glmm overdisp function:
overdisp_fun <- function(model) {  
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))    #code emailed FROM RAPH (EMAILED ~27/2/13)
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mMain_log)                          #==> p=1.0   and quasipoisson glm says param = 0.3 (ie way less than 1; see above)


###DEALING WITH OVERDISPERSION, ALSO FROM RAPH'S E-MAIL:     

#"The suggested way of dealing with overdispersion (see Bolker blogs) is to add an observation-level random vector to the model (in addition to the other random effects). 
#..Practically, this means creating an ?obs? variable, and then adding (1|obs) to the model.
#People often cite Bolker et al 2009 for this (not sure if that?s correct though)...and Jiang 2007 (Jiang, J., 2007. Linear and Generalized Linear Mixed Models and Their 
#..Applications. Springer, New York, USA)....it is supposed to emulate a poisson-lognormal distribution of errors, as per the Elston paper (attached)"
#NB this approach makes little sense in the LMM case, where it would be confounded with the residual error)

#or same thing from http://comments.gmane.org/gmane.comp.lang.r.lme4.devel/7157:
#Dear Ricardo, You can add an observation-level random effect. That will capture the overdispersion. i.e.: 
#..yourData$ID <- seq_len(nrow(yourData)); m1 <- glmer(richness~leafs*galls+(1|Plant) + (1|ID),family=poisson, data = yourData). 
#...MCMCglmm uses the same trick to model overdispersion [another post on this website says use MCMCglmm,but this seems to be equivalent].          

#=============================
#testing autocorrelation:
#=============================
bird_Rescaled$LAT <- bird$XCOORD; bird_Rescaled$LONG <- bird$YCOORD  #couldn't use scaled versions of x & y, so had to re-import form unscaled df. 
CorGram <- spline.correlog(x=bird_Rescaled[,"LAT"], y=bird_Rescaled[,"LONG"], z=resid(mMain_log, type="pearson"))
plot.spline.correlog(CorGram)         
graphics.off()                    #no autocorrelation!

#=============================
# model validation plots:
#=============================
par(mfrow=c(3,3))                                                             
ResD <- resid(mMain_log, type="deviance")
plot(ResD~predict(mMain_log, type = "response"), main = "Deviance Resids vs fitted values")                                                                      
plot(ResD~bird_Rescaled$NBT, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$PatchPC5yr, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Urban.YN, main = "Deviance Resids vs Xvar")                     #all fine!
plot(ResD~bird_Rescaled$DaysSinceNov1, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$MinSinceSunrise, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Rain, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Noise, main = "Deviance Resids vs Xvar")


#______________________________________________________________________________________________________________________________________________________

# 3.c. plotting fitted values    *** Sp.nat ***     
#______________________________________________________________________________________________________________________________________________________

#=============================
# first, plotting the raw data:
#=============================

#untransformed NBT:
op <- par(mfrow=c(3,2))
plot(Sp.nat ~NBT, data=subset(bird, PatchPC5yr=="No control"), main="NC", ylim=c(0,11), xlim=c(0,100))
plot(Sp.nat ~NBT, data=subset(bird, PatchPC5yr=="Cyclic possum"), main="CC", ylim=c(0,11), xlim=c(0,100))
plot(Sp.nat ~NBT, data=subset(bird, PatchPC5yr=="Intermediate RP"), main="LRP", ylim=c(0,11), xlim=c(0,100))
plot(Sp.nat ~NBT, data=subset(bird, PatchPC5yr=="Intensive control"), main="HRP", ylim=c(0,11), xlim=c(0,100))
plot(Sp.nat ~NBT, data=subset(bird, PatchPC5yr=="Predator free"), main="PF", ylim=c(0,11), xlim=c(0,100))
par(op)
graphics.off()

#transformed sqrt(NBT+0.01):
op <- par(mfrow=c(3,2))
plot(Sp.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="No control"), main="NC", ylim=c(0,11), xlim=c(0,sqrt(100+0.01)))
plot(Sp.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Cyclic possum"), main="CC", ylim=c(0,11), xlim=c(0,sqrt(100+0.01)))
plot(Sp.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Intermediate RP"), main="LRP", ylim=c(0,11), xlim=c(0,sqrt(100+0.01)))
plot(Sp.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Intensive control"), main="HRP", ylim=c(0,11), xlim=c(0,sqrt(100+0.01)))
plot(Sp.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Predator free"), main="PF", ylim=c(0,11), xlim=c(0,sqrt(100+0.01)))
par(op)
graphics.off()


#=============================
# run model
#=============================
mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

SpNatVsNBT_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird richness") + ylim(0,10) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=10, label ="A.", hjust=0, size=5, fontface="bold")  #left justified
SpNatVsNBT_incLegend
SpNatVsNBT <- SpNatVsNBT_incLegend + theme(legend.position="none")  #removes legend


#______________________________________________________________________________________________________________________________________________________

# 3.d. repeating 3c, but with PatchPC cf PatchPC5yr         *** Sp.nat ***     
#______________________________________________________________________________________________________________________________________________________

#=============================
# run model
#=============================
mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

SpNatVsNBT_AllPC_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird richness") + ylim(0,10) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=10, label ="A.", hjust=0, size=5, fontface="bold")  #left justified
SpNatVsNBT_AllPC_incLegend
SpNatVsNBT_AllPC <- SpNatVsNBT_AllPC_incLegend + theme(legend.position="none")  #removes legend

png("Sp.nat~NBT+PCall.png",res=800,width=5500,height=4000)                  
SpNatVsNBT_AllPC_incLegend
graphics.off()

#xxx


#______________________________________________________________________________________________________________________________________________________

# 3.e. model validation ? overdispersion, autocorrelation, plots        *** Ab.nat ***     
#______________________________________________________________________________________________________________________________________________________

mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson) 

#updating to include ID random effect, since overdisp. function below showed overdispersion:
mMain_log2 <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson) 

#=============================
###TESTING GLMM OVERDISPERSION
#=============================

# first, try quasi poisson glm:
mQ <- glm(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise, data=bird_Rescaled, family=quasipoisson)
summary(mQ)  #overdispersion param=5!!

#now with glmm overdisp function:
overdisp_fun <- function(model) {  
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))    #code emailed FROM RAPH (EMAILED ~27/2/13)
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mMain_log)                          #==> p<0.0001 i.e. overdispersed!!!
overdisp_fun(mMain_log2)                          #==> FIXED!!!


###DEALING WITH OVERDISPERSION, ALSO FROM RAPH'S E-MAIL:     

#"The suggested way of dealing with overdispersion (see Bolker blogs) is to add an observation-level random vector to the model (in addition to the other random effects). 
#..Practically, this means creating an ?obs? variable, and then adding (1|obs) to the model.
#People often cite Bolker et al 2009 for this (not sure if that?s correct though)...and Jiang 2007 (Jiang, J., 2007. Linear and Generalized Linear Mixed Models and Their 
#..Applications. Springer, New York, USA)....it is supposed to emulate a poisson-lognormal distribution of errors, as per the Elston paper (attached)"
#NB this approach makes little sense in the LMM case, where it would be confounded with the residual error)

#or same thing from http://comments.gmane.org/gmane.comp.lang.r.lme4.devel/7157:
#Dear Ricardo, You can add an observation-level random effect. That will capture the overdispersion. i.e.: 
#..yourData$ID <- seq_len(nrow(yourData)); m1 <- glmer(richness~leafs*galls+(1|Plant) + (1|ID),family=poisson, data = yourData). 
#...MCMCglmm uses the same trick to model overdispersion [another post on this website says use MCMCglmm,but this seems to be equivalent].          

#=============================
#testing autocorrelation:
#=============================
bird_Rescaled$LAT <- bird$XCOORD; bird_Rescaled$LONG <- bird$YCOORD  #couldn't use scaled versions of x & y, so had to re-import form unscaled df. 
CorGram <- spline.correlog(x=bird_Rescaled[,"LAT"], y=bird_Rescaled[,"LONG"], z=resid(mMain_log2, type="pearson"))
plot.spline.correlog(CorGram)         
graphics.off()                    #no autocorrelation!

#=============================
# model validation plots:
#=============================
par(mfrow=c(3,3))                                                             
ResD <- resid(mMain_log2, type="deviance")
plot(ResD~predict(mMain_log2, type = "response"), main = "Deviance Resids vs fitted values")                                                                      
plot(ResD~bird_Rescaled$NBT, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$PatchPC5yr, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Urban.YN, main = "Deviance Resids vs Xvar")                     #all fine!
plot(ResD~bird_Rescaled$DaysSinceNov1, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$MinSinceSunrise, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Rain, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Noise, main = "Deviance Resids vs Xvar")

#______________________________________________________________________________________________________________________________________________________

# 3.f. plotting fitted values         *** Ab.nat ***     
#______________________________________________________________________________________________________________________________________________________

#=============================
# first, plotting the raw data:
#=============================

#untransformed NBT:
op <- par(mfrow=c(3,2))
plot(Ab.nat ~NBT, data=subset(bird, PatchPC5yr=="No control"), main="NC", ylim=c(0,90), xlim=c(0,100))
plot(Ab.nat ~NBT, data=subset(bird, PatchPC5yr=="Cyclic possum"), main="CC", ylim=c(0,90), xlim=c(0,100))
plot(Ab.nat ~NBT, data=subset(bird, PatchPC5yr=="Intermediate RP"), main="LRP", ylim=c(0,90), xlim=c(0,100))
plot(Ab.nat ~NBT, data=subset(bird, PatchPC5yr=="Intensive control"), main="HRP", ylim=c(0,90), xlim=c(0,100))
plot(Ab.nat ~NBT, data=subset(bird, PatchPC5yr=="Predator free"), main="PF", ylim=c(0,90), xlim=c(0,100))
par(op)
graphics.off()

#transformed sqrt(NBT+0.01):
op <- par(mfrow=c(3,2))
plot(Ab.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="No control"), main="NC", ylim=c(0,90), xlim=c(0,sqrt(100+0.01)))
plot(Ab.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Cyclic possum"), main="CC", ylim=c(0,90), xlim=c(0,sqrt(100+0.01)))
plot(Ab.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Intermediate RP"), main="LRP", ylim=c(0,90), xlim=c(0,sqrt(100+0.01)))
plot(Ab.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Intensive control"), main="HRP", ylim=c(0,90), xlim=c(0,sqrt(100+0.01)))
plot(Ab.nat ~sqrt(NBT+0.01), data=subset(bird, PatchPC5yr=="Predator free"), main="PF", ylim=c(0,90), xlim=c(0,sqrt(100+0.01)))
par(op)
graphics.off()


#=============================
# run model
#=============================
mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("E", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

AbNatVsNBT_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird abundance") + ylim(0,50) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=50, label ="B.", hjust=0, size=5, fontface="bold")  #left justified
AbNatVsNBT_incLegend
AbNatVsNBT <- AbNatVsNBT_incLegend + theme(legend.position="none")  #removes legend




#______________________________________________________________________________________________________________________________________________________

# 3.g. repeating 3c, but with PatchPC cf PatchPC5yr         *** Ab.nat ***     
#______________________________________________________________________________________________________________________________________________________

#=============================
# run model
#=============================
mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID) + (1|ID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.01))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))),                #...so shouldnt matter which i choose
                        ID = factor(rep("1",times=nrow(Xvar))))

# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

AbNatVsNBT_AllPC_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird abundance") + ylim(0,50) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=50, label ="B.", hjust=0, size=5, fontface="bold")  #left justified
AbNatVsNBT_AllPC_incLegend
AbNatVsNBT_AllPC <- AbNatVsNBT_AllPC_incLegend + theme(legend.position="none")  #removes legend

png("Ab.nat~NBT+PCall.png",res=800,width=5500,height=4000)                  
AbNatVsNBT_AllPC_incLegend
graphics.off()


#______________________________________________________________________________________________________________________________________________________

# 3.h.  stitching plots together          ***FINAL*** 
#______________________________________________________________________________________________________________________________________________________

# first, single plot to get figure legend:
png("Ab.nat~NBT+PC5yr_incLeg.png",res=800,width=5500,height=4000)                  
AbNatVsNBT_incLegend
graphics.off()

png("Sp.nat+Ab.nat~NBT+PC5yr.png",res=800,width=7000,height=4000)                  
grid.arrange(SpNatVsNBT, AbNatVsNBT, nrow=1, main=textGrob("",gp=gpar(fontsize=8))) # can't use mfrow() with ggplot
graphics.off()

png("Sp.nat+Ab.nat~NBT+PCALL.png",res=800,width=7000,height=4000)                  
grid.arrange(SpNatVsNBT_AllPC, AbNatVsNBT_AllPC, nrow=1, main=textGrob("",gp=gpar(fontsize=8))) # can't use mfrow() with ggplot
graphics.off()

#______________________________________________________________________________________________________________________________________________________

# 3.i.  re-running final models with "No control" as baseline, to get param estimates for tables,    ***FINAL***
#...and calculating sample size and NBT range for each PC category    
#______________________________________________________________________________________________________________________________________________________

#=============================
# rerunning models             ###NB for some reason package 'afex' prevents labelling of factors in output - so don't load!
#=============================

tempDF <- bird_Rescaled
tempDF$PatchPC5yr <- relevel(tempDF$PatchPC5yr, ref="No control")
#tempDF$PatchPC5yr <- factor(tempDF$PatchPC5yr, levels=c("Predator free","Intensive control","Intermediate RP","Cyclic possum",
#                            "No control"), labels=c("Predator free","Intensive control","Intermediate RP","Cyclic possum",
#                            "No control"))
mMain_log2 <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID), data=tempDF, family=poisson)   
summary(mMain_log2)

mMain_log2 <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID) + (1|ID), data=tempDF, family=poisson)   
summary(mMain_log2)


#=============================
# calculate sample size and range
#=============================

NBTnc <- subset(bird, select='NBT', PatchPC5yr=="No control")
nrow(NBTnc); min(NBTnc);max(NBTnc)
NBTlrp <- subset(bird, select='NBT', PatchPC5yr=="Intermediate RP")
nrow(NBTlrp); min(NBTlrp);max(NBTlrp)
NBTpp <- subset(bird, select='NBT', PatchPC5yr=="Cyclic possum")
nrow(NBTpp); min(NBTpp);max(NBTpp)
NBThrp <- subset(bird, select='NBT', PatchPC5yr=="Intensive control")
nrow(NBThrp); min(NBThrp);max(NBThrp)
NBTpf <- subset(bird, select='NBT', PatchPC5yr=="Predator free")
nrow(NBTpf); min(NBTpf);max(NBTpf)


#______________________________________________________________________________________________________________________________________________________

# 3.j. re-plotting fitted values including 'vegClass' in models      ***Sp.nat***
#______________________________________________________________________________________________________________________________________________________

#=============================
# run model                                    # NB have to un-hash veg data at read in (line ~474) for this bit
#=============================
mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + vegClass
                   + (1|Observer) + (1|PatchID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

SpNatVsNBT_vegClass_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird richness") + ylim(0,10) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=10, label ="A.", hjust=0, size=5, fontface="bold")  #left justified
SpNatVsNBT_vegClass_incLegend
SpNatVsNBT_vegClass <- SpNatVsNBT_vegClass_incLegend + theme(legend.position="none")  #removes legend

#______________________________________________________________________________________________________________________________________________________

# 3.k. re-plotting fitted values including 'vegClass' in models      ***Ab.nat*** 
#______________________________________________________________________________________________________________________________________________________

#=============================
# run model                                    # NB have to un-hash veg data at read in (line ~474) for this bit
#=============================
mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + vegClass
                   + (1|Observer) + (1|PatchID), data=bird_Rescaled, family=poisson)   
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        vegClass =factor(rep("stdNATBRO",times=nrow(Xvar)), levels =c("kakik","pohut","scrub","stdNATBRO")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

AbNatVsNBT_vegClass_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird abundance") + ylim(0,50) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=50, label ="B.", hjust=0, size=5, fontface="bold")  #left justified
AbNatVsNBT_vegClass_incLegend
AbNatVsNBT_vegClass <- AbNatVsNBT_vegClass_incLegend + theme(legend.position="none")  #removes legend


#______________________________________________________________________________________________________________________________________________________

# 3.l. Stitching vegClass plots together
#______________________________________________________________________________________________________________________________________________________

png("Sp.nat+Ab.nat~NBT+vegClass.png",res=800,width=7000,height=4000)                  
grid.arrange(SpNatVsNBT_vegClass, AbNatVsNBT_vegClass, nrow=1, main=textGrob("",gp=gpar(fontsize=8))) # can't use mfrow() with ggplot
graphics.off()


#______________________________________________________________________________________________________________________________________________________

# 3.m. re-plotting fitted values excluding reintroduced species      ***Sp.nat*** 
#______________________________________________________________________________________________________________________________________________________

### NB have to read in all data after reading in bird as ***Birds.AllTiers.SUMMEDAcrossCounts_NoReintro*** 
### ...instead of Birds.AllTiers.SUMMEDAcrossCounts

#=============================
# run model                                    
#=============================
mMain_log <- glmer(Sp.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise
                   + (1|Observer) + (1|PatchID), data=bird_Rescaled, family=poisson) 
summary(mMain_log)                                                    #NBT significant, PC no (but NC isnt baseline)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

SpNatVsNBT_NoReint_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird richness") + ylim(0,10) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=10, label ="A.", hjust=0, size=5, fontface="bold")  #left justified
SpNatVsNBT_NoReint_incLegend
SpNatVsNBT_NoReint <- SpNatVsNBT_NoReint_incLegend + theme(legend.position="none")  #removes legend

#______________________________________________________________________________________________________________________________________________________

# 3.n. re-plotting fitted values excluding reintroduced species      ***Ab.nat*** 
#______________________________________________________________________________________________________________________________________________________

### NB have to read in all data after reading in bird as ***Birds.AllTiers.SUMMEDAcrossCounts_NoReintro*** 
### ...instead of Birds.AllTiers.SUMMEDAcrossCounts

#=============================
# run model                                    # NB have to un-hash veg data at read in (line ~474) for this bit
#=============================
mMain_log <- glmer(Ab.nat ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise
                   + (1|Observer) + (1|PatchID), data=bird_Rescaled, family=poisson)
summary(mMain_log)

#=============================
# 1. FITTED VALUES FOR CONTROL=NC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("No control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df1 <- data.frame(Xvars,Means)#,loCI,upCI)
df1$PC <- as.factor(rep("NC", times=nrow(df1))) # name here determines legend labels

#=============================
# 2. FITTED VALUES FOR CONTROL=CC:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Cyclic possum",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df2 <- data.frame(Xvars,Means)#,loCI,upCI)
df2$PC <- as.factor(rep("PP", times=nrow(df2))) # name here determines legend labels

#=============================
# 3. FITTED VALUES FOR CONTROL=LRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intermediate RP",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df3 <- data.frame(Xvars,Means)#,loCI,upCI)
df3$PC <- as.factor(rep("LRP", times=nrow(df3))) # name here determines legend labels

#=============================
# 4. FITTED VALUES FOR CONTROL=HRP:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Intensive control",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df4 <- data.frame(Xvars,Means)#,loCI,upCI)
df4$PC <- as.factor(rep("HRP", times=nrow(df4))) # name here determines legend labels

#=============================
# 5. FITTED VALUES FOR CONTROL=PF:
#=============================

# specify X variable I'm predicting for:
Xvar <- data.frame(logNBT=seq(from=min(bird_Rescaled$logNBT), to=max(bird_Rescaled$logNBT), by=0.001))

# Specify other variables, whose levels I want to hold constant:        
OtherVars <- data.frame(PatchPC5yr = factor(rep("Predator free",times=nrow(Xvar)), levels =c("Predator free", "Intensive control", "Intermediate RP", "Cyclic possum", "No control")),
                        Urban.YN =factor(rep("N",times=nrow(Xvar)), levels =c("N","Y")),
                        DaysSinceNov1 = rep(mean(bird_Rescaled$DaysSinceNov1),times=nrow(Xvar)),   
                        MinSinceSunrise =rep(mean(bird_Rescaled$MinSinceSunrise),times=nrow(Xvar)),
                        Noise =rep(mean(bird_Rescaled$Noise),times=nrow(Xvar)),
                        Rain =rep(mean(bird_Rescaled$Rain),times=nrow(Xvar)),                    
                        Observer = factor(rep("Abigail Forbes",times=nrow(Xvar))),     #just using first value. random effects are estimated at zero, 
                        PatchID = factor(rep("Atiu",times=nrow(Xvar))))                #...so shouldnt matter which i choose
# predict data:
PlotData <- cbind(Xvar, OtherVars)
g <- predict(mMain_log, newdata=PlotData, type='response')#, se.fit=TRUE)
Xvars <- exp(PlotData$logNBT*sd(log(bird$NBT+0.1))+ mean(log(bird$NBT+0.1)))-0.1     #reversing transform process (which was take log, then stdize) 
Means <- g                                                                           #...to get back to orig scale
#loCI <- (g$fit-1.96*g$se)      ### no s.e's provided unless running through MuMIn
#upCI <- (g$fit+1.96*g$se)
df5 <- data.frame(Xvars,Means)#,loCI,upCI)
df5$PC <- as.factor(rep("PF", times=nrow(df5))) # name here determines legend labels

#=============================
### AND PLOTTING:
#=============================

# first, cutting predictions by range measured:
df1 <- subset(df1, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="No control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="No control")))
df2 <- subset(df2, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Cyclic possum")))
df3 <- subset(df3, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intermediate RP")))
df4 <- subset(df4, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Intensive control")))
df5 <- subset(df5, Xvars>=min(subset(bird$NBT, bird$PatchPC5yr=="Predator free")) & Xvars<=max(subset(bird$NBT, bird$PatchPC5yr=="Predator free")))

dfAll <- rbind(df5,df4,df3,df2,df1) #order here determines order of legend
head(dfAll)

AbNatVsNBT_NoReint_incLegend <- ggplot(data=dfAll, aes(x = Xvars, y = Means, lty=PC, col=PC)) + 
  geom_line() + 
  scale_linetype_manual("Pest control", values =c(1,2,3,2,1)) +   #manually assigning linetypes, + legend title          
  scale_color_manual("Pest control", values =c('dark grey','dark grey','black','black','black')) +   #manually assigning colours          
  xlab("Native forest cover (%)") + ylab("Bird abundance") + ylim(0,50) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # gives b&w graph instead of nasty default background  
  theme(axis.title=element_text(size=10)) +
  theme(legend.key = element_blank()) + scale_fill_discrete(name="myname") + 
  annotate("text", x=1, y=50, label ="B.", hjust=0, size=5, fontface="bold")  #left justified
AbNatVsNBT_NoReint_incLegend
AbNatVsNBT_NoReint <- AbNatVsNBT_NoReint_incLegend + theme(legend.position="none")  #removes legend


#______________________________________________________________________________________________________________________________________________________

# 3.o. Stitching ex-reintro plots together
#______________________________________________________________________________________________________________________________________________________

png("Sp.nat+Ab.nat~NBT_noReintro.png",res=800,width=7000,height=4000)                  
grid.arrange(SpNatVsNBT_NoReint, AbNatVsNBT_NoReint, nrow=1, main=textGrob("",gp=gpar(fontsize=8))) # can't use mfrow() with ggplot
graphics.off()

#basic poisson model, no overdispersion or zero inflation accounted for:

#______________________________________________________________________________________________________________________________________________________

#4. Statistical modelling of individual species                 ***FINAL***
#______________________________________________________________________________________________________________________________________________________

#make no control baseline level of PC
bird_Rescaled$PatchPC5yr <- relevel(bird_Rescaled$PatchPC5yr, ref="No control")

#calculating number of sites with non-zero counts (to get a rough idea of whether I need Poiss, NegBin, ZIP, or ZINB)
nrow(subset(bird_Rescaled, kereru>0)) #85 rows
nrow(subset(bird_Rescaled, tui>0)) #169 rows
nrow(subset(bird_Rescaled, fantail>0)) #158 rows
nrow(subset(bird_Rescaled, greywarbler>0)) #188 rows    #==> kaka is only one that is dominated by zeros! even tomtits occur at >25% of sites
nrow(subset(bird_Rescaled, silvereye>0)) #184 rows
nrow(subset(bird_Rescaled, tomtit>0)) #58 rows
nrow(subset(bird_Rescaled, kaka>0)) #15 rows

#getting a rough idea of whether there is zero inflation, by comparing observed data to hypothesised dist:

histVsPois2 <- function(sp){ 
  poissDat <- rpois(n=length(sp), lambda=mean(sp))
  plotDat <- data.frame(abund=c(sp, poissDat), 
                        datType=c(rep('observed', times=length(sp)), 
                                  rep('predicted', times=length(poissDat))))
  #qplot(abund, data=plotDat, fill=datType, title=deparse(substitute(sp)))
  ggplot(data=plotDat, aes(x=abund)) +
    geom_histogram(data=subset(plotDat, datType=='predicted'), fill='green', alpha=0.2) +
    geom_histogram(data=subset(plotDat, datType=='observed'), fill='blue', alpha=0.2) + 
    ggtitle(deparse(substitute(sp)))
  
} #function just simulates poiss data based on measured mean and then overlays hist. 
histVsPois2(bird$kereru)
histVsPois2(bird$tui)
histVsPois2(bird$fantail)       ###looks like overdisp and zero inflation for everything except GW and kaka.
histVsPois2(bird$greywarbler)
histVsPois2(bird$silvereye)    ###BUT take care - 'overall' dist isn't important, but rather dist for each level of xvars.
histVsPois2(bird$tomtit)       #...ie predictors may explain why there are so many zeros.
histVsPois2(bird$kaka)

#________________________________________________________________________________________________________________________________

# 4.a. kereru
#______________________________________________________________________________________________________________________________________________________

#=====================
#Comparing fit of different transformations of NBT, but not accounting for overdispersion or zero inflation yet: 
#=====================

## Note that'mFull' models won't run unless I've merged in the 'veg2' dataframe at the top of the code.

# a) using raw NBT
mFull_raw <- glmer(kereru ~ NBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   

mMain_raw <- glmer(kereru ~ NBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson,
                   glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser

mInt_raw <- glmer(kereru ~ NBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID), data=bird_Rescaled, family=poisson) 

# b) using sqrt NBT - fits best in terms of normalising NBT distribution
mFull_sqrt <- glmer(kereru ~ sqrtNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID), data=bird_Rescaled, family=poisson)

mMain_sqrt <- glmer(kereru ~ sqrtNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                    + (1|PatchID), data=bird_Rescaled, family=poisson,
                    glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser

mInt_sqrt <- glmer(kereru ~ sqrtNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson,
                   glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser   

# c) using log NBT - as per SP-Area curve (ie linear when log-log transformed)
mFull_log <- glmer(kereru ~ logNBT + PatchPC5yr + Urban.YN + climatePC1 + vegClass + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson)   

mMain_log <- glmer(kereru ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                   + (1|PatchID), data=bird_Rescaled, family=poisson,
                   glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser   

mInt_log <- glmer(kereru ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID), data=bird_Rescaled, family=poisson,
                  glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser   

AIC(mMain_raw, mInt_raw, mMain_sqrt, mInt_sqrt, mMain_log, mInt_log)#, mFull_raw, mFull_sqrt, mFull_log)


###==> log models are best; but interaction model is better than main effects only (dAIC~4)

summary(mInt_log)  #interesting... sig effects of all PC except PP (vs baseline and no control), plus sig interactions in which forest loss effect weakens w PC.
summary(glht(mInt_log, mcp(PatchPC5yr="Tukey")))  

#==============
# dAIC is <4, so model averaging:  ***BUT STILL NOT SURE IF MODELS ARE VALID***
#==============
CM <- list()                   
CM[[1]] <-  mMain_log
CM[[2]] <-  mInt_log

model.sel(CM)                                           #===> Need to average models; dAIC=2.51, weight =0.78 (in favour of int model)
r.squaredGLMM(mInt_log)          #highish R2...49% for fixed effects.


#________________________________________________________________________________________________________________________________

# 4.a.i model validation on best model
#______________________________________________________________________________________________________________________________________________________

mInt_log <- glmer(kereru ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                  + (1|PatchID), data=bird_Rescaled, family=poisson,
                  glmerControl(optimizer="bobyqa"))# model fails to converge with default optimiser

###TESTING GLMM OVERDISPERSION

# first, try quasi poisson glm:
mQ <- glm(kereru ~ logNBT + PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise, data=bird_Rescaled, family=quasipoisson)
summary(mQ)  # overdispersion param is 2.8 - think too high??

#now with glmm overdisp function:
overdisp_fun <- function(model) {  
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))    #code emailed FROM RAPH (EMAILED ~27/2/13)
  (rdf <- nrow(model@frame)-model.df)
  rp <- residuals(model)
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(mInt_log)                          #==> p=0.55; ie no sig overdispersion. But what about zero inflation? Also, histograms of observed count freqs
#...vs theoretical (see top of bookmark 4) suggest overdisp and zero inflation. 

#=============================
#testing autocorrelation:
#=============================
bird_Rescaled$LAT <- bird$XCOORD; bird_Rescaled$LONG <- bird$YCOORD  #couldn't use scaled versions of x & y, so had to re-import form unscaled df. 
CorGram <- spline.correlog(x=bird_Rescaled[,"LAT"], y=bird_Rescaled[,"LONG"], z=resid(mInt_log, type="pearson"))
plot.spline.correlog(CorGram)         
graphics.off()                    #no autocorrelation!

#=============================
# model validation plots:
#=============================
par(mfrow=c(3,3))                                                             
ResD <- resid(mInt_log, type="deviance")
plot(ResD~predict(mInt_log, type = "response"), main = "Deviance Resids vs fitted values")                                                                      
plot(ResD~bird_Rescaled$NBT, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$PatchPC5yr, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Urban.YN, main = "Deviance Resids vs Xvar")                     #hmm some strong patterning for several plots, and an outlier
plot(ResD~bird_Rescaled$DaysSinceNov1, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$MinSinceSunrise, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Rain, main = "Deviance Resids vs Xvar")
plot(ResD~bird_Rescaled$Noise, main = "Deviance Resids vs Xvar")

#with boxplots hard to tell if differing spread is real or an artefact of small sample size of some levels. Using jitterplots:
p1 <- qplot(y=ResD, x=bird_Rescaled$Urban.YN, geom="jitter", alpha=I(0.2), position = position_jitter(w = 0.1), main="Deviance resids vs Xvar")
p2 <- qplot(y=ResD, x=bird_Rescaled$PatchPC5yr, geom="jitter", alpha=I(0.2), position = position_jitter(w = 0.1),, main="Deviance resids vs Xvar")
ggplot2.multiplot(p1, p2, cols=2) 
#==> does appear to be differing spread for some factor levels.
#investigating outlier:
tempDF <- bird
tempDF$ResD <- ResD
subset(tempDF, ResD==max(tempDF$ResD)) #KMA site with 12 kereru
plot(bird$kereru)                    #but this isn't super unusual...several other sites with 8; one with 15!
#guess its just unusual in the context of the model.
#==> don't think outlier really is an outlier, could well be real.



#________________________________________________________________________________________________________________________________

# 4.a.ii Patterned resids, so trying to remodel with (1) neg binom model, and (2) LMM with variance modelled
#______________________________________________________________________________________________________________________________________________________

#==============
### Neg binom
#==============

mInt_log_NB <- glmer.nb(kereru ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise + (1|Observer)     #count conditions
                        + (1|PatchID), data=bird_Rescaled,
                        glmerControl(optimizer="bobyqa"))#, optCtrl=list(maxfun=1e9)))# model fails to converge, even with max iterations specified.
summary(mInt_log_NB)
AIC(mInt_log_NB,mInt_log)         #neg binom model is better fit, and highly sig for EVERYTHING.... ***but didn't converge properly so don't think can trust***

#==============
### LMM with variance modelling
#==============

#first, run LMM with no var modelling to check assumptions:
logbirds <- bird_Rescaled
logbirds$logkereru <- log(logbirds$kereru+0.01)
qplot(logkereru, data=logbirds)                  #creating log-transformed response var. [which is not close to normal!!!]

mInt_log_LMM <- lme(logkereru ~ logNBT*PatchPC5yr + Urban.YN + DaysSinceNov1 + MinSinceSunrise + Rain + Noise,# + (1|Observer)     #count conditions
                    random=~1|Observer + 1|PatchID, data=logbirds)
summary(mInt_log_LMM)
