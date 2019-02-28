### Data Aggregation for SHINY Empiric Prescribing App for ATLAS data only 

# Libraries
library(ggplot2); library(reshape2); library(magrittr); library(dplyr); library(fuzzyjoin)

### Where is the data? 
# setwd("")

### Read in the data
data <- read.csv("data/complete_Atlas_Reuse_Data.csv")#,fileEncoding = "UCS-2LE")
dim(data) # Original = 633820 x 111
colnames(data)

# List of resistancs
resistances <- c("Amikacin_I","Amoxycillin.clavulanate_I","Ampicillin_I","Azithromycin_I","Cefepime_I","Cefoxitin_I","Ceftazidime_I","Ceftriaxone_I","Clarithromycin_I",
                 "Clindamycin_I","Erythromycin_I","Imipenem_I","Levofloxacin_I","Linezolid_I","Meropenem_I","Metronidazole_I","Minocycline_I","Penicillin_I",
                 "Piperacillin.tazobactam_I","Tigecycline_I","Vancomycin_I","Ampicillin.sulbactam_I","Aztreonam_I","Aztreonam.avibactam_I","Cefixime_I",
                 "Ceftaroline_I","Ceftaroline.avibactam_I","Ceftazidime.avibactam_I","Ciprofloxacin_I","Colistin_I","ColistinP80_I","Daptomycin_I","Doripenem_I",
                 "Ertapenem_I","Gatifloxacin_I","Gentamicin_I","Moxifloxacin_I","Oxacillin_I","Quinupristin.dalfopristin_I",
                 "Sulbactam_I","Teicoplanin_I","Tetracycline_I","Trimethoprim.sulfa_I","Ceftolozane.tazobactam_I",colnames(data)[102:111])

##### Prevalence by 
# year / country / species / drug
datam <- data[,c("Species","Country","Year","Source",resistances)]
colnames(datam) <- sub("_I", "", colnames(datam)) # remove _I
colnames(datam) <- gsub(".", ' ', colnames(datam), fixed = T) # remove "."

##### Same resistance - so map one column to the other
# Ceftriaxone Cefuroxime
datam$Cefuroxime <- datam$Ceftriaxone
# Rename Amoxycillin.clavulanate to co-amoxiclav
colnames(datam)[colnames(datam)=="Amoxycillin clavulanate"] <- "Co-amoxiclav"
# Ampicillin and amoxicillin resistance shared - add in new column
datam$Amoxicillin <- datam$Ampicillin
# Flucloxacillin and oxacillin resistance shared - add in new column
datam$Flucloxacillin <- datam$Oxacillin

### Melt 
datam <- melt(datam, id.vars = c("Country","Year","Species","Source"))
datam <- datam[-which(datam$value == ""),] # Remove entries with no data
datam <- datam[complete.cases(datam),] # Remove NAs
dim(datam) # 14079156 x 6 

### Group by country / year / species / resistance
# i.e. not by source. Use for plotting
datam_g <- datam %>%
  group_by(Country, Year, Species, variable) %>%
  summarise(resistant=sum(value == "Resistant"), n=n(), 
            intermediate=sum(value == "Intermediate"), rate_r = (resistant + intermediate) / n,
            rate_i = intermediate / n,
            susceptible=sum(value == "Susceptible"),  rate_s = susceptible / n)

##### Species check
# length(unique(data$Species)) # 287
# length(unique(datam_g$Species)) # 287 = some data for all species 

###########################*************************************#####################################################################################
###########################*************************************#####################################################################################
#####********************** Sub-set needed for EMPIRIC DATA ANALYSIS ************************#######################################################
###########################*************************************#####################################################################################

####****** SOURCE & SYNDROME ******####
#### Which sources are in the currently included syndromes?
synd_map <- read.csv("data/syndrome_map.csv", stringsAsFactors = FALSE)[1:27,1:2]
colnames(synd_map)[1] = "data_source"
source_list <- unique(synd_map$data_source)
colnames(synd_map)[colnames(synd_map)=='data_source'] <- "Source"

# Merge with mapping also removes those sources not in the syndromes
# Use datam here, don't group by country yet
datam_source<-merge(datam, synd_map, by = "Source") 


####****** DRUGS ******####
#### Which resistances to drugs are needed to affect treatment of each syndrome?
drg_bkdwn <- read.csv("drug_bkdwn.csv")
drugs_here <- drg_bkdwn[,3:(dim(drg_bkdwn)[2])]
drug_list <- unique(drugs_here[!is.na(drugs_here)])
## Which drugs to include? Only those that are used to treat the syndromes included
drug_list <- drug_list[-which(drug_list == "")]

drg_bkdwn2 <- read.csv("drug_bkdwn_groups.csv")#,fileEncoding = "UCS-2LE")
drugs_here2 <- drg_bkdwn2[,3:(dim(drg_bkdwn2)[2])]
drug_list2 <- unique(drugs_here2[!is.na(drugs_here2)])
## Which drugs to include? Only those that are used to treat the syndromes included
drug_list2 <- drug_list2[-which(drug_list2 == "")]

drug_list = c(drug_list, drug_list2)

####****** SPECIES ******####
#### Which species contribute to each syndrome?
sp_bkdwn <- read.csv("sp_bkdwn.csv")

## remove those with no contribution 
sum_sp <- colSums(sp_bkdwn[,-1], na.rm = TRUE)
w<-c(1,which(sum_sp > 0)+1)
sp_bkdwn <- sp_bkdwn[,(w)] 

colnames(sp_bkdwn)[colnames(sp_bkdwn)=='?..syndrome'] <- "syndrome"
colnames(drg_bkdwn)[colnames(drg_bkdwn)=='?..syndrome'] <- "syndrome"

rowSums(sp_bkdwn[,-1]) # check all 100
msp_bkdwn <- melt(sp_bkdwn, id.vars = c("syndrome"))
msp_bkdwn$value <- msp_bkdwn$value/100
colnames(msp_bkdwn) <- c("syndrome","Variable","Prop_Syn")
msp_bkdwn_list <- msp_bkdwn[which(msp_bkdwn$Prop_Syn > 0),]
msp_bkdwn_list$Variable <- gsub(".", ' ', msp_bkdwn_list$Variable, fixed = T)
## Need "Streptococcus, viridans group"
g <- grep("viri",msp_bkdwn_list$Variable)
msp_bkdwn_list[g,"Variable"] <- "Streptococcus, viridans group"

msp_bkdwn_list$Species = msp_bkdwn_list$Variable
msp_bkdwn_list <- msp_bkdwn_list[,c("syndrome","Species","Prop_Syn")]
## Which species to include? These are the ones that cause the syndromes included 
bug_list <- unique(msp_bkdwn_list$Species)

####****** SUBSET DATA ******####
### Subset of data that is in drug, bug and source list
datam_source_bug <- filter(datam_source, Species %in% bug_list)
datam_source_bug_drug <- filter(datam_source_bug, variable %in% drug_list)

dim(datam) # 14079156
dim(datam_source) # 10956303
dim(datam_source_bug) # 9049994
dim(datam_source_bug_drug) # 4745604

####****** YEAR *******####
datam_source_bug_drug17 <- as.data.frame(subset(datam_source_bug_drug, Year == 2017)) # ONLY 2017?

###******** SPECIES ***####
length(unique(bug_list)) # need 19
length(unique(datam_source_bug_drug$Species)) # have 19
length(unique(datam_source_bug_drug17$Species)) # in 2017 - drops to 17, lose two! One is only in puc. OK to ignore. 
# add in mock line to keep in bar plot of app
dl <- dim(datam_source_bug_drug17)[1] # 435557
datam_source_bug_drug17[dl+1, ] <- datam_source_bug_drug17[dl, ]
datam_source_bug_drug17[dl+1,c("Species","value","syndrome")] <- c("Streptococcus, viridans group","Susceptible","cap") 

length(unique(datam_source_bug_drug17$Species)) 
dim(datam_source_bug_drug17) # 468675 with extra line

####****** GROUP ******####
### Group by country / year (left in but should all be 2017) / species / syndrome / antibiotic
# No longer do by source - don't need, just syndrome
# Need for sub-setting below
datam_g_s <- datam_source_bug_drug17 %>%
  dplyr::group_by(Country, Year, Species, syndrome, variable) %>%
  dplyr::summarise(resistant=sum(value == "Resistant"), n=n(), 
                   intermediate=sum(value == "Intermediate"), 
                   susceptible=sum(value == "Susceptible"),  
                   rate_r = (resistant + intermediate) / n,
                   rate_i = intermediate / n,
                   rate_s = susceptible / n)


####****** MERGE ******####
### Merge with contribution 
datam_gs <- merge(datam_g_s, msp_bkdwn_list, by = c("Species","syndrome")) 
## lose "Streptococcus dysgalactiae"

write.csv(datam_gs, "data_bug_drug_synd_all.csv")

### linkage to cost data ##################################################################################
cost.data <- read.csv("data/MSH_price_access.csv")

# Take from drug list
antib <- unique(melt(drg_bkdwn[,c(1:4,8,9,13)], id.vars = c("syndrome","Age"))[,"value"], na.rm = TRUE)
antib <- antib[-which(antib == "")]
base.df <- data.frame(antib)

cost.data$antib <- cost.data$Variable
## with test set of prelimary list of drugs; no max set the amount is 18,640
## with 0.375 set this goes down to 1,629
# drugs capitalised in both sets so no need to set ignore.case=TRUE

full.data <- stringdist_full_join(base.df, cost.data, by ="antib",
                                  
                                  method = "jw", distance_col = "fuzzy.no",
                                  
                                  max_dist = 0.375)

## split into perfectly matched and non-perfectly matched
## the non-perfectly matched were searched by included antibiotic and a list of reasonable matches
# were made and coded to be kept
# matched.data = 45
matched.data <- subset(full.data, fuzzy.no==0)
# query data = 1106
query.data <- subset(full.data, fuzzy.no!=0)
# NA = 729 (855+45+729 = 1,629)

## scrolled through to see matches (searched for those with no exact match)
query.data <- query.data[((query.data$antib.x=="Clindamycin" & query.data$antib.y=="Clindamycin (Base)")|
                            (query.data$antib.x=="Penicillin" & query.data$antib.y=="Penicillamine") |
                            (query.data$antib.x=="Penicillin" & query.data$antib.y=="Penicillin G")|
                            (query.data$antib.x=="Penicillin" & query.data$antib.y=="Penicillin, Benzyl")|
                            (query.data$antib.x=="Piperacillin tazobactam" &
                               query.data$antib.y=="Piperacillin+Tazobactam")|
                            (query.data$antib.y=="Amoxicillin/Clavulanic Acid")|
                            (query.data$antib.x=="Trimethoprim" & query.data$antib.y=="Sulfamethoxazole/Trimethoprim")
),]

w<-which(query.data$antib.y == "Amoxicillin/Clavulanic Acid")
query.data[w,"antib.x"] <- "Co-amoxiclav"

## combine datasets back together
com.data <- rbind(matched.data, query.data)

## add back in antibiotics that didn't have a cost
com.data[nrow(com.data) + 1,] = c("Fosfomycin" ,rep(NA,(ncol(com.data)-1)))

## inflation to 2017
inflation <- read.csv('data/WB_inflation.csv')
i.2015 <- (inflation$X2015..YR2015.[1])/100
i.2016 <- (inflation$X2016..YR2016.[1])/100
i.2017 <- (inflation$X2017..YR2017.[1])/100

## MSH prices in July 2015 prices, need to inflate to end of 2015 then onwards
i.2015.2 <- ((1+i.2015)^(1/2))-1
com.data$c.inf <- as.numeric(com.data$Supplier.Median..US..)+(as.numeric(com.data$Supplier.Median..US..)*i.2015.2)
com.data$c.inf <- (com.data$c.inf+(com.data$c.inf*i.2016))
com.data$c.inf <- (com.data$c.inf+(com.data$c.inf*i.2017))
com.data$c.inf <- round(com.data$c.inf,3)

# now prices are in Jan 2018 prices

## rename columns and merge some;
com.data$Antibiotic.ATLAS <- com.data$antib.x
com.data$Antibiotic.AR.IA <- com.data$antib.y
com.data$Formulation <- paste(com.data$Dosage.Form, com.data$Strength)
com.data$Cost <- paste(com.data$c.inf, "per", com.data$Comparison.Unit)
com.data$Access <- com.data$WHO.Status

## linking in Q's grouping data csv
drg_grps <- read.csv("data/drug_groups.csv")
drg_grps <- subset(drg_grps, is.single==TRUE)
colnames(drg_grps)[colnames(drg_grps)=='antibiotic'] <- "Antibiotic.ATLAS"

com.grp <- merge(com.data, drg_grps, by=c("Antibiotic.ATLAS"), all=TRUE)
# drop the ones from drg_grps not in com.data
com.grp <- subset(com.grp, !is.na(com.grp$antib.x))

### adding in AWaRe list
## hand entered and changed first letter to capital (to match ATLAS formatting)
## added in Q's grouping spellings from drg_groups.csv
# *Watch group antibiotics included in the EML/EMLc only for specific, limited indications
# i.e. some are in both access and watch for different indications

v.access <- c("Amoxicillin", "Cefotaxime", "Amikacin", "Gentamicin",
                            "Amoxicillin/Clavulanic Acid", "Ceftriaxone", "Azithromycin", "Metronidazole",
            "Ampicillin", "Cloxacillin", "Chloramphenicol", "Nitrofurantoin",
              "Benzathine benzylpenicillin", "phenoxymethylpenicillin", "Ciprofloxacin",
            "Spectinomycin", "Benzylpenicillin", "piperacillin tazobactam", "Clarithromycin",
              "Piperacillin+Tazobactam", "Cefalexin", "Procaine Benzyl Penicillin",
              "Clindamycin", "Vancomycin", "Cefazolin", "Meropenem",
              "Doxycycline", "Cefixime","sulfamethoxazole+trimethoprim","Sulfamethoxazole/Trimethoprim",
              "Trimethoprim-sulfa", "")

v.watch <- c("Quinolones", "Fluoroquinolones","Ciprofloxacin",
             "Levofloxacin", "Moxifloxacin", "Norfloxacin",
             "3rd-generation Cephalosporins", "Cefixime", "Ceftriaxone",
             "Cefotaxime", "Ceftazidime","Macrolides", "Azithromycin",
             "Clarithromycin", "Erythromycin", "Glycopeptides", "Teicoplanin",
             "Vancomycin", "Piperacillin+Tazobactam", "Carbapenems", "meropenem",
             "Imipenem", "Cilastatin", "Penems", "Faropenem","Cephalosporins (3rd gen)")

v.reserve <- c("Aztreonam","Fosfomycin","4th-generation Cephalosporins","cefepime",
               "Oxazolidinones","linezolid", "5-th generation Cephalosporins",
               "Ceftaroline","Tigecycline","Polymyxins","Polymyxin B",
               "Colistin","Daptomycin","Cephalosporins (4th gen)","Cephalosporins (5th gen)")

## code NA groups differently to avoid all getting dumped in "Access"
com.grp$group <- as.character(com.grp$group)
com.grp$group[is.na(com.grp$group)] <- "Not Grouped"

## list if in watch group
com.grp$AWaRe <- NA

com.grp$AWaRe[com.grp$Antibiotic.ATLAS %in% v.access|
                com.grp$Antibiotic.AR.IA %in% v.access|
                com.grp$group %in% v.access] <- c("Access")

com.grp$AWaRe[com.grp$Antibiotic.ATLAS %in% v.watch|
                com.grp$Antibiotic.AR.IA %in% v.watch|
                com.grp$group %in% v.watch] <- c("Watch")

com.grp$AWaRe[com.grp$Antibiotic.ATLAS %in% v.reserve|
                com.grp$Antibiotic.AR.IA %in% v.reserve|
                com.grp$group %in% v.reserve] <- c("Reserve")

## adapt WHO Access terms
com.grp$WHO.Access[com.grp$Access=="E"] <- c("Essential Medicine - Core List")
com.grp$WHO.Access[com.grp$Access=="EP"] <- c("Essential Medicine - Core List")
com.grp$WHO.Access[com.grp$Access=="C"] <- c("Essential Medicine - Complementary List")
com.grp$WHO.Access[com.grp$Access=="P"] <- c("Essential Medicine in different dosage/strength")
com.grp$WHO.Access[com.grp$Access=="T"] <- c("Potential alternative to an Essential Medicine")
com.grp$WHO.Access[com.grp$Access=="N"] <- c("Not on the list")
com.grp$WHO.Access[is.na(com.grp$Access)] <- c("Not on the list")

com.grp$AWaRe[is.na(com.grp$AWaRe)] <- c("Not on the list")
econ <- com.grp[ ,c("Antibiotic.ATLAS", "Antibiotic.AR.IA", "group",
                    "Cost","Formulation","WHO.Access", "AWaRe")]

write.csv(econ, "econ_data.csv")

