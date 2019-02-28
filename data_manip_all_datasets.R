###Data aggregation for SHINY Empiric Prescribing App for ECDC, ResistanceMap and GLASS data

#Libraries
library(reshape2); library(magrittr); library(dplyr)

#Where is the data?
#setwd('')

####****** SOURCES AND SYNDROMES ******####

synd_map <- read.csv("syndrome_map.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")[1:22,1:2]
source_list <- unique(synd_map$data_source)
## Need to alter data to include the source to syndrome mapping
colnames(synd_map)[colnames(synd_map)=="data_source"] <- "Source"


####****** DRUGS ******####

drg_bkdwn <- read.csv("drug_bkdwn.csv", fileEncoding="UTF-8-BOM")
drugs_here <- drg_bkdwn[,3:(dim(drg_bkdwn)[2])]
drug_list <- unique(drugs_here[!is.na(drugs_here)])
## Which drugs to include? Only those that are used to treat the syndromes included
drug_list <- drug_list[-which(drug_list == "")]

drg_bkdwn2 <- read.csv("drug_bkdwn_groups.csv")
drugs_here2 <- drg_bkdwn2[,3:(dim(drg_bkdwn2)[2])]
drug_list2 <- unique(drugs_here2[!is.na(drugs_here2)])
## Which drugs to include? Only those that are used to treat the syndromes included
drug_list2 <- drug_list2[-which(drug_list2 == "")]

drug_list = c(drug_list, drug_list2)


####****** SPECIES ******####

#Which species contribute to each syndrome?
sp_bkdwn <- read.csv("sp_bkdwn.csv", fileEncoding="UTF-8-BOM")

rowSums(sp_bkdwn[,-1]) # check all 100 
msp_bkdwn <- melt(sp_bkdwn, id.vars = c("syndrome"))
msp_bkdwn$value <- msp_bkdwn$value/100
colnames(msp_bkdwn) <- c("syndrome","Variable","Prop_Syn")
msp_bkdwn_list <- msp_bkdwn[which(msp_bkdwn$Prop_Syn > 0),]
msp_bkdwn_list$Variable <- gsub(".", ' ', msp_bkdwn_list$Variable, fixed = T)
msp_bkdwn_list$Species = msp_bkdwn_list$Variable
msp_bkdwn_list <- msp_bkdwn_list[,c("syndrome","Species","Prop_Syn")]

## Which species to include? These are the ones that cause the syndromes included 
bug_list <- unique(msp_bkdwn_list$Species)



####****** Aggregate ECDC data ******####

ecdc_data = read.csv("ECDC_data.csv")

ecdc_data = ecdc_data[which(ecdc_data$Time == 2017),] #only include 2017

#Reformat bacteria and antibiotic names
temp = data.frame(do.call("rbind", strsplit(as.character(ecdc_data$Population), "|", fixed=T))) #split bacteria and antibiotic names
ecdc_data$Population = temp$X1 #reassign bacteria names
ecdc_data$Antibiotic = temp$X2 #reassign antibiotic names
rm(temp)
#drop some resistances (otherwise mismatch below, because don't give I and R proportion):
ecdc_data = ecdc_data[-grep("Combined", ecdc_data$Antibiotic),]
ecdc_data = ecdc_data[-grep("MRSA", ecdc_data$Antibiotic),]
ecdc_data = ecdc_data[-grep("High-level", ecdc_data$Antibiotic),]

#only take rows looking at fraction of resistant isolates and total number of tested isolates:
ecdc_data = ecdc_data[union(which(ecdc_data$Indicator== "Non-susceptible (I and R) isolates proportion  "),
                            which(ecdc_data$Indicator== "Total tested isolates")),]

#collapse dataset to have resistance proportion and sample size on the same rows:
#this should be TRUE (ie the first half of the dataset is in the same order as the second):
all.equal(ecdc_data$Population[1:690], ecdc_data$Population[691:1380])
ecdc_data = as.data.frame(cbind(ecdc_data[1:(dim(ecdc_data)[1]/2),], 
                                sample_size=ecdc_data$NumValue[(dim(ecdc_data)[1]/2+1):dim(ecdc_data)[1]]))

#reformat dataframe keeping only variables of interest:
ecdc_data = data.frame(Species = ecdc_data$Population, 
                       Country = ecdc_data$RegionName, 
                       Year = ecdc_data$Time, 
                       variable = ecdc_data$Antibiotic,
                       n = as.numeric(as.character(ecdc_data$sample_size)),
                       rate_r = as.numeric(as.character(ecdc_data$NumValue))/100)
ecdc_data = ecdc_data[-(which(is.na(ecdc_data$rate_r))),]

#ecdc data is only from invasive isolates, ie sepsis:
ecdc_data$syndrome = "sep"


#assumption to match names of species in sp_bkdwn:
ecdc_data$Species = as.character(ecdc_data$Species)
ecdc_data$Species[which(ecdc_data$Species == "Acinetobacter spp.")] = "Acinetobacter  non speciated"


#Filter and merge:
ecdc_data <- filter(ecdc_data, Species %in% bug_list)
ecdc_data <- filter(ecdc_data, variable %in% drug_list)
ecdc_data <- merge(ecdc_data, msp_bkdwn_list, by = c("Species","syndrome"))


write.csv(ecdc_data, "ECDC_shinyapp.csv", row.names = FALSE)



####****** Aggregate ResistanceMap data ******####

resist_data = read.csv("resist_data.csv")
resist_data = resist_data[resist_data$year == 2015,] #2015 data only (2016 and 2017 data are scarce)

#reformat dataframe keeping only variables of interest:
resist_data = data.frame(Species = resist_data$Organism,
                         Country = resist_data$CountryName,
                         Year = resist_data$year,
                         variable = resist_data$Antimicrobial, 
                         n = resist_data$IsolatesTested,
                         rate_r = resist_data$value/100)

#assumption to match names of species in sp_bkdwn:
resist_data$Species = as.character(resist_data$Species)
resist_data$Species[which(resist_data$Species == "Enterobacter aerogenes/cloacae")] = "Enterobacter aerogenes"

#resist data is only from invasive isolates, ie sepsis:
resist_data$syndrome = "sep"

#Filter and merge
resist_data <- filter(resist_data, Species %in% bug_list)
resist_data <- filter(resist_data, variable %in% drug_list)
resist_data <- merge(resist_data, msp_bkdwn_list, by = c("Species","syndrome"))


write.csv(resist_data, "RESIST_shinyapp.csv", row.names = FALSE)



####****** Aggregate GLASS data ******####

glass = read.csv("glass_data.csv")

#Add Aminopenicillins goup
glass2 = glass[which(glass$ï..Antibiotic=="Ampicillin"),]
glass2$ï..Antibiotic = as.character(glass2$ï..Antibiotic)
glass2$ï..Antibiotic = "Aminopenicillins"

glass = rbind(glass,glass2)
rm(glass2)

#reformat dataframe keeping only variables of interest:
glass_data = as.data.frame(cbind(as.character(glass$Pathogen),
                                 as.character(glass$Country),
                                 rep(2017,length(glass$Country)),
                                 as.character(glass$ï..Antibiotic),
                                 as.character(glass$tested),
                                 as.character(glass$Resistant),
                                 as.character(glass$Intermediate),
                                 as.character(glass$SPECIMEN..Sample.)))
colnames(glass_data) = c("Species", "Country", "Year", "variable", "n", "resistant","intermediate","syndrome")
glass_data$resistant = as.numeric(as.character(glass_data$resistant))
glass_data$intermediate = as.numeric(as.character(glass_data$intermediate))
glass_data$n = as.numeric(as.character(glass_data$n))

#Regroup all age groups
glass_data = glass_data %>%
  group_by(Species, Country, Year, variable, syndrome) %>%
  summarise(n = sum(n), resistant=sum(resistant), intermediate=sum(intermediate))

glass_data$rate_r = (glass_data$resistant+glass_data$intermediate)/glass_data$n

glass_data = glass_data[-which(is.na(glass_data$rate_r)),]

#match code based on synd_map.csv
glass_data$syndrome = unlist(sapply(glass_data$syndrome, FUN=function(x){
  if(x == "blood") x="sep"
  else if(x == "urine") x="c.uti"
  else if(x == "genital") x="puc"
  else if(x == "stool") x="sto"
}))

#Remove isolates from stool (not included in our symptoms)
glass_data = glass_data[-(which(glass_data$syndrome == "sto")),]


#Filter and merge
glass_data <- filter(glass_data, Species %in% bug_list)
glass_data <- filter(glass_data, variable %in% drug_list)
glass_data <- merge(glass_data, msp_bkdwn_list, by = c("Species","syndrome"))

#Clean up
glass_data = glass_data[,-c(7,8)]


write.csv(glass_data, "GLASS_shinyapp.csv", row.names = F)


####****** Combine all datasets into one ******####

atlas_data = read.csv("data_bug_drug_synd_all.csv")[-1]

atlas_data = data.frame(Species = atlas_data$Species,
                        syndrome = atlas_data$syndrome,
                        Country = atlas_data$Country,
                        Year = atlas_data$Year,
                        variable = atlas_data$variable,
                        n = atlas_data$n,
                        rate_r = atlas_data$rate_r,
                        Prop_Syn = atlas_data$Prop_Syn,
                        Dataset = "ATLAS")

#Identify the origin of the datasets
ecdc_data$Dataset = "ECDC"
resist_data$Dataset = "RESIST"
glass_data$Dataset = "GLASS"

#Combine
all_datasets = rbind(atlas_data,ecdc_data,resist_data,glass_data)

write.csv(all_datasets, "all_datasets.csv", row.names = F)
