## Shiny app for Empiric Prescribing calculations
# National level, with global data

# Required libraries
library(shiny)
library(ggplot2)
library(plyr)
library(dplyr)
library(rworldmap)
library(RColorBrewer)
library(reshape2)
library(DT)

###********************************************************************************************************************************************************************#####
#### REQUIRED DATA AND FUNCTIONS ####
###********************************************************************************************************************************************************************#####

datam_map_all <- read.csv("all_datasets.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM") # aggregated data
drg_bkdwn <- read.csv("drug_bkdwn.csv", fileEncoding="UTF-8-BOM") # baseline empiric therapy recommendations
drg_bkdwn_group <- read.csv("drug_bkdwn_groups.csv", fileEncoding="UTF-8-BOM") # baseline empiric therapy recommendations with drug groupings instead of individual drugs
sp_bkdwn <- read.csv("sp_bkdwn.csv", fileEncoding="UTF-8-BOM") # baseline contributing bacteria distributions
sp_all <- read.csv("sp_all.csv", fileEncoding="UTF-8-BOM") # has just the list of species to include
econ <- read.csv("econ_data.csv", fileEncoding="UTF-8-BOM")[,-1] # economic data: WHO EML and AWaRE

# Function to keep NAs if all NA
suma = function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)

# To map just ECDC data countries
ecdc_countries <- c("Austria","Belgium","Bulgaria","Croatia","Cyprus","Czech Rep.","Denmark","Estonia","Finland","France",
                    "Germany","Greece","Hungary","Iceland","Ireland","Italy","Latvia",
                    "Lithuania","Luxembourg","Malta","Netherlands","Norway","Poland",
                    "Portugal","Romania","Slovakia","Slovenia","Spain","Sweden","United Kingdom")

###********************************************************************************************************************************************************************#####
#### SHINY CODE ####
###********************************************************************************************************************************************************************#####
ui <- fluidPage(  
  
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  
  titlePanel("Empiric Prescribing"),
  
  sidebarPanel(h1("Syndrome"),
               
               # Syndrome input
               selectInput("variable", "Choose syndrome:",
                           c("Sepsis" = "sep",
                             "Complicated UTI" = "c.uti",
                             "Hospital acquired pneumonia" = "hap",
                             "Community acquired pneumonia" = "cap",
                             "Cellulitis / skin abcess" = "skin",
                             #"Purulent urthritis / cervicitis" = "puc",
                             "Upper respiratory tract infection" = "urti",
                             "Bacterial meningitis" = "meng",
                             "Septic arthritis" = "sepa")),
               
               # Age input
               selectInput("age", "Age", 
                           c("Adult" = "adult",
                             "Child" = "child",
                             "Neonate" = "infant")),
               
               hr(), 
               
               h4("Optional inputs"),
               
               # Barchart of bug composition
               plotOutput("inputplot"),
               
               # Resistance cutoff input
               numericInput("res_cut", em("Resistance cutoff (%)"), 15, min = 0, max = 100),
               
               # Which data? 
               checkboxGroupInput("checkGroup", label = em("Data to include"),
                                  c("ATLAS", "GLASS", "RESISTANCEMAP" = "RESIST",
                                    "ECDC"),
                                  selected = "ATLAS")
  ),
  
  mainPanel(
    
    # Output: Tabset w/ plot, summary, and table
    tabsetPanel(type = "tabs", 
                
                # First tab: bacterial distribution 
                tabPanel("Input: Contributing pathogen distribution", 
                         textOutput("sliders_total"),
                         verbatimTextOutput("text_warning_glas"), # warnings if use data other than ATLAS - FUTURE area to update: e.g. grey out syndromes in input potential
                         verbatimTextOutput("text_warning_resm"),
                         verbatimTextOutput("text_warning_ecdc"),
                         h2("\n"),
                         uiOutput("sliders"),
                         h6("Proportions must sum to 100")
                ),
                # Second tab: input of empiric therapy and economics
                tabPanel("Input: Empiric therapy guidelines",
                         tableOutput("empirictherapy"),
                         h2(""),
                         h4("WHO Essential Medicines List (EML) and AWaRE data"),
                         h6("Costs: Costs presented are in Dec. 2017 USD."),
                         h6("WHO EML: Based on The International Medical Products Price Guide (2015)."),
                         h6("AWaRe: Based on the 2017 World Health Organisation Essential Medicines List, for those in multiple categories for different indications the most severe AWaRe grouping was assigned"),
                         dataTableOutput('econ')),
                # Third tab: map of resistance prevalence to first line
                tabPanel("Output: Map of resistance prevalence",
                         h4("Maps of prevalence of resistance to drugs in first line therapy"),
                         textOutput("text1"),
                         plotOutput("coolplot"),
                         h6("Countries shaded grey had no data. Countries with hatching had less than 10 isolates to inform prevalence."),
                         hr(),
                         textOutput("text2"),
                         plotOutput("coolplot2")),
                # Fourth tab: Underlying data table
                tabPanel("Output: Table of data", 
                         dataTableOutput("results")),
                # Fifth tab: Map and table of therapy recommendations
                tabPanel("Output: Therapy recommendations",
                         verbatimTextOutput("text_warning_glas2"), # data warnings 
                         verbatimTextOutput("text_warning_resm2"),
                         verbatimTextOutput("text_warning_ecdc2"),
                         # textOutput("text1") #,
                         plotOutput("trafficplot"),
                         h6("Linked to Resistance cutoff value. Countries shaded grey had no data.
                            If no image, then there is not enough data for recommendations"),
                         h4("Recommendations for therapy"),
                         tableOutput("recc")
                )
    )
  )
)

###********************************************************************************************************************************************************************#####
#### FUNCTIONS TO POWER SHINY ####
###********************************************************************************************************************************************************************#####

server <- function(input, output) {
  
  ### ****************************************************************************************************************************************###
  # MODEL - takes in original data and does many of the basic manipulations to get the data in the correct format for the rest of the analysis
  ### ****************************************************************************************************************************************###
  
  # Reactive dependencies - if these change then MODEL will run again and update values
  xxchange <- reactive({
    paste(input$variable, input$age, input$checkGroup, 
          input$range1, input$range2, input$range3, input$range4, 
          input$range5, input$range6, input$range7, input$range8,
          input$range9, input$range10, input$range11, input$range12,
          input$range13, input$range14, input$range15, input$range16,
          input$range17, input$range18)
  })
  
  # New data subset
  model <- eventReactive(xxchange(), {
    
    # Which dataset? 
    data2 <- datam_map_all[datam_map_all$Dataset %in% input$checkGroup,]
    
    # needed function
    suma = function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)
    
    # Take averages over datasets in checkgroup
    datam_map <- data2 %>%
      dplyr::group_by(Country,syndrome,Species,variable) %>%
      dplyr::summarise(rate_r=mean(rate_r), Prop_Syn = mean(Prop_Syn),n = sum(n))
    
    # Which syndrome? 
    synd <- input$variable
    # What age?
    age <- input$age
    
    # Subset of data for this syndrome
    datam_map2 <- datam_map[which(datam_map$syndrome == synd),]
    
    # If not ATLAS then have to use grouped drugs
    if(length(input$checkGroup) > 1 || input$checkGroup != "ATLAS"){
      drg_bkdwn <- drg_bkdwn_group
    }
    
    # Drugs for this synd
    d11 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"First.line.1st.drug"])
    d12 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"First.line.2nd.drug"])
    d21 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Second.line.1st.drug"])
    d22 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Second.line.2nd.drug"])
    d3 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Third.line.1st.drug"])
    
    pvars <- as.character(sp_all[,1]) # all species to be included
    
    # Species for this syndrome  
    sp_b <- melt(sp_bkdwn[which(sp_bkdwn[,1] == synd),],id.vars = "syndrome")
    sp_b$variable <- gsub(".", ' ', sp_b$variable, fixed = T)
    ## Need "Streptococcus, viridans group"
    g <- grep("viri",sp_b$variable)
    sp_b[g,"variable"] <- "Streptococcus, viridans group"
    
    # Values for the proportions are those given by the sliders 
    values_slid <- filter(sp_b, variable %in% pvars)
    
    pvars_use <- pvars
    names <- c()
    for(i in 1:length(pvars_use)){names <- c(names,paste0("range", i))}
    
    values <- c()
    for(i in 1:length(pvars_use)){values <- c(values,input[[names[i]]])}
    
    # Data for barchart of proportions
    bar_plot_data <- as.data.frame(cbind(rep(1, length(values)),values))
    colnames(bar_plot_data) <- c("x","value")
    bar_plot_data$Species <- values_slid$variable
    
    ## New levels from input bars in bar_plot_data
    # merge with resistance data 
    datam_map_new <- merge(datam_map2, bar_plot_data, by = "Species")
    # update with new values
    datam_map_new$Prop_Syn <- datam_map_new$value/100 
    datam_map_new$res_prop <- 100*datam_map_new$Prop_Syn * datam_map_new$rate_r
    
    ## Group by country
    datam_c <- datam_map_new %>%
      dplyr::group_by(Country, syndrome, variable) %>%
      dplyr::summarise(res_perc=sum(res_prop), n = sum(n))
    
    # Add in column - if less than 10 isolates
    w<-which(datam_c$n > 10)
    datam_c$n10 <- 100
    datam_c[w,"n10"] <- 40
    
    # Data for drug one and potentially two in first line therapy
    datam_map_new_d11 <- as.data.frame(datam_c[which(datam_c$variable == d11),])
    datam_map_new_d12 <- as.data.frame(datam_c[which(datam_c$variable == d12),])
    
    ### Map resistance proportion in this syndrome
    europe_11 <- 0;europe_12 <- 0
    mapped_data_d11 <-0; mapped_data_d12 <- 0;
    new_world_11 <-0; new_world_12 <- 0;
    
    # indicators for printing plots - only under certain conditions - FUTURE area to update: more sophisticated warnings around plot output
    p1 <- 0; p2 <- 0; p3 <- 0; 
    
    # Only if data will map data be generated and plotted
    if(dim(datam_map_new_d11)[1] > 1){
      p1 <- 1; # PRINT
      mapped_data_d11 <- joinCountryData2Map(datam_map_new_d11, joinCode = "NAME", nameJoinColumn = "Country")
      
      new_world_11 <- subset(mapped_data_d11, continent != "Antarctica") # remove Antartica
      if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){europe_11 <- subset(mapped_data_d11, NAME%in%ecdc_countries)}} # if only ECDC data
    
    if(dim(datam_map_new_d12)[1] > 1){
      p2 <- 1; # PRINT
      mapped_data_d12 <- joinCountryData2Map(datam_map_new_d12, joinCode = "NAME", nameJoinColumn = "Country")
      
      new_world_12 <- subset(mapped_data_d12, continent != "Antarctica")
      if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){europe_12 <- subset(mapped_data_d12, NAME%in%ecdc_countries)}}
    
    ## CHECK IF ANY data
    datam_any <- which(unique(datam_c$variable) %in% c(d11,d12,d21,d22,d3))
    if(length(datam_any) > 0){p3 <- 1} # DON"T PRINT ANYTHING IF NO DATA - FUTURE area to update: earlier warnings / prevent user going further 
    
    # return all object as a list
    list(datam_map = datam_map, datam_c = datam_c, pvars = pvars, 
         d11 = d11, d12 = d12, d21 = d21, d22 = d22, d3 = d3,
         mapped_data_d11 = new_world_11, europe_map_11 = europe_11,
         mapped_data_d12 = new_world_12, europe_map_12 = europe_12, p1 = p1, p2 = p2, p3 = p3)
    
  }
  )
  
  
  ### ****************************************************************************************************************************************###
  # TREAT - formulates the treatment decisions 
  ### ****************************************************************************************************************************************###
  
  # Reactive dependencies - if these change then TREAT will run again 
  xxchanget <- reactive({
    paste(model()$datam_map, model()$datam_c, model()$d11,model()$d12, model()$d21, model()$d22, model()$d3, input$res_cut,
          input$range1, input$range2, input$range3, input$range4, 
          input$range5, input$range6, input$range7, input$range8,
          input$range9, input$range10, input$range11, input$range12,
          input$range13, input$range14, input$range15, input$range16,
          input$range17, input$range18)
  })
  
  # Treatment data subset
  treat <- eventReactive(xxchanget(), {
    
    ### Data
    datam_map_new <- model()$datam_c
    
    ###### THERAPY DECISIONS #########
    ### Resistance cut off
    res_cut <- input$res_cut
    
    ### Column for each therapy 
    datam_map_new$first <- NA
    datam_map_new$second <- NA
    datam_map_new$third <- NA
    
    ######## where should use first line? ###########################
    # first line therapy for this synd = d11 / d12 combination 
    if(model()$p1 == 1||model()$p2 == 1){ # if data
      # first line therapy for this synd = d11 / d12 combination
      
      d11 <- model()$d11
      d12 <- model()$d12
      
      # Combination therapy?
      if(d12 != ""){
        ## Subset data to be just those with both drugs
        dd <- datam_map_new
        dd <- dd[which(dd$variable %in% c(d11,d12)),]
        dd <- dd %>% dplyr::group_by(Country) %>% dplyr::summarise(nn = sum(res_perc < res_cut), 
                                                                   count = n(), cisol1 = sum(n), nh = sum(res_perc > res_cut)) 
        if(dim(dd)[1]>0){
          datam_map_new <- merge(dd,datam_map_new, by = "Country",all.y=TRUE)}
        
        # If both lower than threshold then... use first line therapy
        w<-which(datam_map_new$nn == 2)
        datam_map_new[w,"first"] <- 1
        
        # If data on both drugs and not low resistance then don't use first line
        w2<-setdiff(which(datam_map_new$count == 2),w)
        datam_map_new[w2,"first"] <- 0
        
        # If only data for one or the other then leave as NA
        # UNLESS one drug > res_cut: then can't use first
        w3 <- which(datam_map_new$nh > 0)
        datam_map_new[w3,"first"] <- 0
        
      } else { # if only one drug
        dd <- datam_map_new
        dd <- dd[which(dd$variable %in% c(d11)),]
        dd <- dd %>% dplyr::group_by(Country) %>% dplyr::summarise(nn = sum(res_perc < res_cut), count = n(), cisol1 = sum(n), nh = sum(res_perc > res_cut)) 
        if(dim(dd)[1]>0){
          datam_map_new <- merge(dd,datam_map_new, by = "Country",all.y=TRUE)}
        
        ### Where resistance high? 
        wr <- which(datam_map_new$res_perc >  res_cut) # resistance higher than res_cut input (20%)
        wd1 <- which(datam_map_new$variable == d11)
        w1n <- intersect(wd1, wr) # where information and resistance to d1 > res_cut
        w1y <- setdiff(wd1, wr) # where information and resistance to d1 < res_cut
        
        datam_map_new[w1y,"first"] <- 1 # Use first line, resistance low
        datam_map_new[w1n,"first"] <- 0 # Don't use first line 
      }
      
      ## remove nn and count
      if("nn" %in% colnames(datam_map_new)){
        w<-which(colnames(datam_map_new) == "nn"); w2<-which(colnames(datam_map_new) == "count"); w3<-which(colnames(datam_map_new) == "nh")
        datam_map_new <- datam_map_new[,-c(w,w2,w3)]}
      
      ######## where should use second line? ###########################
      # where should use second line, if available?
      # second line therapy for this synd = d2 ### NOT DONE FOR COMBINATION THERAPY
      d21 <- model()$d21 #as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),8])
      d22 <- model()$d22 #as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),8])
      
      # Combination therapy?
      if(d22 != ""){
        ## Subset data to be just those with both drugs
        dd <- datam_map_new
        dd <- dd[which(dd$variable %in% c(d21,d22)),]
        dd <- dd %>% dplyr::group_by(Country) %>% dplyr::summarise(nn = sum(res_perc < res_cut), count = n(), cisol2 = sum(n), nh = sum(res_perc > res_cut)) 
        if(dim(dd)[1]>0){
          datam_map_new <- merge(dd,datam_map_new, by = "Country",all.y=TRUE)}
        w<-which(datam_map_new$nn == 2) # which countries have both res less than res_cut
        datam_map_new[w,"second"] <- 1
        
        # If data on both drugs and not low resistance then don't use
        # first says which has data on both, second which has resistance < res_cut
        # setdiff takes whats in first and not in second
        w2<-setdiff(which(datam_map_new$count == 2),w) 
        datam_map_new[w2,"second"] <- 0
        
        # UNLESS one drug > res_cut: then can't use second
        w3 <- which(datam_map_new$nh > 0)
        datam_map_new[w3,"second"] <- 0
        
      } else { # if only one drug
        dd <- datam_map_new
        dd <- dd[which(dd$variable %in% c(d21)),]
        dd <- dd %>% dplyr::group_by(Country) %>% dplyr::summarise(nn = sum(res_perc < res_cut), count = n(), cisol2 = sum(n), nh = sum(res_perc > res_cut)) 
        if(dim(dd)[1]>0){
          datam_map_new <- merge(dd,datam_map_new, by = "Country",all.y=TRUE)}
        
        ### Where resistance high? 
        wr <- which(datam_map_new$res_perc >  res_cut) # resistance higher than res_cut input (20%)
        wd2 <- which(datam_map_new$variable == d21)
        w2n <- intersect(wd2, wr) # where information and resistance to d2 > res_cut
        w2y <- setdiff(wd2, wr) # where information and resistance to d2 < res_cut
        
        datam_map_new[w2y,"second"] <- 1 # Could use second line, resistance low
        datam_map_new[w2n,"second"] <- 0 # Can't use second line, resistance to d2 too high
      }
      
      ## remove nn and count
      if("nn" %in% colnames(datam_map_new)){
        w<-which(colnames(datam_map_new) == "nn"); w2<-which(colnames(datam_map_new) == "count"); w3<-which(colnames(datam_map_new) == "nh")
        datam_map_new <- datam_map_new[,-c(w,w2,w3)]}
    }
    
    if(!("cisol1" %in% colnames(datam_map_new))){
      datam_map_new$cisol1 <- NA
    } 
    
    if(!("cisol2" %in% colnames(datam_map_new))){
      datam_map_new$cisol2 <- NA
    }
    
    ######## where should use third line? ###########################
    # where should use third line, if available?
    # third line therapy for this synd = d3 ### NOT DONE FOR COMBINATION THERAPY
    d3 <- model()$d3 #as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),8])
    
    dd <- datam_map_new # just do to get cisol... 
    dd <- dd[which(dd$variable %in% c(d3)),]
    dd <- dd %>% dplyr::group_by(Country) %>% dplyr::summarise(nn = sum(res_perc < res_cut), 
                                                               count = n(), cisol3 = sum(n), nh = sum(res_perc > res_cut)) 
    if(dim(dd)[1]>0){
      datam_map_new <- merge(dd,datam_map_new, by = "Country",all.y=TRUE)}
    
    ### Where resistance high? 
    wr <- which(datam_map_new$res_perc >  res_cut) # resistance higher than res_cut input (20%)
    wd3 <- which(datam_map_new$variable == d3)
    w3n <- intersect(wd3, wr) # where information and resistance to d3 > res_cut
    w3y <- setdiff(wd3, wr) # where information and resistance to d3 < res_cut
    
    datam_map_new[w3y,"third"] <- 1 # Could use third line, resistance low
    datam_map_new[w3n,"third"] <- 0 # Can't use third line, resistance to d3 too high
    
    if(!("cisol3" %in% colnames(datam_map_new))){
      datam_map_new$cisol3 <- NA
    } 
    
    ### Need to group the above information by country
    suma = function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)
    data_treat <- ddply(datam_map_new, .(Country), summarise, 
                        first_c = suma(first),  
                        second_c = suma(second),  
                        third_c = suma(third),
                        nisol = ifelse(sum(is.na(cisol1),is.na(cisol2),is.na(cisol3))< 3, min(cisol1,cisol2,cisol3,na.rm = TRUE) , 11))# Only want to say if result and < 10. #min(cisol1, cisol2, cisol3, na.rm = TRUE))
    
    ### Need at least 10 isolates
    data_treat$"Number of isolates" <- ""
    data_treat[which(data_treat$nisol < 10),"Number of isolates"] <- "< 10"
    
    # Sum over all the information given above
    suma = function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)
    data_treat$sum_t <- apply(data_treat[,c("first_c","second_c","third_c")],1,suma)
    
    #### New column for reccomendation 
    data_treat$recc <- 0
    ## Look at sum_t
    # if NA then no data on any of the three
    w<-which(!is.na(data_treat$sum_t)) # Those that are not NA
    
    # if 1 then one of columns had data 
    w11 <- which(data_treat$first_c >= 1) # Here can use first line ... 
    data_treat[w11,"recc"] <- 1 # so use it! 
    
    # These have information but can't use first 
    w10 <- which(data_treat$first_c == 0) # Here can't use first line ... 
    data_treat[w10,"recc"] <- 2.1 # so recommend second, based on first
    
    # But some will have data on second line
    w21 <- intersect(w10, which(data_treat$second_c >= 1)) # Here can't use first line ... 
    data_treat[w21,"recc"] <- 2 # and can use second so use it
    
    w20 <- intersect(w10, which(data_treat$second_c == 0)) # Here can't use first or second line ... 
    data_treat[w20,"recc"] <- 3.12 # so recommend third
    
    # But some will have data on third line
    w31 <- intersect(w20, which(data_treat$third_c >= 1)) # Here can't use first or second line ... 
    data_treat[w31,"recc"] <- 3 # and can use third so use it
    
    w30 <- intersect(w20, which(data_treat$third_c == 0)) # Here can't use first, second or third line ... 
    data_treat[w30,"recc"] <- 4.123 # so recommend fourth
    
    ### FOR NOW - IF NO DATA ON FIRST THEN CAN'T SAY DON'T USE
    suma = function(x) if (all(is.na(x))) x[NA_integer_] else max(x, na.rm = TRUE)
    ## But label where have data on others 
    data_treat$nofirstd <- as.numeric(apply(data_treat[,c("first_c","second_c","third_c")],1,suma))
    w0<-which(data_treat$recc == 0)
    wd<-which(!is.na(data_treat$nofirstd))
    data_treat[intersect(w0,wd),"recc"] <- 0.6 # where currently recc = 0, but data on 2nd or 3rd... round up to 1
    
    if(any(is.na(data_treat[w,"recc"]))){print("STOP")} # if any not assigned
    
    # Remove those whose recc is zero - don't want to plot
    wrc<-which(data_treat$recc == 0)
    #print(data_treat[wrc,]) # those with no data on the drugs used
    if(length(wrc) > 0 ){data_treat <- data_treat[-wrc,]} # remove so they are grey in the map
    
    # return all object as a list
    list(data_treat = data_treat)
  }
  )
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$TEXTX - Warnings texts if no data
  ### ****************************************************************************************************************************************###
  
  ### Text if no prevalence map
  output$text1 <- renderText({
    if(model()$p1 ==0) {paste0("No data for first line drug in this database")
    }else{""}
  })
  
  output$text2 <- renderText({
    if(model()$p2 ==0){paste0("No data for second drug in first line therapy in this database")
    }else{""}
  })
  
  output$text3 <- renderText({
    if(model()$p3 ==0) {paste0("No data for any of the therapy options in this database")
    }else{""}
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$TEXT_WARNING_DATA - Warnings texts around syndromes included in each dataset
  ### ****************************************************************************************************************************************###
  
  output$text_warning_ecdc <- renderText({
    if("ECDC" %in% input$checkGroup){if(input$variable != "sep"){
      "No data for this syndrome in the ECDC database: only sepsis is included."}else{""}
    }else{""}
  })
  output$text_warning_resm <- renderText({
    if("RESIST" %in% input$checkGroup){if(input$variable != "sep"){
      "No data for this syndrome in the RESISTANCEMAP database: only sepsis is included."}else{""}
    }else{""}
  })
  output$text_warning_glas <- renderText({
    if("GLASS" %in% input$checkGroup){if(!input$variable %in% c("sep","c.uti","puc")){
      "No data for this syndrome in the GLASS database: only sepsis,
      complicated UTI and purulent urthritis / cervicitis are included."}else{""}
    }else{""}
  })
  
  output$text_warning_ecdc2 <- renderText({
    if("ECDC" %in% input$checkGroup){if(input$variable != "sep"){
      "No data for this syndrome in the ECDC database: only sepsis is included."}else{""}
    }else{""}
  })
  output$text_warning_resm2 <- renderText({
    if("RESIST" %in% input$checkGroup){if(input$variable != "sep"){
      "No data for this syndrome in the RESISTANCEMAP database: only sepsis is included."}else{""}
    }else{""}
  })
  output$text_warning_glas2 <- renderText({
    if("GLASS" %in% input$checkGroup){if(!input$variable %in% c("sep","c.uti","puc")){
      "No data for this syndrome in the GLASS database: only sepsis,
      complicated UTI and purulent urthritis / cervicitis are included."}else{""}
    }else{""}
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$SLIDERS - generates sliders for contributing pathogen distribution 
  ### ****************************************************************************************************************************************###
  
  output$sliders <- renderUI({
    
    synd <- input$variable
    pvars <- as.character(sp_all[,1]) #as.character(unique(datam_map$Species)) # all species in data used 
    
    sp_b <- melt(sp_bkdwn[which(sp_bkdwn[,1] == synd),],id.vars = "syndrome")
    sp_b$variable <- gsub(".", ' ', sp_b$variable, fixed = T)
    ## Need "Streptococcus, viridans group"
    g <- grep("viri",sp_b$variable)
    sp_b[g,"variable"] <- "Streptococcus, viridans group"
    
    values_slid <- filter(sp_b, variable %in% pvars)
    #values_in_syn <- c(100*unique(datam_map2$Prop_Syn), rep(0,length(pvars_diff)))
    
    lapply(seq(values_slid$variable), function(i) {
      sliderInput(inputId = paste0("range",i),
                  label = em(values_slid$variable[i]),
                  min = 0, max = 100, value = values_slid$value[i])
    })
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$SLIDERS_TOTAL - generates total of inputted sliders for contributing pathogen distribution to check sums to 100
  ### ****************************************************************************************************************************************###
  
  output$sliders_total <- renderText({
    
    synd <- input$variable
    pvars <- as.character(sp_all[,1]) #as.character(unique(datam_map$Species)) # all species in data used 
    
    sp_b <- melt(sp_bkdwn[which(sp_bkdwn[,1] == synd),],id.vars = "syndrome")
    sp_b$variable <- gsub(".", ' ', sp_b$variable, fixed = T)
    values_slid <- filter(sp_b, variable %in% pvars)
    
    pvars_use <- pvars
    names <- c()
    for(i in 1:length(pvars_use)){names <- c(names,paste0("range", i))}
    
    values <- c()
    for(i in 1:length(pvars_use)){values <- c(values,input[[names[i]]])}
    
    if(sum(values)<101){paste0("Current total is: ", sum(values))
    }else{paste0("Current total is: ", sum(values), ".\n This is not 100. PLEASE ADJUST")
    }
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$INPUTPLOT - creates barchart for left panel of inputted (or baseline) contributing bacterial species
  ### ****************************************************************************************************************************************###
  
  output$inputplot <- renderPlot({
    
    cols <- colorRampPalette(brewer.pal(6,"Reds"), bias = 2)(13)
    
    synd <- input$variable
    pvars <- as.character(sp_all[,1]) #as.character(unique(datam_map$Species)) # all species in data used 
    
    sp_b <- melt(sp_bkdwn[which(sp_bkdwn[,1] == synd),],id.vars = "syndrome")
    sp_b$variable <- gsub(".", ' ', sp_b$variable, fixed = T)
    ## Need "Streptococcus, viridans group"
    g <- grep("viri",sp_b$variable)
    sp_b[g,"variable"] <- "Streptococcus, viridans group"
    
    values_slid <- filter(sp_b, variable %in% pvars)
    
    pvars_use <- pvars
    names <- c()
    for(i in 1:length(pvars_use)){names <- c(names,paste0("range", i))}
    
    values <- c()
    for(i in 1:length(pvars_use)){values <- c(values,input[[names[i]]])}
    
    bar_plot_data <- as.data.frame(cbind(rep(1, length(values)),values))
    colnames(bar_plot_data) <- c("x","value")
    bar_plot_data$Bacteria <- values_slid$variable
    
    p = ifelse(sum(bar_plot_data$value)!=100, 
               'atop(bold("Proportions do not equal 100,\nplease adjust"))', "")
    
    ggplot(bar_plot_data,aes(x = x, y = value, fill = Bacteria) ) + 
      scale_y_continuous(lim = c(0,100),"Percentage of syndrome (%)")+ 
      geom_bar(stat="identity") +  coord_flip() +
      guides(fill=guide_legend(nrow=length(pvars),byrow=TRUE)) +
      guides(fill=guide_legend(ncol=2)) + 
      theme(legend.text = element_text(face = rep("italic", length(pvars_use)))) +
      theme(legend.position = "top") + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + 
      annotate("text", x = 0.5, y = 50, size = 5, label = p, colour = "black", parse = TRUE)
    
  })
  
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$COOLPLOTX - Functions for creating MAP OF PREVALENCE of resistance to first line drugs in first line therapy
  ### ****************************************************************************************************************************************###
  
  output$coolplot <- renderPlot({
    
    # Data
    mapped_data <- model()$mapped_data_d11
    if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){mapped_data <- model()$europe_map_11}
    
    if(model()$p1 == 1){
      cols <- rev(colorRampPalette(brewer.pal(11,"Spectral"), bias = 2)(13))
      
      ### Map resistance proportion in this syndrome
      mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "res_perc", 
                                  catMethod = seq(0,100,5),
                                  colourPalette = cols, 
                                  addLegend = FALSE, 
                                  mapTitle = paste0("Percentage of syndrome isolates that are
                                \nresistant to first line empirical therapy: ",model()$d11),
                                  missingCountryCol = gray(.8),
                                  nameColumnToHatch = "n10")
      do.call( addMapLegend, c( mapParams
                                , legendLabels="all"
                                , legendWidth=0.5))
    }
  })
  
  ### MAP OF PREVALENCE of resistance to second line drug in first line therapy
  output$coolplot2 <- renderPlot({
    
    if(model()$d12 != ""){
      # Data
      mapped_data <- model()$mapped_data_d12
      if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){mapped_data <- model()$europe_map_12}
      
      if(model()$p2 == 1){
        cols <- rev(colorRampPalette(brewer.pal(11,"Spectral"), bias = 2)(13))
        
        titlet <- paste0("Percentage of syndrome isolates that are
                      \nresistant to the second drug in the first line empirical therapy: ",model()$d12)
        
        ### Map resistance proportion in this syndrome
        mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "res_perc", 
                                    catMethod = seq(0,100,5),
                                    colourPalette = cols, 
                                    addLegend = FALSE, 
                                    mapTitle = titlet,
                                    missingCountryCol = gray(.8),
                                    nameColumnToHatch = "n10")
        do.call( addMapLegend, c( mapParams
                                  , legendLabels="all"
                                  , legendWidth=0.5))
        
      } 
    }
    
  })
  
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$TRAFFICPLOT - Functions for creating MAP OF WHERE RESISTANCE TO EACH THERAPY above the cutoff 
  ### ****************************************************************************************************************************************###
  
  output$trafficplot <- renderPlot({
    
    datam_map_new <- treat()$data_treat
    
    if(model()$p3 == 1 & dim(datam_map_new)[1]>0){
      
      cols <- c("green","yellow","orange","red")
      
      ### Map resistance proportion in this syndrome
      # round recc for map
      datam_map_new$recc <- round(datam_map_new$recc,0)
      datam_map_new <- as.data.frame(datam_map_new)
      
      mapped_data <- joinCountryData2Map(datam_map_new, joinCode = "NAME", nameJoinColumn = "Country")
      mapped_data <- subset(mapped_data, continent != "Antarctica")
      if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){mapped_data <- subset(mapped_data, NAME%in%ecdc_countries)}
      
      t1 <- model()$d11
      if(model()$d12 != ""){t1 <- paste0(model()$d11,"&",model()$d12)}
      t2 <- model()$d21
      if(model()$d22 != ""){t2 <- paste0(model()$d21,"&",model()$d22)}
      t3 <- model()$d3
      
      mapParams <- mapCountryData(mapped_data, nameColumnToPlot = "recc", 
                                  catMethod = seq(0,4,1),
                                  colourPalette = cols, 
                                  addLegend = 'FALSE', 
                                  mapTitle = ifelse(length(model()$d3)>0, 
                                                    paste0("Can first (", t1,", green), second (", t2,", yellow) or\n third line therapy (", t3,", orange) be used?\nOr is resistance to all recommended therapies seen (red)?"),
                                                    paste0("Can first (", t1,", green), second (", t2,", yellow)\nOr is resistance to all recommended therapies seen (red)?")),
                                  missingCountryCol = gray(.8))
      
      #mapParams$legendText <- c('first','second','third','none')
      
      do.call( addMapLegend, c(mapParams)) #, list(legendText=c('all','first','second','third','none'))))#legendLabels="all")) #, x='bottom',title="Region",horiz=TRUE))
    }
    
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$EMPIRICTHERAPY - Table in first input of our baseline therapy
  ### ****************************************************************************************************************************************###
  
  ## FUTURE area to update: allow users to change this and input their own drug choices (limit to certain set for each therapy? setting? patient group?)
  output$empirictherapy <- renderTable({
    
    synd = input$variable
    
    # If not ATLAS then have to use grouped drugs
    if(length(input$checkGroup) > 1 || input$checkGroup != "ATLAS"){
      drg_bkdwn <- drg_bkdwn_group
    }
    
    # What is the empiric therapy for this syndrome? 
    emp_data <- drg_bkdwn[which(drg_bkdwn$syndrome == synd),]
    
    emp_data <- emp_data[,colSums(is.na(emp_data))<nrow(emp_data)]
    
    emp_data[,-1] 
    
  })
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$ECON - Table in first input of economic availability and usage of drugs
  ### ****************************************************************************************************************************************###
  
  output$econ <- renderDataTable({
    
    # Data - as in empiric therapy table
    synd = input$variable
    
    # If not ATLAS then have to use grouped drugs
    if(length(input$checkGroup) > 1 || input$checkGroup != "ATLAS"){
      drg_bkdwn <- drg_bkdwn_group
    }
    
    emp_data <- drg_bkdwn[which(drg_bkdwn$syndrome == synd),]
    emp_data <- emp_data[,colSums(is.na(emp_data))<nrow(emp_data)]
    
    # Grab from econ data
    which_drugs <- unique(unlist(emp_data[,-c(1,2)]))
    # If not ATLAS then have to use grouped 
    if(length(input$checkGroup) > 1 || input$checkGroup != "ATLAS"){
      econ_drugs <- filter(econ, group %in% which_drugs)
      colnames(econ_drugs)[colnames(econ_drugs) == "group"] <- "Antibiotic"
    }else{ 
      econ_drugs <- filter(econ, Antibiotic.AR.IA %in% which_drugs)
      colnames(econ_drugs)[colnames(econ_drugs) == "Antibiotic.AR.IA"] <- "Antibiotic"
    }
    
    colnames(econ_drugs)[colnames(econ_drugs) == "WHO.Access"] <- "WHO EML"
    colnames(econ_drugs)[colnames(econ_drugs) == "AWaRe"] <- "WHO AWaRe"
    
    econ_drugs <- econ_drugs[,c("Antibiotic","Cost","Formulation","WHO EML", "WHO AWaRe")]
    ee <- datatable(econ_drugs) %>% formatStyle(
      'WHO EML',
      backgroundColor = styleEqual(c("Essential Medicine - Core List","Not on the list","Essential Medicine in different dosage/strength","Essential Medicine - Complementary List",
                                     "Potential alternative to an Essential Medicine"), c('red', 'green','orange','orange','green'))
    ) %>% formatStyle(
      'WHO AWaRe',
      backgroundColor = styleEqual(c("Watch","Reserve","Access"), c('orange','red','green'))
    )
    
    ee
  })
  
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$RESULTS - Table of underlying data
  ### ****************************************************************************************************************************************###
  
  output$results <- renderDataTable({
    
    synd = input$variable
    datam_map <- model()$datam_map
    datam_map <- datam_map[which(datam_map$syndrome == synd),]
    
    if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){
      datam_map <- subset(datam_map, Country%in%ecdc_countries)}
    
    age <- input$age
    
    # If not ATLAS then have to use grouped drugs
    if(length(input$checkGroup) > 1 || input$checkGroup != "ATLAS"){
      drg_bkdwn <- drg_bkdwn_group
    }
    
    # first line drug for this synd
    # Drugs for this synd
    d11 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"First.line.1st.drug"])
    d12 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"First.line.2nd.drug"])
    w1 <- c(which(datam_map$variable == d11),which(datam_map$variable == d12))
    
    d21 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Second.line.1st.drug"])
    d22 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Second.line.2nd.drug"])
    w2 <- c(which(datam_map$variable == d21),which(datam_map$variable == d22))
    
    d3 <- as.character(drg_bkdwn[intersect(which(drg_bkdwn$syndrome == synd),which(drg_bkdwn$Age == age)),"Third.line.1st.drug"])
    w3 <- which(datam_map$variable == d3)
    
    datam_map$Therapy <- ""
    datam_map[w1,"Therapy"] <- "First"
    datam_map[w2,"Therapy"] <- "Second"
    datam_map[w3,"Therapy"] <- "Third"
    if(length(intersect(w1,w2))>1){datam_map[intersect(w1,w2),"Therapy"] <- "First/Second"}
    if(length(intersect(w2,w3))>1){datam_map[intersect(w2,w3),"Therapy"] <- "Second/Third"}
    if(length(intersect(w1,w3))>1){datam_map[intersect(w1,w3),"Therapy"] <- "First/Third"}
    
    datam_map <- datam_map[c(w1,w2,w3),]
    
    datam_map$Antibiotic <- datam_map$variable
    datam_map$"Resistance prevalence (percentage)" <- round(100*datam_map$rate_r,2)
    dd <- datam_map[,c("Country","Species","Therapy","Antibiotic","n","Resistance prevalence (percentage)")]
    dd <- dd[order(dd$Country),]
    dd[,c("Country","Species","Therapy","Antibiotic","n","Resistance prevalence (percentage)")]
    
  }
  )
  
  ### ****************************************************************************************************************************************###
  # OUTPUT$RECC - Recommendations table 
  ### ****************************************************************************************************************************************###
  
  output$recc <- renderTable({
    
    if(model()$p3 == 1 & dim(treat()$data_treat)[1]>0){
      # Data
      drec <- ""
      drec <- treat()$data_treat
      
      if(length(input$checkGroup) == 1 && input$checkGroup == "ECDC"){
        drec <- subset(drec, Country%in%ecdc_countries)}
      
      # Relabel columns
      t1 <- paste0("First line: ",model()$d11) 
      if(model()$d12 !=""){t1 <- paste0("First line: ",model()$d11," & ",model()$d12)} 
      colnames(drec)[colnames(drec) == "first_c"] <- t1
      t2 <- paste0("Second line: ",model()$d21)
      if(model()$d22 !=""){t2 <- paste0("Second line: ",model()$d21," & ",model()$d22)}
      colnames(drec)[colnames(drec) == "second_c"] <- t2
      colnames(drec)[colnames(drec) == "third_c"] <- paste0("Third line: ",model()$d3)
      
      # Recode rows
      w1<-which(drec$recc == 1)
      w2<-which(drec$recc == 2)
      w3<-which(drec$recc == 3)
      w4<-which(drec$recc == 2.1)
      w5<-which(drec$recc == 3.12)
      w6<-which(drec$recc == 0.6)
      w7<-which(drec$recc == 4.123)
      
      drec$Recommendation = ""
      drec[w1,"Recommendation"] = "Use first line"
      drec[w2,"Recommendation"] = "Use second line"
      drec[w3,"Recommendation"] = "Use third line (if exists)"
      drec[w4,"Recommendation"] = "Use second line, as resistance to first (but no data on resistance to second)"
      drec[w5,"Recommendation"] = "Use third line, as resistance to both first and second (but no data on resistance to third)"
      drec[w6,"Recommendation"] = "Use first line, but no data to inform - consider second or third if data"
      drec[w7,"Recommendation"] = "Consider alternatives! Resistance to all recommended therapies seen"
      
      #### Add in proportion of syndrome allocated
      # drugs for this syndrome
      drug_list_this_synd <- c(model()$d11,model()$d12,model()$d21,model()$d22,model()$d3)
      
      w<-which(drug_list_this_synd == "")
      if(length(w) > 0){drug_list_this_synd <- drug_list_this_synd[-w]} 
      
      # Averages across the data for all checkgroups
      data2 <- datam_map_all[datam_map_all$Dataset %in% input$checkGroup,] # those in this group
      datam_map <- data2 %>%dplyr::group_by(Country,syndrome,Species,variable) %>%
        dplyr::summarise(rate_r=mean(rate_r), Prop_Syn = mean(Prop_Syn),n = sum(n)) # mean over all datasets
      
      # sums up coverage by syndrome and drug for each country. Only those drugs in this syndrome's treatment
      dd <- datam_map %>% group_by(Country, syndrome, variable) %>% dplyr::summarise(total_prop = sum(Prop_Syn)) %>%
        subset(variable %in% drug_list_this_synd) 
      
      dd2 <- subset(dd,syndrome == input$variable) %>% group_by(Country) %>% dplyr::summarise(mm = mean(total_prop)) # Takes mean across the variables in this syndrome
      drec2 <- merge(drec, dd2, by = "Country")
      colnames(drec2)[colnames(drec2) == "mm"] <- "Mean syndrome coverage"
      
      drec2[,c("Country",t1,t2,
               paste0("Third line: ",model()$d3),"Recommendation","Number of isolates","Mean syndrome coverage")]

      
    }
  },
  caption = "For the different therapy options (First, Second or Third line), the values in the table 
             signify if data is available for all drugs in the therapy (NA if not) and that the therapy 
             can be used (1) or not (0). The number of isolates is blank unless there were fewer than 
              10 isolates to inform the decision, where a '< 10' is indicated. The final column gives the mean syndrome coverage
            across the relevant drugs for this syndrome for this country. A value of 0.5 means that only 50% of the contributing pathogens
            were available from sources linked to this syndrome and had been tested for the relevant antibiotics.",
  caption.placement = getOption("xtable.caption.placement", "top"), 
  caption.width = getOption("xtable.caption.width", NULL)
  )
  
  
}

### ****************************************************************************************************************************************###
# RUN SHINY APP
### ****************************************************************************************************************************************###

shinyApp(ui = ui, server = server)