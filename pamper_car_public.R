############################## Description ##############################

## This is a commented file containing code to perform secondary analyses of PAMPer trial data (Sperry et al., 2018, NEJM)
## 
## CAR: Clustering and Regression
##
## Public code associated with Gruen et al., 2020 JCI Insight
##
## Last updated: March 2020

############################## Workspace Setup ##############################
## Set Up Workspace
    rm(list=ls()) # Clear all variables; delete all the objects in the workspace and start with a clean slate
    setwd("/Users/dgruen/Documents/Academia/Pitt_T32/projects/pamper/data/") # set the working directory
    set.seed(386) # make code reproducible

## Dependencies (these libraries must be loaded for code and analyses)
    library(haven)  # reads in SPSS files
    library(gdata)  # reads excel files and more
    library(RColorBrewer) # colors for graphs
    library(scales) # pretty plots
    library(reshape2) # used for plotting with ggplot2
    library(survminer) # creates survival curves
    library(survival) # required for survminer
    library(caTools) # classification & splitting data for validation etc. 
    library(broom) # way of dealing with model output results in a tidy way, need for GLM
    library(clustertend) # for validating clusters
    library(dendextend) # dendrogram building and visualization
    library(fpc) # Flexible Procedures for Clustering
    library(clValid) # Statistical and biological validation of clustering results
    library(pvclust) # hierarchical clustering with p-values
    library(mclust) # model-based clustering, classification, and density estimation
    library(NbClust) # determining the number of clusters
    library(RCurl) # cross validation
    library(prettyR) # cross validation
    library(lsmeans) # Least-squares means
    library(multcompView) # Summarize multiple paired comparisons
    library(rcompanion) # Variety of Functions to Support Extension Education Program Evaluation
    library(nlme) # Anova
    library(corrplot) # visualize results of PCA
    library(robustbase) # for dendrograms for hierarchiacal clustering (*choose "n")
    library(scatterplot3d) # 3D scatterplots}
    library(factoextra) # pca & hierarchical clustering
    library(FactoMineR) # pca & hierarchical clustering
    library(ggfortify) # pca & hierarchical clustering
    library(devtools) # developer tools
    library(tableone) # to make a summmary table
    library(np) # nonparametric linear regression
    library(xtable) # make things go to latex e.g. tables
    library(jtools) # good for stats summaries
    library(ggstance) # required for jtools / summ
    library(sandwich) # robust standard errors for linear regression
    library(lmtest) # for lin reg
    library(multiwayvcov)  # for lin reg
    library(wesanderson) # Wes Anderson color palette for ggplot :)
    library(gplots) # for making nice heatmaps / heatmap.2
    library(dotwhisker) # dot and whisker plot for regression coefficients
    library(tidyverse)  # for pretty visuals (also includes things like haven and ggplot2)
    library(ggdendro) # for building heatmaps plus dendros 
    library(grid)

############################## File Setup ##############################

# Redacted

############################## Clean Data ##############################

# Redacted

############################## Manual Data Inspection ############################## 

# Redacted

############################## Basic Stats, Tables, Figures ##############################

## Table: Comparison of sampled and unsampled cohorts
    df = df_master_filtered_wide
    df$sampled = NA  # designate a new column
    
    df$sampled[which(is.na(df$IL_6_0))] = 0 # not sampled
    df$sampled[which(!is.na(df$IL_6_0))] = 1 # sampled
    df$sampled = as.factor(df$sampled) # make a factor
    
    myVars <- c("initial_GCS", "iss", "ais_head", "age", "gender", "alive_at_30", "TBI", "PH_intubation", "PH_sbp_70", "any_blunt", "PH_CPR", "PH_crystalloid", "PH_time", "PH_prbc", "PH_blood", "initial_GCS_8")
    catVars <- c("gender", "alive_at_30", "TBI", "PH_intubation", "PH_sbp_70", "any_blunt", "PH_CPR", "initial_GCS_8", "PH_blood")
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "sampled", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
    xtable(tab1)



## Table: Sampled biomarker cohort plasma vs. no plasma   
    df <- df_unique_pamper

    # df = filter(df, alive_at_30==1|alive_at_30==2)
    df = filter(df, !is.na(alive_at_30))
    
    # Convert head_ais from factor to number ! WARNING !
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    myVars <- c("age", "gender", "race",
                
                "initial_GCS", "initial_GCS_8", "total_gcs_score", "iss",
                "PH_sbp_70", "TBI", "chest_injury", "abdominal_injury", 
                "extremity_injury", "spinal_cord_injury", "any_blunt", "blunt_mechanism", 
                "any_penetrating", "penetrating_mechanism",
                "severe_head",  "ais_head", "vitals_hr", "vitals_sbp",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "transfer",
                
                "INR", "alpha", "K", "MA", "LY30", "hyperfibrinolysis",
                "transfusion_24h", "prbc_24h", "plasma_24h",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                
                "alive_at_30", "mortality_24h", "ed_coagulopathy", 
                "MOF", "ALI", "NI", "icu_los", "hospital_los", "mech_vent_days")
    
    catVars <- c("gender", "race", "alive_at_30", "TBI", "severe_head", "ed_coagulopathy", "chest_injury",
                 "abdominal_injury", "extremity_injury", "spinal_cord_injury",
                 "ALI", "NI", "MOF", "hyperfibrinolysis", "PH_intubation", "PH_sbp_70",
                 "any_blunt", "blunt_mechanism",  "any_penetrating", "penetrating_mechanism",
                 "PH_CPR", "initial_GCS_8", "plasma_on_helicopter", "transfer",
                 "PH_blood", "ed_trali", "mortality_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h"
    )
    
    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "FFP", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, exact = "crani_24", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    # xtable(tab1)
    View(tab1)


## Table: Comparison of marker values for plasma vs. no plasma
    df = df_master_filtered_wide
    # dput(names(df)) # output names
    # df = filter(df, iss>=30)
    myVars <- c("adiponectin_0", "adiponectin_24", "celldeath_0", "celldeath_24", 
                "GM_CSF_0", "GM_CSF_24", "GM_CSF_72", "IL_10_0", "IL_10_24", 
                "IL_10_72", "IL_17A_0", "IL_17A_24", "IL_17A_72", "IL_17E_0", 
                "IL_17E_24", "IL_17E_72", "IL_1b_0", "IL_1b_24", "IL_1b_72", 
                "IL_2_0", "IL_2_24", "IL_2_72", "IL_21_0", "IL_21_24", "IL_21_72", 
                "IL_22_0", "IL_22_24", "IL_22_72", "IL_23_0", "IL_23_24", "IL_23_72", 
                "IL_27_0", "IL_27_24", "IL_27_72", "IL_33_0", "IL_33_24", "IL_33_72", 
                "IL_4_0", "IL_4_24", "IL_4_72", "IL_5_0", "IL_5_24", "IL_5_72", 
                "IL_6_0", "IL_6_24", "IL_6_72", "IL_7_0", "IL_7_24", "IL_7_72", 
                "IL_8_0", "IL_8_24", "IL_8_72", "IL_9_0", "IL_9_24", "IL_9_72", 
                "IP_10_0", "IP_10_24", "IP_10_72", "MCP_1_0", "MCP_1_24", "MCP_1_72", 
                "MIG_0", "MIG_24", "MIG_72", "S100A10_0", "S100A10_24", "suPAR_0", 
                "suPAR_24", "syn_0", "syn_24", "TM_0", "TM_24", "TNFa_0", "TNFa_24", 
                "TNFa_72", "VEGF_0", "VEGF_24")
    
    tab1 <- CreateTableOne(vars = myVars, strata = "FFP", data = df, testNonNormal = kruskal.test)
    tab1 <- print(tab1, nonnormal = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
    xtable(tab1)
    tab1


## Table: comparison of marker values for survivors vs. nonsurvivors at 30 days
    df <- df_master_filtered_wide
    df <- filter(df, alive_at_30==1 | alive_at_30==2)
    df$alive_at_30 <- droplevels(df$alive_at_30)
    
    myVars <- c("adiponectin_0", "adiponectin_24", "celldeath_0", "celldeath_24", 
                "GM_CSF_0", "GM_CSF_24", "GM_CSF_72", "IL_10_0", "IL_10_24", 
                "IL_10_72", "IL_17A_0", "IL_17A_24", "IL_17A_72", "IL_17E_0", 
                "IL_17E_24", "IL_17E_72", "IL_1b_0", "IL_1b_24", "IL_1b_72", 
                "IL_2_0", "IL_2_24", "IL_2_72", "IL_21_0", "IL_21_24", "IL_21_72", 
                "IL_22_0", "IL_22_24", "IL_22_72", "IL_23_0", "IL_23_24", "IL_23_72", 
                "IL_27_0", "IL_27_24", "IL_27_72", "IL_33_0", "IL_33_24", "IL_33_72", 
                "IL_4_0", "IL_4_24", "IL_4_72", "IL_5_0", "IL_5_24", "IL_5_72", 
                "IL_6_0", "IL_6_24", "IL_6_72", "IL_7_0", "IL_7_24", "IL_7_72", 
                "IL_8_0", "IL_8_24", "IL_8_72", "IL_9_0", "IL_9_24", "IL_9_72", 
                "IP_10_0", "IP_10_24", "IP_10_72", "MCP_1_0", "MCP_1_24", "MCP_1_72", 
                "MIG_0", "MIG_24", "MIG_72", "S100A10_0", "S100A10_24", "suPAR_0", 
                "suPAR_24", "syn_0", "syn_24", "TM_0", "TM_24", "TNFa_0", "TNFa_24", 
                "TNFa_72", "VEGF_0", "VEGF_24")
    
    tab1 <- CreateTableOne(vars = myVars, strata = "alive_at_30", data = df, testNonNormal = kruskal.test)
    tab1
    tab1 <- print(tab1, nonnormal = FALSE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing = TRUE)
    xtable(tab1)
    # write.csv(tab1, file = "myTable_dsg.csv") #csv


## Table: comparison of marker values for survivors vs. nonsurvivors at 24 hours
    df = df_master_filtered_wide
    df = filter(df, !is.na(mortality_24h))
    myVars <- c("GM_CSF_0", "GM_CSF_24", "IL_10_0", "IL_10_24",
                "IL_17A_0", "IL_17A_24", "IL_17E_0", 
                "IL_17E_24", "IL_1b_0", "IL_1b_24", 
                "IL_2_0", "IL_2_24", "IL_21_0", "IL_21_24", 
                "IL_22_0", "IL_22_24", "IL_23_0", "IL_23_24", 
                "IL_27_0", "IL_27_24", "IL_33_0", "IL_33_24", 
                "IL_4_0", "IL_4_24", "IL_5_0", "IL_5_24", 
                "IL_6_0", "IL_6_24", "IL_7_0", "IL_7_24", 
                "IL_8_0", "IL_8_24", "IL_9_0", "IL_9_24", 
                "IP_10_0", "IP_10_24", "MCP_1_0", "MCP_1_24", 
                "MIG_0", "MIG_24", "syn_0", "syn_24", "TM_0", "TM_24", 
                "TNFa_0", "TNFa_24", "VEGF_0", "VEGF_24", "S100A10_0", "S100A10_24",
                "suPAR_0", "suPAR_24", "adiponectin_0", "adiponectin_24",
                "celldeath_0", "celldeath_24")
    
    tab1 <- CreateTableOne(vars = myVars, strata = "t_24h_censor", data = df, testNonNormal = kruskal.test)
    tab1 <- print(tab1, nonnormal = FALSE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
    tab1
    xtable(tab1)


## Figure: Box plot of all T0 markers plasma vs. no plasma
    df = df_master_filtered_long
    df = filter(df, !is.na(Biomarker))
    df = filter(df, hour==0)
    ggplot(df,
           aes(x=factor(plasma_on_helicopter),
               y=Concentration,
               group = factor(plasma_on_helicopter),
               color = factor(plasma_on_helicopter))) +
      geom_boxplot(outlier.size=0) +
      geom_jitter(width=0.3, alpha=0.15) +
      facet_wrap(~Biomarker, scales = "free")
    
    ## keep outliers but zoom in with facet wrap
    ## issue described here: https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
    ## function
    calc_boxplot_stat <- function(x) {
      coef <- 1.5
      n <- sum(!is.na(x))
      # calculate quantiles
      stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
      names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
      iqr <- diff(stats[c(2, 4)])
      # set whiskers
      outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
      if (any(outliers)) {
        stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
      }
      return(stats)
    }
    
    ## df and plot
    df <- df_master_filtered_long
    df <- filter(df, hour==0|hour==24)
    df <- filter(df, !is.na(Biomarker))
    df <- filter(df, !is.na(FFP))
    df <- filter(df, !is.na(hour))
    df <- filter(df, alive_at_30==1|alive_at_30==2)
    ggplot(df %>% 
             unite(twoby, hour, alive_at_30, remove = F), 
           aes(x = twoby,
               y = Concentration,
               fill=as.factor(alive_at_30))) +
      stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
      facet_wrap(~Biomarker, scales = "free") +
      theme_classic() +
      scale_fill_manual(values = c("#ffffff", "#969696")) +
      theme(aspect.ratio=1)


## Figure: T0 concentration vs. ISS for plasma vs. no plasma    
## Loess Curve
## Small version: only IL-6, syn, vegf, tm
    df = df_master_filtered_long
    df = filter(df, !is.na(Biomarker))
    df = filter(df, alive_at_30 == 1)
    df = filter(df, hour==0)
    df <- filter(df, Biomarker == "IL_6" |
                   #Biomarker == "IL_8" |
                   Biomarker == "MCP_1" |
                   #Biomarker == "TNF_a" |
                   #Biomarker == "celldeath" |
                   Biomarker == "syn" |
                   Biomarker == "VEGF")
    
    ggplot(df,
           aes(x = iss,
               y = Concentration,
               color=as.factor(plasma_on_helicopter),
               shape=as.factor(plasma_on_helicopter))) +
      scale_color_manual(values = c("#cb181d", "#737373")) +
      geom_point(alpha = 0.2, size = 2) +
      geom_smooth() +
      theme_minimal() +
      xlab("ISS") + 
      ylab("Concentration") +
      facet_wrap(~Biomarker, scales = "free")


## Compare markers: Survivors vs. nonsurvivors
    ## Compare markers: Survivors vs. nonsurvivors (24 hour mortality)
    df <- df_master_filtered_long
    df <- filter(df, !is.na(mortality_24h))
    
    # Do this in a loop for all biomarkers
    ubio <- unique(df$Biomarker) # vector of unique biomarkers
    resList <- list() # a list to store the results of each t test in
    for (i in 1:length(ubio)){
      grp1 <- filter(df, Biomarker == ubio[i]) %>% # subset df by relevant biomarker
        filter(mortality_24h == 0) 
      
      grp2 <- filter(df, Biomarker == ubio[i]) %>%
        filter(mortality_24h == 1) 
      
      resList[[i]] <- wilcox.test(grp1$Concentration, grp2$Concentration) # run test and store result
    }
    
    resList


############################## Regression Analysis ##############################
## 24  Hour Model
## Data Frame
    df <- df_master_filtered_long
    df <- filter(df, !is.na(Biomarker))
    df <- filter(df, hour == 24)
    df <- filter(df, iss > 30)

## Loop & results setup
    term <- unique(df$Biomarker) # vector of unique biomarkers
    res_df <- NULL;

## Model 
    for (i in 1:length(term)){
      
      df.i <- filter(df, Biomarker == term[i]) # subset df by relevant biomarker
      
      glm.results <- glm(Concentration ~ 
                           FFP + iss + INR + 
                           ais_head +
                           initial_GCS + PH_sbp_70  +
                           PH_crystalloid + PH_prbc
                         ,
                         data = df.i)
      
      ## Coeff. Results and store p values 
      coeffs <- tidy(coeftest(glm.results, vcov = sandwich)) # robust; sandwich
      FFP_coeff = coeffs[2,]
      res_df = rbind(res_df, FFP_coeff)
      
    }
    
    res_df <- cbind(term, res_df)
    res_df

## Table: Beta Coefficients for 24 hour model
    df <- res_df
    df <- df[,-c(2,5)] # specify a df and drop the second column
    # xtable(df, digits=c(0,2,2,2,2,-2))
    # xtable(df, digits=c(0,2,2,2,-2))
    xtable(df, digits=c(0,2,2,2,6))

## 0 Hour Model
## Data Frame
    df <- df_master_filtered_long
    df <- filter(df, !is.na(Biomarker))
    df <- filter(df, hour == 0)

## Loop & results setup
    term <- unique(df$Biomarker) # vector of unique biomarkers
    res_df <- NULL;

## Model 
    for (i in 1:length(term)){
      
      df.i <- filter(df, Biomarker == term[i]) # subset df by relevant biomarker
      
      glm.results <- glm(Concentration ~ 
                           FFP + iss + INR + 
                           ais_head +
                           initial_GCS + PH_sbp_70  +
                           PH_crystalloid + PH_prbc
                         ,
                         data = df.i)
      
      ## Coeff. Results and store p values 
      coeffs <- tidy(coeftest(glm.results, vcov = sandwich)) # robust; sandwich
      FFP_coeff = coeffs[2,]
      res_df = rbind(res_df, FFP_coeff)
      
    }
    
    res_df <- cbind(term, res_df)
    res_df

## Table: Beta Coefficients for 0 hour model
    df <- res_df
    df <- df[,-c(2,5)] # specify a df and drop the second column
    # xtable(df, digits=c(0,2,2,2,2,-2))
    # xtable(df, digits=c(0,2,2,2,-2))
    xtable(df, digits=c(0,2,2,2,6))

############################## PCA and HCA ##############################
## DF
    df <- biomarkers_for_PCA_select
    
    tally_pca <- df %>%
      tally()
    tally_pca
    
    df <- as.data.frame(df)
    df <- scale(df)

#PCA
    pca <- prcomp(t(df))
    
    ## Visualize
    summary(pca)

    loadings <- pca$rotation
    a <- loadings[,1]
    b <- loadings[,2]
    c <- loadings[,3]  
    d <- loadings[,4]  
    e <- loadings[,5] 
    f <- loadings[,6] 
    g <- loadings[,7] # keep 7 (75% variance)
    
    bound_loadings <- cbind(a,b,c,d,e,f,g)

## HCA
    df <- bound_loadings
    dist <- get_dist(df, method = "spearman")
    hc <- hclust(dist, method = 'ward.D2')
    plot(hc, axes=F, xlab='', ylab='', sub ='')
    rect.hclust(hc, k=2, border='red')
    
    ## get cluster#
    groups <- cutree(hc, k = 2)
    x<-cbind(biomarkers_for_PCA_select, groups)
    
    ## make a DF for coloring
    df <- df_master_filtered_wide_T0
    df$color_variable <- df$plasma_on_helicopter
    df <- select(df, "stnum", "color_variable")
    df <-  merge(df, cytokines_for_PCA_select_with_PAMPerID, by="stnum", all.y = TRUE)
    df <- select(df, "pamperID", "color_variable")
    df$color_variable <- as.integer(df$color_variable)
    whatyouwant <- setNames(df$color_variable, df$pamperID)

    ## Heatmap
    df <- biomarkers_for_PCA_select
    df <- scale(df)
    
    # Colors=rev(brewer.pal(10, "Spectral")) ## 11 and 9 also look good
    Colors= rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)) # also looks nice / blues & reds
    col1 <- brewer.pal(8, "Set1")
    display.brewer.pal(3, "BrBG")
    
    # make a color pal thats dark grey for plasma, light grey for SC and white for NA
    plas_color_pal <- c("#4d4d4d", "#e0e0e0", "#ffffff")
    display.brewer.pal(3, plas_color_pal)

    heatmap(as.matrix(df),
            Rowv=as.dendrogram(hc),
            col=Colors,
            xlab = "Biomarker", 
            ylab =  "PAMPer Patients",
            main = "HCPC Dendrogram and Biomarker Heatmap",
            RowSideColors=plas_color_pal[whatyouwant])

## Assess clusters
    ## Validate Clusters
    ## 30 indices for determining the optimal number of clusters  
    
    ## Raw Data
    df=biomarkers_for_PCA_select
    df=scale(df)
    nb <- NbClust(df, distance = "euclidean", min.nc = 2,
                  max.nc = 10, method = "ward.D2")
    ## 8 propose 2, 7 propose 4 splits. So at least 2 splits.
    ## According to the majority rule, the best number of clusters is  2

    ## PCA / precalculated distance matrix
    df=bound_loadings
    nb <- NbClust(df, diss=dist, distance = NULL, min.nc = 2,
                  max.nc = 10, method = "ward.D2")
    ## 6 propose 2, 4 propose 5 splits. So at least 2 splits.
    ## According to the majority rule, the best number of clusters is  2

## Assessment
    df = x
    df = rownames_to_column(df, var = "Pamp_num")
    
    df$pamperID = NA
    df$pamperID = substr(df$Pamp_num, 6, nchar(df$Pamp_num))
    df = merge (df, ID_key, by = "pamperID")
    
    df = select(df, -"Pamp_num") # delete this column, don't need
    df <- df %>%
      select(stnum, everything()) # moves stnum column to the beginning
    
    PAMPer_clustered = df
    
    df = PAMPer_clustered
    tally_clusters <- df %>%
      group_by(groups) %>%   # group by Biomarker
      tally()
    tally_clusters


## Figure: Overall survival (plasma vs standard care)
    df <- PAMPer_all_data_clustered
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df)
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "FFP",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               # conf.int = TRUE, # Add confidence interval
               title = "Survival Curves",
               pval = TRUE, # add p value
               pval.method = TRUE)    # Add p-value method name


## Figure: Survival vs. Clusters
    df <- PAMPer_all_data_clustered
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ cluster, data = df)
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Clusters",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               # conf.int = TRUE, # Add confidence interval
               title = "Survival Curves",
               pval = TRUE, # add p value
               pval.method = TRUE)    # Add p-value method name
    
## Figure: Survival Plasma vs. SOC  (facet: clusters 1&2)
    df<-PAMPer_all_data_clustered
    
    fit <- survfit(Surv(t_30d_mort_h, t_30d_censor_h) ~ plasma_on_helicopter, data = df)
    
    ggsurvplot(fit,
               legend = "bottom", 
               legend.title = "Prehospital Plasma",
               break.time.by = 250, # break time axis by 168 hours (make it weeks)
               legend.labs = c("Plasma", "No Plasma"),
               #conf.int = TRUE, # Add confidence interval
               title = "Survival Curves",
               pval = TRUE, # add pvalue
               pval.method = TRUE,
               facet.by = "cluster")    # Add p-value method names

## Assess group info
    df <- PAMPer_all_data_clustered
    # dput(df)
    
    tally_clusters <- df %>%
      group_by(cluster, FFP) %>%   # group by Biomarker
      tally()
    tally_clusters
    
    # Convert head_ais from factor to number
    df$ais_head <- as.numeric(levels(df$ais_head))[df$ais_head]
    
    df <- filter(df, alive_at_30==1 | alive_at_30==2) # note! must do this to calc p for alive@30
    
    myVars <- c("age", "gender",  "race",
                
                "initial_GCS", "initial_GCS_8", "total_gcs_score", "iss",
                "PH_sbp_70", "TBI", "chest_injury", "abdominal_injury", 
                "extremity_injury", "spinal_cord_injury", "any_blunt", "blunt_mechanism",
                "any_penetrating", "penetrating_mechanism",
                "severe_head", "ais_head", "vitals_hr", "vitals_sbp",
                
                "PH_intubation", "PH_CPR", "plasma_on_helicopter", "PH_crystalloid",
                "PH_prbc", "PH_blood", "PH_time", "PH_time_high", "transfer",
                
                "INR", "alpha", "K", "MA", "LY30", "hyperfibrinolysis",
                "transfusion_24h", "prbc_24h", "plasma_24h",
                "platelets_24h", "crystalloid_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                
                "alive_at_30", "mortality_24h", "ed_coagulopathy", 
                "MOF", "ALI", "NI", "icu_los", "hospital_los", "mech_vent_days"
                
                # "adiponectin_0", "adiponectin_24", "celldeath_0", 
                # "celldeath_24", "GM_CSF_0", "GM_CSF_24", 
                # "GM_CSF_72", "IL_10_0", "IL_10_24", "IL_10_72", "IL_17A_0", "IL_17A_24", 
                # "IL_17A_72", "IL_17E_0", "IL_17E_24", "IL_17E_72", "IL_1b_0", 
                # "IL_1b_24", "IL_1b_72", "IL_2_0", "IL_2_24", "IL_2_72", "IL_21_0", 
                # "IL_21_24", "IL_21_72", "IL_22_0", "IL_22_24", "IL_22_72", "IL_23_0", 
                # "IL_23_24", "IL_23_72", "IL_27_0", "IL_27_24", "IL_27_72", "IL_33_0", 
                # "IL_33_24", "IL_33_72", "IL_4_0", "IL_4_24", "IL_4_72", "IL_5_0", 
                # "IL_5_24", "IL_5_72", "IL_6_0", "IL_6_24", "IL_6_72", "IL_7_0", 
                # "IL_7_24", "IL_7_72", "IL_8_0", "IL_8_24", "IL_8_72", "IL_9_0", 
                # "IL_9_24", "IL_9_72", "IP_10_0", "IP_10_24", "IP_10_72", "MCP_1_0", 
                # "MCP_1_24", "MCP_1_72", "MIG_0", "MIG_24", "MIG_72", "S100A10_0", 
                # "S100A10_24", "suPAR_0", "suPAR_24", "syn_0", "syn_24", "TM_0", 
                # "TM_24", "TNFa_0", "TNFa_24", "TNFa_72", "VEGF_0", "VEGF_24"
    )
    
    catVars <- c("gender","race", "alive_at_30", "TBI", "severe_head", "ed_coagulopathy", "chest_injury",
                 "abdominal_injury", "extremity_injury", "spinal_cord_injury",
                 "ALI", "NI", "MOF", "hyperfibrinolysis", "PH_intubation", "PH_sbp_70",
                 "any_blunt", "blunt_mechanism", "any_penetrating", "penetrating_mechanism",
                 "PH_CPR", "initial_GCS_8", "plasma_on_helicopter", "transfer",
                 "PH_blood", "ed_trali", "mortality_24h", "vaso_24h", "prbc_10_24h", "prbc_4_24h",
                 "PH_time_high"
    )

    tab1 <- CreateTableOne(vars = myVars, data = df, strata = "groups", factorVars = catVars)
    tab1 <- print(tab1, nonnormal = TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, missing=TRUE, showAllLevels = TRUE)
    # xtable(tab1)
    View(tab1)
