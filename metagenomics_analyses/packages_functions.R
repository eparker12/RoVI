####################################
### Script to load RoVI packages ###
####################################

library(phyloseq)
library(ggplot2)
library(knitr)
library(RColorBrewer)
library(gridExtra)
library(kableExtra)
library(plyr)
library(magrittr)
library(dplyr)
library(tidyr)
library(cluster)
library(plotly)
library(DT)
library(pheatmap)
library(reshape2)
library(lme4)
library(matrixStats)
library(randomForest)
library(vegan)
library(corrplot)
library(data.table)
library(ggpubr)
library(labdsv)
library(crossval)
library(randomcoloR)
library(shiny)
library(binom)
library(ggsignif)
library(DescTools)
library(RVAideMemoire)
library(scales)
library(cowplot)
library(DECIPHER)
library(phangorn)
library(wesanderson)
library(inlmisc)
library(formattable)
library(ggExtra)
library(sjstats)
library(ALDEx2)
library(FSA)
library(Maaslin2)

theme_set(theme_bw())

# set country colours
India_col = "#CC79A7"  #"#3f9d7f" #"#66c2a5"
India_OPV_col = "#3f9d7f"
India_IPV_col = "#045a8d"
India_col_2 = "#66c2a5" #"#045a8d"
Malawi_col = "#E69F00"
UK_col = "#0072B2"
r_col = "#045a8d"
nr_col = "#d95f02"

########################################
### Functions for RoVI analyses work ###
########################################

# function to summarise phyloseq object at various stages of filtering
filtering_stats = function(input_ps) {
  nsamples = nsamples(input_ps)
  ntaxa = sum(taxa_sums(input_ps)>0)
  total_count = sum(otu_table(input_ps))
  min_count = min(sample_sums(input_ps))
  mean_count = round(mean(sample_sums(input_ps)))
  sd_count = round(sd(sample_sums(input_ps)))
  return(c(nsamples, ntaxa, total_count, min_count, mean_count, sd_count))
}

# function to extract n characters from right of string
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }

# function to plot prevalence vs mean abundance (when present) for sample subset
feature_abundance_function = function(input_rps, input_country="India", input_colour=India_col, yannotation=0.5) {
  
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps_t_no0 = input_rps
  otu_table(rps_t_no0)[otu_table(rps_t_no0)==0] = NA
  
  # calculate feature means
  tax_abund = data.frame(round(rowMeans(otu_table(rps_t_no0), na.rm=T),6), round(matrixStats::rowSds(otu_table(rps_t_no0), na.rm=T),6), rowSums(otu_table(input_rps)>0))
  names(tax_abund) = c("mean_abundance", "SD_abundance", "present")
  tax_abund$prevalence = round(tax_abund$present/ncol(otu_table(input_rps)),3)
  
  # append feature names
  tax_abund$feature = data.frame(as(tax_table(input_rps),"matrix"))$Feature
  
  # order by decreasing abundance, add rank
  tax_abund = tax_abund[order(-tax_abund$mean_abundance),]
  tax_abund$rank = 1:nrow(tax_abund)
  N_over_20_prev = sum(tax_abund$prevalence>=0.20)
  N_total = nrow(tax_abund)
  
  # Calculate cumulative abundance distribution
  tax_abund$cumulative = NA
  tax_abund$cumulative[1] = tax_abund$mean_abundance[1]
  for (i in 2:nrow(tax_abund)) { tax_abund$cumulative[i] = tax_abund$cumulative[i-1] + tax_abund$mean_abundance[i] }
  tax_abund$label = tax_abund$feature
  tax_abund$label[tax_abund$prevalence<0.5] = NA
  ggplot(tax_abund, aes(x=prevalence, y=mean_abundance, label=label)) +
    geom_point(alpha=0.5, size=2, colour=input_colour) +
    geom_text(size = 3, colour=input_colour, check_overlap = TRUE, nudge_y=0.01) +
    ylab("Mean relative abundance\nwhen present") + xlab("Prevalence") +
    scale_x_continuous(limits = c(0,1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    annotate("text", x = 0.55, y = yannotation, label = paste0(N_over_20_prev, "/", N_total, " with prevalence ≥20%"), colour=input_colour,  fontface = 'italic') +
    theme(axis.text = element_text(size=13), axis.title = element_text(size=13))
}

# longitudinal plot of alpha diversity
alpha_longitudinal = function(t1, outcome, outcome_title) {
  # select alpha diversity metric and outcome measure
  t1$outcome = t1[,outcome]
  
  # create plot
  ggplot(t1, aes(x = factor(country), y = Observed, colour=factor(outcome))) +
    geom_point(size = 1, alpha = 0.15, position = position_jitterdodge(jitter.width=0.35)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + scale_color_manual(outcome_title, values = c(nr_col, r_col)) +
    ggtitle("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Richness") + xlab(" ") +
    theme_bw() + theme(legend.position = "none") +
    theme(strip.background = element_blank()) +
    theme(axis.text = element_text(size=14), strip.text=element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=16, hjust=0.5))
}

# function to collate beta diversity results
collated_beta_function = function(rps, outcome, unweighted_logical=TRUE) {
  collated_beta = data.frame(t(data.frame(
    IND_BS3 = beta_outcome(rps, "India", outcome, unweighted=unweighted_logical),
    INDexp_BS3 = beta_outcome(rps, "India (exposed)", outcome, unweighted=unweighted_logical),
    INDunexp_BS3 = beta_outcome(rps, "India (unexposed)", outcome, unweighted=unweighted_logical),
    MLW_BS3 = beta_outcome(rps, "Malawi", outcome, unweighted=unweighted_logical)
  )))
  names(collated_beta) = c("R2", "p", "n")
  collated_beta$country = c("IND", "IND (neo+)", "IND (neo-)", "MLW")
  collated_beta$country =  factor(collated_beta$country, levels =  rev(c("IND", "IND (neo+)", "IND (neo-)", "MLW")))
  collated_beta$full = c(1,0,0,1)
  collated_beta
}

# function to calculate Bray-Curtis R-sq and p value for specified outcome and sample group
beta_outcome = function(input_rps, country_subset, outcome_selected, unweighted=TRUE) {
  if (outcome_selected == "BB2_IgA") {
    sample_data(input_rps)$outcome = log(as.numeric(as(sample_data(input_rps), "matrix")[,outcome_selected]))
  } else {
    sample_data(input_rps)$outcome = as(sample_data(input_rps), "matrix")[,outcome_selected]
  }
  
  input_m = data.frame(sample_data(input_rps))
  
  if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") {
    samplelist = subset(input_m, country_exposure==country_subset & !is.na(outcome))$sample_ID
  } else {
    samplelist = subset(input_m, country==country_subset & !is.na(outcome))$sample_ID
  }
  
  ps_t = prune_samples(samplelist, input_rps)
  t1 = data.frame(sample_data(ps_t))
  
  if (unweighted == TRUE) {
    beta_dist = phyloseq::distance(ps_t, method = "bray", binary=TRUE)
  } else {
    beta_dist = phyloseq::distance(ps_t, method = "bray")
  }
  
  adon = adonis(unname(beta_dist) ~ outcome, data = t1, permutations = 999)$aov.tab
  c(adon$R2[1]*100, adon$`Pr(>F)`[1], nrow(t1))
}

### function to run collection of Random Forests
collate_rf_binary = function(output_name) {
  
  # create dataframe for collated accuracy results
  collated_rf = data.frame(country = rep(NA,length(country_list)))
  collated_rf$responders = collated_rf$nonresponders = NA
  
  # run RF in loop across country strata
  for (c in 1:length(country_list)) {
    input_m$outcome = input_m[,outcome]
    country_subset = country_list[c]
    
    if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
      samplelist = subset(input_m, country_exposure==country_subset & !is.na(outcome))$sample_ID
    } else {
      samplelist = subset(input_m, country==country_subset & !is.na(outcome))$sample_ID 
    }
    input_rps = prune_samples(samplelist, rps)
    sample_data(input_rps)$outcome = as(sample_data(input_rps), "matrix")[,outcome]
    rf = run_rf_binary_outcome(input_rps, country_subset, prev_threshold=0.05, folds=5, iterations=20, n_per_group=50)
    
    collated_rf$country[c] = country_subset
    collated_rf$nonresponders[c] = as.character("nonresponders")
    collated_rf$responders[c] = as.character("responders")
    collated_rf$crossval_median[c] = rf$gini.cv$crossval_median[1]
    collated_rf$crossval_q1[c] = rf$gini.cv$crossval_quartile_1[1]
    collated_rf$crossval_q3[c] = rf$gini.cv$crossval_quartile_3[1]
  }
  write.csv(collated_rf,file=paste0("RF_outputs/",output_name))
}

### function to run Random Forests
run_rf_binary_outcome = function(input_rps, country_subset, prev_threshold=0.05, folds=5, iterations=iter, n_per_group=50) {
  
  # determine taxa present above threshold in at least one response subset
  ps1list = sample_names(input_rps)[sample_data(input_rps)$outcome==0] 
  ps1 = prune_samples(ps1list, input_rps)
  ps1 = filter_taxa(ps1, function(x) { sum(x > 0) > prev_threshold*nsamples(ps1) }, TRUE)
  ps2list = sample_names(input_rps)[sample_data(input_rps)$outcome==1]
  ps2 = prune_samples(ps2list, input_rps)
  ps2 = filter_taxa(ps2, function(x) { sum(x > 0) > prev_threshold*nsamples(ps2) }, TRUE)
  taxlist = unique(c(rownames(tax_table(ps1)),rownames(tax_table(ps2))))
  
  # if <50 responders or <50 non-responders, amend n_per_group to match minimum group size
  mingroup = min(c(length(ps1list), length(ps2list)))
  if (mingroup < 50) { n_per_group = mingroup }
  
  # create RF input
  rf_reference_rps = prune_taxa(taxlist, input_rps)
  rf_input = data.frame(t(otu_table(rf_reference_rps))) # Create matrix of OTU abundances 
  
  # create response-specific subsets
  rf_input$outcome = factor(sample_data(rf_reference_rps)$outcome)
  rf_input_1 = subset(rf_input, outcome==0)
  rf_input_2 = subset(rf_input, outcome==1)
  
  # perform crossval on whole dataset
  crossval_stats = list()
  crossval_importance = list()
  
  for (i in 1:iter) {
    # random draw for each iteration
    rf_s = rbind(rf_input_1[sample(nrow(rf_input_1),n_per_group),],rf_input_2[sample(nrow(rf_input_2),n_per_group),])
    # perform crossval for iteration above
    cv.rf = crossval_update(predfun, X=rf_s[,!grepl("outcome", names(rf_s))], Y=factor(rf_s$outcome), K=folds, B=1, verbose=FALSE)
    # write outputs
    crossval_stats[[i]] = cv.rf$stat.cv
    crossval_importance[[i]] = cv.rf$Gini_importance_output
  }
  
  # collate cross-validation statistics - stats
  stat.cv = as.data.frame(crossval_stats[[1]])
  for (i in 2:iter) { stat.cv = rbind(stat.cv, as.data.frame(crossval_stats[[i]])) }
  stat.cv$accuracy = (stat.cv$TP+stat.cv$TN)/(stat.cv$TP+stat.cv$TN+stat.cv$FP+stat.cv$FN)
  
  # collate crossvalidation statistics - gini
  gini.cv = as.data.frame(crossval_importance)
  gini.cv$mean = rowMeans(gini.cv[,1:(iter*folds)])
  gini.cv$sd = rowSds(as.matrix(gini.cv[,1:(iter*folds)]))
  gini.cv = gini.cv[,(ncol(gini.cv)-1):ncol(gini.cv)]
  gini.cv$feature_id = rownames(gini.cv)
  
  # append clean feature names
  features = data.frame(as(tax_table(rf_reference_rps), "matrix"))
  features$feature_id = rownames(features)
  gini.cv = merge(gini.cv, features, by = "feature_id")
  gini.cv <- gini.cv[order(-gini.cv$mean),]
  rownames(gini.cv) = gini.cv$feature_id
  
  # drop mode_delivery for Malawi
  if (country_subset=="Malawi") { covariate_subset = maaslin_covariates[maaslin_covariates!="mode_delivery"] } else { covariate_subset = maaslin_covariates }
  
  # append prevalence, abundance data and Fisher's test p values for presence/absence
  for (i in 1:nrow(gini.cv)) {
    feature = gini.cv$feature_id[i]
    gini.cv$nonresponder_present[i] = sum(rf_input_1[,feature]>0)
    gini.cv$nonresponder_n[i] = nrow(rf_input_1)
    gini.cv$nonresponder_prev[i] = round(sum(rf_input_1[,feature]>0)/nrow(rf_input_1)*100,1)
    gini.cv$nonresponder_mean[i] = mean(rf_input_1[,feature])*100
    gini.cv$nonresponder_sd[i] = sd(rf_input_1[,feature])*100
    
    gini.cv$responder_present[i] = sum(rf_input_2[,feature]>0)
    gini.cv$responder_n[i] = nrow(rf_input_2)
    gini.cv$responder_prev[i] = round(sum(rf_input_2[,feature]>0)/nrow(rf_input_2)*100,1)
    gini.cv$responder_mean[i] = mean(rf_input_2[,feature])*100
    gini.cv$responder_sd[i] = sd(rf_input_2[,feature])*100
    if (gini.cv$nonresponder_prev[i]>95 & gini.cv$responder_prev[i]>95) {
      gini.cv$fisher_p[i] = NA
      gini.cv$lr_unadjusted_p[i] = NA
      gini.cv$lr_adjusted_p[i] = NA
    } else {
      gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,feature]>0), sum(rf_input_1[,feature]==0),
                                                 sum(rf_input_2[,feature]>0), sum(rf_input_2[,feature]==0)),nrow=2))$p.value
      rf_input$feature_present = as.numeric(rf_input[,feature]>0)
      
      # create glm dataframe
      rf_input$age_at_first_dose = as.numeric(sample_data(rf_reference_rps)$age_at_first_dose)
      rf_input$antibiotic_exposure = factor(sample_data(rf_reference_rps)$antibiotic_exposure)
      rf_input$birth_weight = as.numeric(sample_data(rf_reference_rps)$birth_weight)
      rf_input$breastfed_child = factor(sample_data(rf_reference_rps)$breastfed_child)
      rf_input$mode_delivery = factor(sample_data(rf_reference_rps)$mode_delivery)
      
      # calculate outputs
      gini.cv$lr_unadjusted_p[i] = summary(glm(feature_present ~ outcome, data = rf_input, family = 'binomial'))$coefficients[2,4]
      gini.cv$lr_adjusted_p[i] = summary(glm(as.formula(paste0("feature_present ~ outcome+",paste(covariate_subset, collapse="+"))), data = rf_input, family = 'binomial'))$coefficients[2,4]
    }
    gini.cv$wilcox_p[i] = wilcox.test(rf_input[,feature] ~ rf_input$outcome)$p.value
  }
  
  # calculate fdr-adjusted p values - presence or absence
  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj[!is.na(gini.cv$fisher_p)] = p.adjust(gini.cv$fisher_p[!is.na(gini.cv$fisher_p)], method="BH")
  gini.cv$lr_unadjusted_padj = NA
  gini.cv$lr_unadjusted_padj[!is.na(gini.cv$lr_unadjusted_p)] = p.adjust(gini.cv$lr_unadjusted_p[!is.na(gini.cv$lr_unadjusted_p)], method="BH")
  gini.cv$lr_adjusted_padj = NA
  gini.cv$lr_adjusted_padj[!is.na(gini.cv$lr_adjusted_p)] = p.adjust(gini.cv$lr_adjusted_p[!is.na(gini.cv$lr_adjusted_p)], method="BH")
  gini.cv$wilcox_padj = NA
  gini.cv$wilcox_padj = p.adjust(gini.cv$wilcox_p, method="BH")
  
  # calculate prevalence differences
  gini.cv$mean_diff = gini.cv$responder_mean-gini.cv$nonresponder_mean
  gini.cv$prev_diff = gini.cv$responder_prev-gini.cv$nonresponder_prev
  gini.cv$enriched_responder_mean = gini.cv$mean_diff>0
  gini.cv$enriched_responder_prev = gini.cv$prev_diff>0
  
  # calculate RF rank and top 20 features by importance score
  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  gini.cv$rf_n = n_per_group
  
  # add cross-val median + IQR
  gini.cv$crossval_median = median(stat.cv$accuracy)
  gini.cv$crossval_quartile_1 = quantile(stat.cv$accuracy)[2]
  gini.cv$crossval_quartile_3 = quantile(stat.cv$accuracy)[4]
  gini.cv$feature_id = rownames(gini.cv) 
  
  # add maaslin2 analyses
  maaslin_features = data.frame(t(otu_table(rf_reference_rps)))
  maaslin_metadata = data.frame(sample_data(rf_reference_rps))
  
  # run unadjusted analysis
  fit_data_unadjusted = Maaslin2(
    input_data = maaslin_features, 
    input_metadata = maaslin_metadata, 
    output = paste0("maaslin_outputs/",outcome,"_",country_subset,"_unadjusted"), 
    fixed_effects = outcome,
    standardize = FALSE, normalization = "NONE", transform = "AST",
    min_prevalence = 0,
    max_significance = 0.2,
    plot_scatter = FALSE
  )
  results_unadjusted = fit_data_unadjusted$results
  names(results_unadjusted)[2:ncol(results_unadjusted)] = paste0("maaslin_", names(results_unadjusted)[2:ncol(results_unadjusted)])
  
  # run adjusted analysis
  fit_data_adjusted = Maaslin2(
    input_data = maaslin_features, 
    input_metadata = maaslin_metadata, 
    output = paste0("maaslin_outputs/",outcome,"_",country_subset,"_adjusted"), 
    fixed_effects = c(outcome, covariate_subset),
    standardize = FALSE, normalization = "NONE", transform = "AST",
    min_prevalence = 0,
    max_significance = 0.2,
    plot_scatter = FALSE
  )
  results_adjusted = fit_data_adjusted$results %>% filter(metadata == outcome)
  names(results_adjusted)[2:ncol(results_adjusted)] = paste0("maaslin_", names(results_adjusted)[2:ncol(results_adjusted)])
  names(results_adjusted)[2:ncol(results_adjusted)] = paste0(names(results_adjusted)[2:ncol(results_adjusted)],"_adjusted")
  
  # merge and write results
  if (all(results_unadjusted$feature %in% results_adjusted$feature)) { results_collated = left_join(results_unadjusted, results_adjusted, by = "feature") } else { print("error matching feature names") }
  
  # merge with RF results
  names(results_collated)[1] = "feature_id"
  results_collated = results_collated %>% dplyr::select(-c(maaslin_metadata, maaslin_value, maaslin_metadata_adjusted, maaslin_value_adjusted, 
                                                           maaslin_name_adjusted, maaslin_N_adjusted, maaslin_N.not.zero_adjusted))
  if (all(results_collated$feature_id %in% gini.cv$feature_id)) { gini.cv = merge(gini.cv, results_collated, by = "feature_id") }
  
  # reapply FDR adjustment to maaslin2 after filtering
  gini.cv$maaslin_qval_filt = p.adjust(gini.cv$maaslin_pval, method="BH")
  gini.cv$maaslin_qval_adjusted_filt = p.adjust(gini.cv$maaslin_pval_adjusted, method="BH")
  
  # add module and country data
  gini.cv$country = country_subset
  gini.cv$country_clean = factor(gini.cv$country, levels=c("Malawi", "India (unexposed)", "India (exposed)","India"))
  gini.cv$country_clean = revalue(as.factor(gini.cv$country_clean), c("Malawi" = "MLW", "India" = "IND", "India (exposed)" = "IND (neo+)", "India (unexposed)" = "IND (neo-)"))
  gini.cv$module = module
  
  # save RF results
  write.csv(stat.cv,file=paste0("RF_outputs/crossval_stats_",outcome,"_",country_subset,".csv"))
  write.csv(gini.cv, file=paste0("RF_outputs/crossval_gini_",outcome,"_",country_subset,".csv"))
  return(list(gini.cv=gini.cv))
} 

##########################
### CrossVal functions ###
##########################

group.samples = function(y) {
  # split samples into groups
  if(is.factor(y))
  {
    ygrouped = split(seq(y), y)
  }
  else
  {
    ygrouped = list(all=seq(length(y)))
  }
  
  return( ygrouped )
}

get.folds = function(ygrouped, K) {
  # makes no sense to have more folds than entries in largest group
  groupsize = sapply(ygrouped, length)
  nfolds = min(K, max(groupsize))
  if (K != nfolds) cat("Number of folds:", nfolds, "\n")
  
  # assign the members of each group evenly to the folds
  m = NULL
  for (i in 1:length(ygrouped) )
  {
    a = ceiling(groupsize[i]/nfolds)
    ridx = sample.int(groupsize[i], groupsize[i])
    v = c( rep(NA, nfolds*a-groupsize[i]), ygrouped[[i]][ridx] ) # pad with NAs
    ridx = sample.int(nfolds, nfolds) # reshuffle column containing the NAs
    v[1:nfolds] = v[ridx]
    m = c(m,v)
  }
  m = matrix(m, nrow=nfolds) # note that all NAs of a group are in one column
  
  folds =  vector("list", nfolds)
  for(j in 1:nfolds)
  {
    keep = !is.na(m[j , ])
    folds[[j]] = m[j, keep]
  }
  
  return( folds )
}

predfun = function(Xtrain, Ytrain, Xtest, Ytest){
  rf.fit = randomForest(Xtrain, y=Ytrain, ntree=500, importance=TRUE)
  ynew = predict(rf.fit, Xtest)
  # count false and true positives/negatives
  negative = levels(Ytrain)[2] 
  cm = confusionMatrix(Ytest, ynew, negative=negative)
  importance_Gini <- rf.fit$importance # update to original code: create dataframe during RF iteration that contains importance scores from RF model
  return(list(cm = cm, importance_Gini = importance_Gini)) # update to original code: return both confusion matrix and importance scores
}

crossval_update = function (predfun, X, Y, K = 5, B = 10, verbose = TRUE, ...) {
  ygrouped = group.samples(Y)
  groupsize = sapply(ygrouped, length)
  nfolds = min(K, max(groupsize))
  if (verbose)
    cat("Number of folds:", nfolds, "\n")
  allfolds = B * nfolds
  if (verbose)
    cat("Total number of CV fits:", allfolds, "\n")
  stat.cv = NULL
  Gini_importance_output = NULL # update to original code: output for importance scores
  i = 1
  for (b in 1:B) {
    if (verbose)
      cat("\nRound #", b, "of", B, "\n")
    folds = get.folds(ygrouped, K = nfolds)
    for (f in 1:nfolds) {
      if (verbose)
        cat("CV Fit #", i, "of", allfolds, "\n")
      test.idx = folds[[f]]
      train.x = X[-test.idx, , drop = FALSE]
      train.y = Y[-test.idx]
      test.x = X[test.idx, , drop = FALSE]
      test.y = Y[test.idx]
      stat.new = predfun(train.x, train.y, test.x, test.y,
                         ...)
      stat.cv = rbind(stat.cv, stat.new$cm) # update to original code: since predfun now returns both cm and importance, specify cm output here
      Gini_importance_output = cbind(Gini_importance_output, stat.new$importance_Gini) # update to original code: print importance scores from predfun output
      rownames(stat.cv)[i] = paste0("B", b, ".F", f)
      i = i + 1
    }
  }
  stat = apply(stat.cv, 2, mean)
  stat.se = apply(stat.cv, 2, sd)/sqrt(allfolds)
  return(list(stat.cv = stat.cv, stat = stat, stat.se = stat.se, Gini_importance_output = Gini_importance_output))
}



crossval_regression = function (predfun_regression, X, Y, K = 10, B = 20, verbose = TRUE, ...) {
  ygrouped = group.samples(Y)
  groupsize = sapply(ygrouped, length)
  nfolds = min(K, max(groupsize))
  if (verbose)
    cat("Number of folds:", nfolds, "\n")
  allfolds = B * nfolds
  if (verbose)
    cat("Total number of CV fits:", allfolds, "\n")
  stat.cv = NULL
  IncMSE_importance_output = NULL # update to original code: output for importance scores
  i = 1
  for (b in 1:B) {
    if (verbose)
      cat("\nRound #", b, "of", B, "\n")
    folds = get.folds(ygrouped, K = nfolds)
    for (f in 1:nfolds) {
      if (verbose)
        cat("CV Fit #", i, "of", allfolds, "\n")
      test.idx = folds[[f]]
      train.x = X[-test.idx, , drop = FALSE]
      train.y = Y[-test.idx]
      test.x = X[test.idx, , drop = FALSE]
      test.y = Y[test.idx]
      stat.new = predfun_regression(train.x, train.y, test.x, test.y)
      stat.cv = rbind(stat.cv, stat.new)
      IncMSE_importance_output = cbind(IncMSE_importance_output, stat.new$importance_IncMSE) # update to original code: print importance scores from predfun_regression output
      rownames(stat.cv)[i] = paste0("B", b, ".F", f)
      i = i + 1
    }
  }
  return(list(stat.cv = stat.cv, IncMSE_importance_output = IncMSE_importance_output)) # update to original code
}

# RF regression function
predfun_regression = function(Xtrain, Ytrain, Xtest, Ytest)
{
  rf.fit = randomForest(x = Xtrain, y = Ytrain, ntree=500, importance=TRUE)
  
  #statistics for training set
  tr_msr = mean((rf.fit$pred-rf.fit$y)^2)
  tr_perc_var_expl = 1-sum((rf.fit$y-rf.fit$pred)^2) /sum((rf.fit$y-mean(rf.fit$y))^2)
  tr_lm_R2  = round(summary(lm(rf.fit$pred~rf.fit$y))$r.squared,4)
  tr_pearson = cor(rf.fit$y, rf.fit$pred)
  tr_spearman = cor(rf.fit$y, rf.fit$pred, method="spearman")
  tr_n = length(Ytrain)
  
  #statistics for test set
  ynew = predict(rf.fit, Xtest)
  tst_msr = mean((Ytest - ynew)^2)
  tst_perc_var_expl = 1-sum((Ytest - ynew)^2) /sum((Ytest - mean(Ytest))^2)
  tst_lm_R2  = round(summary(lm(ynew ~ Ytest))$r.squared,4)
  tst_pearson = cor(Ytest, ynew)
  tst_spearman = cor(Ytest, ynew, method="spearman")
  tst_n = length(Ytest)
  
  # importance scores
  importance_IncMSE <- rf.fit$importance[,1]
  
  return(list(tr_msr = tr_msr, tr_perc_var_expl = tr_perc_var_expl, tr_lm_R2 = tr_lm_R2, tr_pearson = tr_pearson, tr_spearman = tr_spearman, tr_n = tr_n,
              tst_msr = tst_msr, tst_perc_var_expl = tst_perc_var_expl, tst_lm_R2 = tst_lm_R2, tst_pearson = tst_pearson, tst_spearman = tst_spearman, tst_n = tst_n,
              importance_IncMSE = importance_IncMSE))
}

# function to run Random Forests (run after deseq results so that list of taxa can be transferred)
collate_rf_IgA = function(output_name) {
  # create dataframe for collated accuracy results
  collated_rf = data.frame(country = rep(NA,length(country_list)))
  collated_rf$q1_tst_lm_R2 = collated_rf$q3_tst_lm_R2  = collated_rf$median_tst_lm_R2 = NA
  
  # run RF in loop across sample types and country-country comparisons
  for (c in 1:length(country_list)) {
    input_m$outcome = log(input_m[,"BB2_IgA"])
    country_subset = country_list[c]
    
    if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
      samplelist = subset(input_m, country_exposure==country_subset & !is.na(outcome))$sample_ID
    } else {
      samplelist = subset(input_m, country==country_subset & !is.na(outcome))$sample_ID
    }
    
    input_rps_prefiltering = prune_samples(samplelist, rps)
    input_rps = filter_taxa(input_rps_prefiltering, function(x) { sum(x > 0) > 0.05*nsamples(input_rps_prefiltering) }, TRUE)
    rf_input = data.frame(t(otu_table(input_rps))) # Create matrix of OTU abundances 

    # run RF using updated crossval_regression function
    log_IgA = log(as.numeric(as(sample_data(input_rps), "matrix")[,"BB2_IgA"]))
    K = 5 # number of folds
    B = 20 # number of repetitions
    cv.out = crossval_regression(predfun = predfun_regression, X = rf_input, Y = log_IgA, K=K, B=B, verbose = FALSE)
    
    # collate cross-validation statistics - stats
    stats = data.frame(cv.out$stat.cv[,1:(ncol(cv.out$stat.cv)-1)])
    stats_df = data.frame(tr_lm_R2 = unlist(stats$tr_lm_R2), tr_spearman = unlist(stats$tr_spearman), tr_n = unlist(stats$tr_n),
                          tst_lm_R2 = unlist(stats$tst_lm_R2), tst_spearman = unlist(stats$tst_spearman), tst_n = unlist(stats$tst_n))
    write.csv(stats_df,file=paste0("RF_outputs/crossval_stats_IgA_",country_subset,".csv"))
    
    # calculate importance score mean and SD across 1000 iterations
    IncMSE_output <- data.frame(cv.out$IncMSE_importance_output)
    for (i in 1:(length(rf_input))) {
      IncMSE_output$IncMSE_means[i] <- mean(as.numeric(IncMSE_output[i,1:B*K]))
      IncMSE_output$IncMSE_sds[i] <- sd(as.numeric(IncMSE_output[i,1:B*K]))
    }
    
    # drop mode_delivery for Malawi
    if (country_subset=="Malawi") { covariate_subset = maaslin_covariates[maaslin_covariates!="mode_delivery"] } else { covariate_subset = maaslin_covariates }
    
    # append prevalence, abundance data and Fisher's test p values for presence/absence
    for (i in 1:nrow(IncMSE_output)) {
      feature = rownames(IncMSE_output)[i]
      IncMSE_output$present[i] = sum(rf_input[,feature]>0)
      IncMSE_output$n[i] = nrow(rf_input)
      IncMSE_output$prev[i] = round(sum(rf_input[,feature]>0)/nrow(rf_input)*100,1)
      IncMSE_output$mean[i] = mean(rf_input[,feature])*100
      IncMSE_output$rho[i] = cor.test(log_IgA, rf_input[,feature], method="spearman")$estimate
      IncMSE_output$spearman_p[i] = cor.test(log_IgA, rf_input[,feature], method="spearman")$p.value
      if (IncMSE_output$prev[i]<5 | IncMSE_output$prev[i]>95) {
       IncMSE_output$wilcox_p[i] = NA
       IncMSE_output$mean_IgA_present[i] = NA
       IncMSE_output$CI_IgA_present[i] = NA
       IncMSE_output$mean_IgA_absent[i] = NA
       IncMSE_output$CI_IgA_absent[i] = NA
       IncMSE_output$lr_unadjusted_p[i] = NA
       IncMSE_output$lr_adjusted_p[i] = NA
       IncMSE_output$GMR[i] = NA
      } else {
       IncMSE_output$wilcox_p[i] = wilcox.test(log_IgA ~ rf_input[,feature]>0)$p.value
       
       # create glm dataframe
       rf_input$feature_present = as.numeric(rf_input[,feature]>0)
       rf_input$logIgA = as.numeric(as.numeric(as(sample_data(input_rps), "matrix")[,"BB2_IgA"]))
       rf_input$age_at_first_dose = as.numeric(sample_data(input_rps)$age_at_first_dose)
       rf_input$antibiotic_exposure = factor(sample_data(input_rps)$antibiotic_exposure)
       rf_input$birth_weight = as.numeric(sample_data(input_rps)$birth_weight)
       rf_input$breastfed_child = factor(sample_data(input_rps)$breastfed_child)
       rf_input$mode_delivery = factor(sample_data(input_rps)$mode_delivery)
       
       # calculate outputs
       gmean_present = Gmean(rf_input$logIgA[rf_input$feature_present==1], conf.level =0.95)
       IncMSE_output$mean_IgA_present[i] = round(as.numeric(gmean_present[1]),3)
       IncMSE_output$CI_IgA_present[i] = paste0(round(as.numeric(gmean_present[2]),3),"-",round(as.numeric(gmean_present[3]),3))
       
       gmean_absent = Gmean(rf_input$logIgA[rf_input$feature_present==0], conf.level=0.95)
       IncMSE_output$mean_IgA_absent[i] = round(as.numeric(gmean_absent[1]),3)
       IncMSE_output$CI_IgA_absent[i] = paste0(round(as.numeric(gmean_absent[2]),3),"-",round(as.numeric(gmean_absent[3]),3))

       IncMSE_output$GMR[i] = IncMSE_output$mean_IgA_present[i]/IncMSE_output$mean_IgA_absent[i]
       
       IncMSE_output$lr_unadjusted_p[i] = summary(glm(feature_present ~ logIgA, data = rf_input, family = 'binomial'))$coefficients[2,4]
       IncMSE_output$lr_adjusted_p[i] = summary(glm(as.formula(paste0("feature_present ~ logIgA+",paste(covariate_subset, collapse="+"))), data = rf_input, family = 'binomial'))$coefficients[2,4]
      }
    }
    
    # calculate fdr-adjusted p values - presence or absence
    IncMSE_output$spearman_padj = p.adjust(IncMSE_output$spearman_p, method="BH")
    IncMSE_output$wilcox_padj = NA
    IncMSE_output$wilcox_padj[!is.na(IncMSE_output$wilcox_p)] = p.adjust(IncMSE_output$wilcox_p[!is.na(IncMSE_output$wilcox_p)], method="BH")
    IncMSE_output$lr_unadjusted_padj = NA
    IncMSE_output$lr_unadjusted_padj[!is.na(IncMSE_output$lr_unadjusted_p)] = p.adjust(IncMSE_output$lr_unadjusted_p[!is.na(IncMSE_output$lr_unadjusted_p)], method="BH")
    IncMSE_output$lr_adjusted_padj = NA
    IncMSE_output$lr_adjusted_padj[!is.na(IncMSE_output$lr_adjusted_p)] = p.adjust(IncMSE_output$lr_adjusted_p[!is.na(IncMSE_output$lr_adjusted_p)], method="BH")
    IncMSE_output$feature_id = rownames(IncMSE_output)
    
    # add pathway details
    features = data.frame(tax_table(input_rps))
    features$feature_id = rownames(features)
    IncMSE_output = merge(features, IncMSE_output, by = "feature_id")
    IncMSE_output_sorted <- IncMSE_output[order(-IncMSE_output$IncMSE_means),] 
    
    # calculate RF rank and top 20 features by importance score
    IncMSE_output_sorted$rf_rank = 1:nrow(IncMSE_output_sorted)
    IncMSE_output_sorted$rf_top20 = as.numeric(IncMSE_output_sorted$rf_rank<=20)

    # add maaslin2 analyses
    maaslin_features = data.frame(t(otu_table(input_rps)))
    maaslin_metadata = data.frame(sample_data(input_rps))

    # run unadjusted analysis
    fit_data_unadjusted = Maaslin2(
      input_data = maaslin_features, 
      input_metadata = maaslin_metadata, 
      output = paste0("maaslin_outputs/IgA_",country_subset,"_unadjusted"), 
      fixed_effects = "IgA",
      standardize = FALSE, normalization = "NONE", transform = "AST",
      min_prevalence = 0,
      max_significance = 0.2,
      plot_scatter = FALSE
    )
    results_unadjusted = fit_data_unadjusted$results
    names(results_unadjusted)[2:ncol(results_unadjusted)] = paste0("maaslin_",names(results_unadjusted)[2:ncol(results_unadjusted)])
    
    # run adjusted analysis
    fit_data_adjusted = Maaslin2(
      input_data = maaslin_features, 
      input_metadata = maaslin_metadata, 
      output = paste0("maaslin_outputs/IgA_",country_subset,"_adjusted"), 
      fixed_effects = c("IgA", covariate_subset),
      standardize = FALSE, normalization = "NONE", transform = "AST",
      min_prevalence = 0,
      max_significance = 0.2,
      plot_scatter = FALSE
    )
    results_adjusted = fit_data_adjusted$results %>% filter(metadata == "IgA")
    names(results_adjusted)[2:ncol(results_adjusted)] = paste0("maaslin_",names(results_adjusted)[2:ncol(results_adjusted)])
    names(results_adjusted)[2:ncol(results_adjusted)] = paste0(names(results_adjusted)[2:ncol(results_adjusted)],"_adjusted")
    
    # merge and write results
    if (all(results_unadjusted$feature %in% results_adjusted$feature)) { results_collated = left_join(results_unadjusted, results_adjusted, by = "feature") } else { print("error matching feature names") }
    
    # merge with RF results
    names(results_collated)[1] = "feature_id"
    results_collated = results_collated %>% dplyr::select(-c(maaslin_metadata, maaslin_value, maaslin_metadata_adjusted, maaslin_value_adjusted, 
                                                             maaslin_name_adjusted, maaslin_N_adjusted, maaslin_N.not.zero_adjusted))
    if (all(results_collated$feature_id %in% IncMSE_output_sorted$feature_id)) { IncMSE_output_sorted = merge(IncMSE_output_sorted, results_collated, by = "feature_id") }
    
    # reapply FDR adjustment to maaslin2 after filtering
    IncMSE_output_sorted$maaslin_qval_filt = p.adjust(IncMSE_output_sorted$maaslin_pval, method="BH")
    IncMSE_output_sorted$maaslin_qval_adjusted_filt = p.adjust(IncMSE_output_sorted$maaslin_pval_adjusted, method="BH")
    
    # add module and country data
    IncMSE_output_sorted$country = country_subset
    IncMSE_output_sorted$country_clean = factor(IncMSE_output_sorted$country, levels=c("Malawi", "India (unexposed)", "India (exposed)","India"))
    IncMSE_output_sorted$country_clean = revalue(as.factor(IncMSE_output_sorted$country_clean), c("Malawi" = "MLW", "India" = "IND", "India (exposed)" = "IND (neo+)", "India (unexposed)" = "IND (neo-)"))
    IncMSE_output_sorted$module = module
    
    # save RF results
    write.csv(IncMSE_output_sorted,file=paste0("RF_outputs/crossval_regression_IgA_",country_subset,".csv"))
    
    # add RF statistics to collated_rf dataframe
    collated_rf$country[c] = country_subset
    collated_rf$median_tst_lm_R2[c] = median(unlist(stats$tst_lm_R2))
    collated_rf$q1_tst_lm_R2[c] = quantile(unlist(stats$tst_lm_R2))[2]
    collated_rf$q3_tst_lm_R2[c] = quantile(unlist(stats$tst_lm_R2))[4]
    
  }
  write.csv(collated_rf,file=paste0("RF_outputs/",output_name))
}

# function to generate reverse log10 scale
# from: https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
