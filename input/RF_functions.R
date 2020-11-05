###################################
### Script to load RF functions ###
###################################

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
    rf.fit = randomForest(Xtrain, y=Ytrain, ntree=5000, importance=TRUE)
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





###########################
### country RF function ###
###########################

run_rf = function(input_rps, prev_threshold=0.01, folds=5, iterations=iter, n_per_group=50) {

  # determine level - genus or rsv
  if (all(is.na(tax_table(input_rps)[,"tax_id"]))) { level = "genus"}
  if (all(!is.na(tax_table(input_rps)[,"tax_id"]))) { level = "rsv"}
  
  # determine taxa present above threshold in at least one country subset
  ps1list = sample_names(input_rps)[sample_data(input_rps)$country_arm==country1]
  ps1 = prune_samples(ps1list, input_rps)
  ps1 = filter_taxa(ps1, function(x) { sum(x > 0) > prev_threshold*nsamples(ps1) }, TRUE)
  ps2list = sample_names(input_rps)[sample_data(input_rps)$country_arm==country2]
  ps2 = prune_samples(ps2list, input_rps)
  ps2 = filter_taxa(ps2, function(x) { sum(x > 0) > prev_threshold*nsamples(ps2) }, TRUE)
  taxlist = unique(c(rownames(tax_table(ps1)),rownames(tax_table(ps2))))
  
  # if <50 in either country subset, amend n_per_group to match minimum group size
  mingroup = min(c(length(ps1list), length(ps2list)))
  if (mingroup < 50) { n_per_group = mingroup }
  
  # create RF input
  rf_reference_ps = prune_taxa(taxlist, input_rps)
  rf_input = data.frame(t(otu_table(rf_reference_ps))) 

  # append taxonomy assignment
  if (level == "genus") { names(rf_input) = tax_table(rf_reference_ps)[,"Genus"] }
  if (level == "rsv") { names(rf_input) = tax_table(rf_reference_ps)[,"tax_id"] }
  
  # create country-specific subsets
  rf_input$country_arm = factor(sample_data(rf_reference_ps)$country_arm)
  rf_input_1 = subset(rf_input, country_arm==as.character(country1))
  rf_input_2 = subset(rf_input, country_arm==as.character(country2))
  
  # perform crossval on whole dataset
  crossval_stats = list()
  crossval_importance = list()
  
  for (i in 1:iter) {
    # random draw for each iteration
    rf_s = rbind(rf_input_1[sample(nrow(rf_input_1),n_per_group),],rf_input_2[sample(nrow(rf_input_2),n_per_group),])
    # perform crossval for iteration above
    cv.rf = crossval_update(predfun, X=rf_s[,!grepl("country_arm", names(rf_s))], Y=factor(rf_s$country_arm), K=folds, B=1, verbose=FALSE)
    # write outputs
    crossval_stats[[i]] = cv.rf$stat.cv
    crossval_importance[[i]] = cv.rf$Gini_importance_output
  }
  
  # collate cross-validation statistics
  stat.cv = as.data.frame(crossval_stats[[1]])
  for (i in 2:iter) { stat.cv = rbind(stat.cv, as.data.frame(crossval_stats[[i]])) }
  stat.cv$accuracy = (stat.cv$TP+stat.cv$TN)/(stat.cv$TP+stat.cv$TN+stat.cv$FP+stat.cv$FN)
    
  # collate cross-validation statistics - stats
  gini.cv = as.data.frame(crossval_importance)
  gini.cv$mean = rowMeans(gini.cv[,1:(iter*folds)])
  gini.cv$sd = rowSds(as.matrix(gini.cv[,1:(iter*folds)]))
  gini.cv = gini.cv[,(ncol(gini.cv)-1):ncol(gini.cv)]
  
  # append taxonomic assignments - gini
  if (level == "genus") { gini.cv$taxonomy = gini.cv$tax_id = rownames(gini.cv) }
  if (level == "rsv") { 
    gini.cv$tax_id = rownames(gini.cv) 
    gini.cv = merge(gini.cv, as(tax_table(rf_reference_ps), "matrix")[,c("tax_id", "taxonomy")], by = "tax_id", keep.x = TRUE) }
  gini.cv <- gini.cv[order(-gini.cv$mean),]
  gini.cv$taxonomy <- factor(gini.cv$taxonomy, levels = gini.cv$taxonomy)
  
  # append prevalence, abundance data and Fisher's test p values for presence/absence
  for (i in 1:nrow(gini.cv)) {
    taxon = gini.cv$tax_id[i]
    gini.cv$group1_country[i] = as.character(country1)
    gini.cv$group1_present[i] = sum(rf_input_1[,taxon]>0)
    gini.cv$group1_present_strict[i] = sum(rf_input_1[,taxon]>=0.001)
    gini.cv$group1_n[i] = nrow(rf_input_1)
    gini.cv$group1_prev[i] = sum(rf_input_1[,taxon]>0)/nrow(rf_input_1)
    gini.cv$group1_prev_strict[i] = sum(rf_input_1[,taxon]>=0.001)/nrow(rf_input_1)
    gini.cv$group1_mean[i] = mean(rf_input_1[,taxon])
    gini.cv$group1_sd[i] = sd(rf_input_1[,taxon])
    
    gini.cv$group2_country[i] = as.character(country2)
    gini.cv$group2_present[i] = sum(rf_input_2[,taxon]>0)
    gini.cv$group2_present_strict[i] = sum(rf_input_2[,taxon]>=0.001)
    gini.cv$group2_n[i] = nrow(rf_input_2)
    gini.cv$group2_prev[i] = sum(rf_input_2[,taxon]>0)/nrow(rf_input_2)
    gini.cv$group2_prev_strict[i] = sum(rf_input_2[,taxon]>=0.001)/nrow(rf_input_2)
    gini.cv$group2_mean[i] = mean(rf_input_2[,taxon])
    gini.cv$group2_sd[i] = sd(rf_input_2[,taxon])
    gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>0), sum(rf_input_1[,taxon]==0),
                                               sum(rf_input_2[,taxon]>0), sum(rf_input_2[,taxon]==0)),nrow=2))$p.value
    gini.cv$fisher_p_strict[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>=0.001), sum(rf_input_1[,taxon]<0.001),
                                               sum(rf_input_2[,taxon]>=0.001), sum(rf_input_2[,taxon]<0.001)),nrow=2))$p.value
    
  }
  
  # calculate fdr-adjusted p values - presence or absence
  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj[!(gini.cv$group1_prev<0.05 & gini.cv$group2_prev<0.05) & !(gini.cv$group1_prev>0.95 & gini.cv$group2_prev>0.95)] = 
    p.adjust(gini.cv$fisher_p[!(gini.cv$group1_prev<0.05 & gini.cv$group2_prev<0.05) & !(gini.cv$group1_prev>0.95 & gini.cv$group2_prev>0.95)], method="BH")
  gini.cv$fisher_signif = as.numeric(gini.cv$fisher_padj<=0.1)
  
  # calculate fdr-adjusted p values - presence or absence at ≥0.1% relative abundance
  gini.cv$fisher_padj_strict = NA
  gini.cv$fisher_padj_strict[!(gini.cv$group1_prev_strict<0.05 & gini.cv$group2_prev_strict<0.05) & !(gini.cv$group1_prev_strict>0.95 & gini.cv$group2_prev_strict>0.95)] = 
    p.adjust(gini.cv$fisher_p_strict[!(gini.cv$group1_prev_strict<0.05 & gini.cv$group2_prev_strict<0.05) & !(gini.cv$group1_prev_strict>0.95 & gini.cv$group2_prev_strict>0.95)], method="BH")
  gini.cv$fisher_signif_strict = as.numeric(gini.cv$fisher_padj_strict<=0.1)
  
  # calculate prevalence differences
  gini.cv$prev_diff = gini.cv$group1_prev-gini.cv$group2_prev
  gini.cv$prev_diff_strict = gini.cv$group1_prev_strict-gini.cv$group2_prev_strict
  gini.cv$enriched_group1 = gini.cv$prev_diff>0
  gini.cv$enriched_group1_strict = gini.cv$prev_diff_strict>0
  
  # calculate RF rank and top 20 features by importance score
  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  
  # add plot colour scheme
  gini.cv$fisher_col = "ns"
  gini.cv$fisher_col[gini.cv$prev_diff>0 & gini.cv$fisher_padj<0.1] = country1
  gini.cv$fisher_col[gini.cv$prev_diff<0 & gini.cv$fisher_padj<0.1] = country2
  
  gini.cv$fisher_col_strict = "ns"
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict>0 & gini.cv$fisher_padj_strict<0.1] = country1
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict<0 & gini.cv$fisher_padj_strict<0.1] = country2
  
  # add cross-val mean + sd
  gini.cv$crossval_mean = mean(stat.cv$accuracy)
  gini.cv$crossval_sd = sd(stat.cv$accuracy)
  
  # save RF results
  write.csv(stat.cv,file=paste0("output_module3/RF_output_stats/crossval_stats_",country1,"_",country2,"_",subset,"_",level,".csv"))
  write.csv(gini.cv, file=paste0("output_module3/RF_output_importance/crossval_gini_",country1,"_",country2,"_",subset,"_",level,".csv"))
  return(list(gini.cv=gini.cv))
} 





#######################################
### ORV binary outcome RF functions ###
#######################################

### function to run Random Forests (run after deseq results so that list of taxa can be transferred)
run_rf_binary_outcome = function(input_rps, country_arm, sample_sub, prev_threshold=0.01, folds=5, iterations=iter, n_per_group=50) {
  
  # determine level - genus or rsv
  if (all(is.na(tax_table(input_rps)[,"tax_id"]))) { level = "genus"}
  if (all(!is.na(tax_table(input_rps)[,"tax_id"]))) { level = "rsv"}
  
  # determine taxa present above threshold in at least one response subset
  # nonresponders
  ps1list = sample_names(input_rps)[sample_data(input_rps)$outcome==0] 
  ps1 = prune_samples(ps1list, input_rps)
  ps1 = filter_taxa(ps1, function(x) { sum(x > 0) > prev_threshold*nsamples(ps1) }, TRUE)
  
  # responders
  ps2list = sample_names(input_rps)[sample_data(input_rps)$outcome==1]
  ps2 = prune_samples(ps2list, input_rps)
  ps2 = filter_taxa(ps2, function(x) { sum(x > 0) > prev_threshold*nsamples(ps2) }, TRUE)
  
  # pick out taxa present in at least one of the lists above
  taxlist = unique(c(rownames(tax_table(ps1)),rownames(tax_table(ps2))))
  
  # if <50 responders or <50 non-responders, amend n_per_group to match minimum group size
  mingroup = min(c(length(ps1list), length(ps2list)))
  if (mingroup < 50) { n_per_group = mingroup }
  
  # create RF input
  rf_reference_ps = prune_taxa(taxlist, input_rps)
  rf_input = data.frame(t(otu_table(rf_reference_ps))) # Create matrix of OTU abundances 

  # append rsv or genus assignment
  if (level == "genus") { names(rf_input) = tax_table(rf_reference_ps)[,"Genus"] }
  if (level == "rsv") { names(rf_input) = tax_table(rf_reference_ps)[,"tax_id"] }
  
  # create response-specific subsets
  rf_input$outcome = factor(sample_data(rf_reference_ps)$outcome)
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
  
  # append taxonomic assignments
  if (level == "genus") { gini.cv$taxonomy = gini.cv$tax_id = rownames(gini.cv) }
  if (level == "rsv") { 
    gini.cv$tax_id = rownames(gini.cv) 
    gini.cv = merge(gini.cv, as(tax_table(rf_reference_ps), "matrix")[,c("tax_id", "taxonomy")], by = "tax_id", keep.x = TRUE) }
  gini.cv <- gini.cv[order(-gini.cv$mean),]
  gini.cv$taxonomy <- factor(gini.cv$taxonomy, levels = gini.cv$taxonomy)
  
  # append prevalence, abundance data and Fisher's test p values for presence/absence
  for (i in 1:nrow(gini.cv)) {
    taxon = gini.cv$tax_id[i]
    gini.cv$nonresponder_present[i] = sum(rf_input_1[,taxon]>0)
    gini.cv$nonresponder_present_strict[i] = sum(rf_input_1[,taxon]>=0.001)
    gini.cv$nonresponder_n[i] = nrow(rf_input_1)
    gini.cv$nonresponder_prev[i] = sum(rf_input_1[,taxon]>0)/nrow(rf_input_1)
    gini.cv$nonresponder_prev_strict[i] = sum(rf_input_1[,taxon]>=0.001)/nrow(rf_input_1)
    gini.cv$nonresponder_mean[i] = mean(rf_input_1[,taxon])
    gini.cv$nonresponder_sd[i] = sd(rf_input_1[,taxon])
    
    gini.cv$responder_present[i] = sum(rf_input_2[,taxon]>0)
    gini.cv$responder_present_strict[i] = sum(rf_input_2[,taxon]>=0.001)
    gini.cv$responder_n[i] = nrow(rf_input_2)
    gini.cv$responder_prev[i] = sum(rf_input_2[,taxon]>0)/nrow(rf_input_2)
    gini.cv$responder_prev_strict[i] = sum(rf_input_2[,taxon]>=0.001)/nrow(rf_input_2)
    gini.cv$responder_mean[i] = mean(rf_input_2[,taxon])
    gini.cv$responder_sd[i] = sd(rf_input_2[,taxon])
    gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>0), sum(rf_input_1[,taxon]==0),
                                               sum(rf_input_2[,taxon]>0), sum(rf_input_2[,taxon]==0)),nrow=2))$p.value
    gini.cv$fisher_p_strict[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>=0.001), sum(rf_input_1[,taxon]<0.001),
                                                      sum(rf_input_2[,taxon]>=0.001), sum(rf_input_2[,taxon]<0.001)),nrow=2))$p.value
  }
  
  # calculate fdr-adjusted p values - presence or absence
  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj[!(gini.cv$nonresponder_prev<0.05 & gini.cv$responder_prev<0.05) & !(gini.cv$nonresponder_prev>0.95 & gini.cv$responder_prev>0.95)] = 
    p.adjust(gini.cv$fisher_p[!(gini.cv$nonresponder_prev<0.05 & gini.cv$responder_prev<0.05) & !(gini.cv$nonresponder_prev>0.95 & gini.cv$responder_prev>0.95)], method="fdr")
  gini.cv$fisher_signif = as.numeric(gini.cv$fisher_padj<=0.1)
  
  # calculate fdr-adjusted p values - presence or absence at ≥0.1% relative abundance
  gini.cv$fisher_padj_strict = NA
  gini.cv$fisher_padj_strict[!(gini.cv$nonresponder_prev_strict<0.05 & gini.cv$responder_prev_strict<0.05) & !(gini.cv$nonresponder_prev_strict>0.95 & gini.cv$responder_prev_strict>0.95)] = 
    p.adjust(gini.cv$fisher_p_strict[!(gini.cv$nonresponder_prev_strict<0.05 & gini.cv$responder_prev_strict<0.05) & !(gini.cv$nonresponder_prev_strict>0.95 & gini.cv$responder_prev_strict>0.95)], method="fdr")
  gini.cv$fisher_signif_strict = as.numeric(gini.cv$fisher_padj_strict<=0.1)
  
  # calculate prevalence differences
  gini.cv$prev_diff = gini.cv$nonresponder_prev-gini.cv$responder_prev
  gini.cv$prev_diff_strict = gini.cv$nonresponder_prev_strict-gini.cv$responder_prev_strict
  gini.cv$enriched_nonresponder = gini.cv$prev_diff>0
  gini.cv$enriched_nonresponder_strict = gini.cv$prev_diff_strict>0
  
  # calculate RF rank and top 20 features by importance score
  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  gini.cv$rf_n = n_per_group
  
  # add plot colour scheme
  gini.cv$fisher_col = "ns"
  gini.cv$fisher_col[gini.cv$prev_diff>0 & gini.cv$fisher_padj<0.1] = "nonresponder"
  gini.cv$fisher_col[gini.cv$prev_diff<0 & gini.cv$fisher_padj<0.1] = "responder"
  
  gini.cv$fisher_col_strict = "ns"
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict>0 & gini.cv$fisher_padj_strict<0.1] = "nonresponder"
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict<0 & gini.cv$fisher_padj_strict<0.1] = "responder"
  
  # add cross-val mean + sd
  gini.cv$crossval_mean = mean(stat.cv$accuracy)
  gini.cv$crossval_sd = sd(stat.cv$accuracy)
  
  # save RF results
  write.csv(stat.cv,file=paste0("output_module4/RF_output_stats/crossval_stats_",outcome,"_",country_arm,"_",sample_sub,"_",level,".csv"))
  write.csv(gini.cv, file=paste0("output_module4/RF_output_importance/crossval_gini_",outcome,"_",country_arm,"_",sample_sub,"_",level,".csv"))
  return(list(gini.cv=gini.cv))
} 

### function to run collection of Random Forests
collate_rf_binary = function(output_name) {
  # create dataframe for collated accuracy results
  collated_rf = data.frame(sample_subset = rep(NA,length(subset_list)*length(country_list)))
  collated_rf$responders = collated_rf$nonresponders = collated_rf$country = collated_rf$crossval_mean = collated_rf$crossval_sd = NA
  
  # run RF in loop across sample types and country-country comparisons
  for (s in 1:length(subset_list)) {
    for (c in 1:length(country_list)) {
      input_m$outcome = input_m[,outcome]
      sample_subset = subset_list[s]
      country_subset = country_list[c]
      
      if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
        samplelist = subset(input_m, country_exposure==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID_full
      } else {
        samplelist = subset(input_m, country==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID_full 
      }
      
      input_rps = prune_samples(samplelist, rps)
      sample_data(input_rps)$outcome = as(sample_data(input_rps), "matrix")[,outcome]
      rf = run_rf_binary_outcome(input_rps, country_arm=country_subset, sample_sub=sample_subset, 
                                 prev_threshold=0.01, folds=5, iterations=20, n_per_group=50)
      
      collated_rf$sample_subset[length(country_list)*(s-1)+c] = sample_subset
      collated_rf$country[length(country_list)*(s-1)+c] = country_subset
      collated_rf$nonresponders[length(country_list)*(s-1)+c] = as.character("nonresponders")
      collated_rf$responders[length(country_list)*(s-1)+c] = as.character("responders")
      collated_rf$crossval_mean[length(country_list)*(s-1)+c] = rf$gini.cv$crossval_mean[1]
      collated_rf$crossval_sd[length(country_list)*(s-1)+c] = rf$gini.cv$crossval_sd[1]
    }
  }
  write.csv(collated_rf,file=paste0("output_module4/",output_name))
}






###########################################
### ORV continuous outcome RF functions ###
###########################################

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
  rf.fit = randomForest(x = Xtrain, y = Ytrain, ntree=5000, importance=TRUE)
  
  #statistics for training set
  tr_msr = mean((rf.fit$pred-rf.fit$y)^2)
  tr_perc_var_expl = 1-sum((rf.fit$y-rf.fit$pred)^2) /sum((rf.fit$y-mean(rf.fit$y))^2)
  tr_lm_R2  = round(summary(lm(rf.fit$pred~rf.fit$y))$r.squared,4)
  tr_pearson = cor(rf.fit$y, rf.fit$pred)^2
  tr_spearman = cor(rf.fit$y, rf.fit$pred, method="s")^2
  tr_n = length(Ytrain)
  
  #statistics for test set
  ynew = predict(rf.fit, Xtest)
  tst_msr = mean((Ytest - ynew)^2)
  tst_perc_var_expl = 1-sum((Ytest - ynew)^2) /sum((Ytest - mean(Ytest))^2)
  tst_lm_R2  = round(summary(lm(ynew ~ Ytest))$r.squared,4)
  tst_pearson = cor(Ytest, ynew)^2
  tst_spearman = cor(Ytest, ynew, method="s")^2
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
  collated_rf = data.frame(sample_subset = rep(NA,length(subset_list)*length(country_list)))
  collated_rf$country = collated_rf$sd_tst_lm_R2  = collated_rf$mean_tst_lm_R2 = NA
  
  # run RF in loop across sample types and country-country comparisons
  for (s in 1:length(subset_list)) {
    for (c in 1:length(country_list)) {
      input_m$outcome = log(input_m[,"BB2_IgA"])
      sample_subset = subset_list[s]
      country_subset = country_list[c]
      
      if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
        samplelist = subset(input_m, country_exposure==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID_full
      } else {
        samplelist = subset(input_m, country==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID_full 
      }
      
      input_rps = prune_samples(samplelist, rps)
      input_rps = filter_taxa(input_rps, function(x) { sum(x > 0) > 0.01*nsamples(input_rps) }, TRUE)
      rf_input = data.frame(t(otu_table(input_rps))) # Create matrix of OTU abundances 

      # append rsv or genus assignment
      if (level=="rsv") { names(rf_input) = tax_table(input_rps)[,"taxonomy"] } else { names(rf_input) = tax_table(input_rps)[,"Genus"]  }
      
      # run RF using updated crossval_regression function
      log_IgA = log(as.numeric(as(sample_data(input_rps), "matrix")[,"BB2_IgA"]))
      K = 5 # number of folds
      B = 20 # number of repititions
      cv.out = crossval_regression(predfun = predfun_regression, X = rf_input, Y = log_IgA, K=K, B=B, verbose = FALSE)
      
      # collate cross-validation statistics - stats
      stats = data.frame(cv.out$stat.cv[,1:(ncol(cv.out$stat.cv)-1)])
      stats_df = data.frame(tr_lm_R2 = unlist(stats$tr_lm_R2), tr_spearman = unlist(stats$tr_lm_R2), tr_n = unlist(stats$tr_n),
                            tst_lm_R2 = unlist(stats$tst_lm_R2), tst_spearman = unlist(stats$tst_spearman), tst_n = unlist(stats$tst_n))
      write.csv(stats_df,file=paste0("output_module4/RF_output_stats/crossval_stats_IgA_",country_subset,"_",sample_subset,"_",level,".csv"))
      
      # calculate importance score mean and SD across 1000 iterations
      IncMSE_output <- data.frame(cv.out$IncMSE_importance_output)
      for (i in 1:(length(rf_input))) {
        IncMSE_output$IncMSE_means[i] <- mean(as.numeric(IncMSE_output[i,1:B*K]))
        IncMSE_output$IncMSE_sds[i] <- sd(as.numeric(IncMSE_output[i,1:B*K]))
      }
      IncMSE_output_sorted <- IncMSE_output[order(-IncMSE_output$IncMSE_means),] 
      write.csv(IncMSE_output_sorted,file=paste0("output_module4/RF_output_importance/crossval_regression_IgA_",country_subset,"_",sample_subset,"_",level,".csv"))
      
      # add RF statistics to collated_rf dataframe
      collated_rf$sample_subset[length(country_list)*(s-1)+c] = sample_subset
      collated_rf$country[length(country_list)*(s-1)+c] = country_subset
      collated_rf$mean_tst_lm_R2[length(country_list)*(s-1)+c] = mean(unlist(stats$tst_lm_R2))
      collated_rf$sd_tst_lm_R2[length(country_list)*(s-1)+c] = sd(unlist(stats$tst_lm_R2))
    }
  }
  write.csv(collated_rf,file=paste0("output_module4/",output_name))
}





############################
### Modular RF functions ###
############################

# binary outcomes
rf_modules_binary = function(module, module_file, outcome, outcome_file, iterations = iter, n_per_group=50, folds=5) {
  module$outcome = outcome
  
  # convert non-numerics and non-integers into factors
  for (i in 1:ncol(module)) {
    if(!(is.double(module[,i]) | is.integer(module[,i]))) { module[,i]<-as.factor(module[,i]) }
  }
  
  ### split into responders and non-responders
  rf_input_1 = subset(module, outcome==0)
  rf_input_2 = subset(module, outcome==1)
  
  # perform crossval on whole dataset
  crossval_stats = list()
  crossval_importance = list()
  
  # if <50 responders or <50 non-responders, amend n_per_group to match minimum group size
  mingroup = min(c(nrow(rf_input_1), nrow(rf_input_2)))
  if (mingroup < 50) { n_per_group = mingroup }
  
  for (i in 1:iter) {
    # random draw for each iteration
    rf_s = rbind(rf_input_1[sample(nrow(rf_input_1),n_per_group),],rf_input_2[sample(nrow(rf_input_2),n_per_group),])
    # perform crossval for iteration above
    cv.rf = crossval_update(predfun, X=rf_s[,!grepl("outcome", names(rf_s))], Y=factor(rf_s$outcome), K=folds, B=1, verbose=FALSE)
    # write outputs
    crossval_stats[[i]] = cv.rf$stat.cv
    crossval_importance[[i]] = cv.rf$Gini_importance_output
  }
  
  # collate crossvalidation statistics
  stat.cv = as.data.frame(crossval_stats[[1]])
  for (i in 2:iter) { stat.cv = rbind(stat.cv, as.data.frame(crossval_stats[[i]])) }
  stat.cv$accuracy = (stat.cv$TP+stat.cv$TN)/(stat.cv$TP+stat.cv$TN+stat.cv$FP+stat.cv$FN)
  write.csv(stat.cv,file=paste0("output_module5/RF_output_stats/crossval_stats_",outcome_file,"_",module_file,".csv"))
  
  # collate crossvalidation statistics
  gini.cv = as.data.frame(crossval_importance)
  gini.cv$mean = rowMeans(gini.cv[,1:(iter*folds)])
  gini.cv$sd = rowSds(as.matrix(gini.cv[,1:(iter*folds)]))
  gini.cv = gini.cv[,(ncol(gini.cv)-1):ncol(gini.cv)]
  gini.cv = gini.cv[order(-gini.cv$mean),]
  
  # add cross-val mean + sd
  gini.cv$crossval_mean = mean(stat.cv$accuracy)
  gini.cv$crossval_sd = sd(stat.cv$accuracy)
  gini.cv$module = module_file
  write.csv(gini.cv,file=paste0("output_module5/RF_output_importance/importance_",outcome_file,"_",module_file,".csv"))
  
  return(gini.cv)
}

# RV-IgA
rf_modules_IgA = function(module, module_file, K=5, B=20) {
  cv.out = crossval_regression(predfun = predfun_regression, X = module, Y = IgA, K=K, B=B, verbose = TRUE) 
  stats = data.frame(cv.out$stat.cv[,1:(ncol(cv.out$stat.cv)-1)])
  stats_df = data.frame(tr_lm_R2 = unlist(stats$tr_lm_R2), tr_spearman = unlist(stats$tr_lm_R2), tr_n = unlist(stats$tr_n),
                        tst_lm_R2 = unlist(stats$tst_lm_R2), tst_spearman = unlist(stats$tst_spearman), tst_n = unlist(stats$tst_n))
  stats_df$module = module_file
  write.csv(stats_df,file=paste0("output_module5/RF_output_stats/crossval_stats_IgA_",module_file,".csv"))
  
  # calculate importance score mean and SD across 1000 iterations
  IncMSE_output <- data.frame(cv.out$IncMSE_importance_output)
  for (i in 1:(length(module))) {
    IncMSE_output$IncMSE_means[i] <- mean(as.numeric(IncMSE_output[i,1:B*K]))
    IncMSE_output$IncMSE_sds[i] <- sd(as.numeric(IncMSE_output[i,1:B*K]))
  }
  IncMSE_output_sorted <- IncMSE_output[order(-IncMSE_output$IncMSE_means),] 
  #head(IncMSE_output_sorted, n = 10)
  write.csv(IncMSE_output_sorted,file=paste0("output_module5/RF_output_importance/importance_regression_IgA_",module_file,".csv"))
  return(stats_df)
}

# function to create gini summary plot
gini_plot = function(input_df, colour1=col1, colour2=col2, strict=TRUE, delivery=FALSE, species=FALSE) {
  if (delivery==TRUE) {
    input_df$enriched_group1 = input_df$enriched_vaginal
    input_df$enriched_group1_strict = input_df$enriched_vaginal_strict
  }
  
  input_df = input_df[order(input_df$rf_rank),]
  input_top20 = input_df[1:20,]
  
  if (species==TRUE) {
    input_top20$taxonomy = input_top20$taxonomy_species
  }
  
  input_top20$taxonomy <- factor(input_top20$taxonomy, levels = input_top20$taxonomy)
  
  if (strict==TRUE) { input_top20$enriched_group1 = input_top20$enriched_group1_strict }
  
  if (all(input_top20$enriched_group1==TRUE)) {
    ggplot(input_top20, aes(x = taxonomy, y = mean, fill=enriched_group1)) +
      scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
      xlab("") + ylab("importance (Gini)") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + scale_fill_manual(values=c(colour1))
  } else {
    ggplot(input_top20, aes(x = taxonomy, y = mean, fill=enriched_group1)) +
      scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
      xlab("") + ylab("importance (Gini)") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + scale_fill_manual(values=c(colour2,colour1))
  }
}

# function to create prevalence/abundance plot
rf_prev_plot = function(input_df=gini.cv, label1=country1, colour1=col1, label2=country2, colour2=col2, strict=TRUE, species=FALSE) {
  input_df = input_df[order(input_df$rf_rank),]
  input_top20 = input_df[1:20,]
  
  if (species==TRUE) {
    input_top20$taxonomy = input_top20$taxonomy_species
  }
  input_top20$taxonomy = factor(input_top20$taxonomy, levels = input_top20$taxonomy)
  
  if (strict==TRUE) {
    input_top20$group1_prev = input_top20$group1_prev_strict
    input_top20$group2_prev = input_top20$group2_prev_strict
  }
  
  ggplot(input_top20, aes(x = taxonomy, y = group1_prev)) +
    scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
    xlab("") + ylab("prevalence") +
    geom_point(aes(colour = label1, size = group1_mean), stat = "identity", alpha=0.7) +
    geom_point(aes(y = group2_prev, colour = label2, size = group2_mean), stat = "identity", alpha=0.7) +
    coord_flip() + theme_bw() + ylim(0,1) + scale_y_continuous(breaks=c(0,0.5,1)) +
    theme(axis.text.y = element_blank(), legend.position = "right") + labs(colour = "enriched group", size = "mean \nabundance") +
    scale_colour_manual(values=c(colour1,colour2)) + scale_size(limits = c(0, 0.5)) +
    guides(colour = guide_legend(order=1), size = guide_legend(order=2))
  #guides(colour = "none")
  
}































### function to run Random Forests (run after deseq results so that list of taxa can be transferred)
# run_rf_delivery = function(input_rps, prev_threshold=0.01, folds=5, iterations=iter, n_per_group=50) {
#   
#   # determine level
#   if (all(is.na(tax_table(input_rps)[,"tax_id"]))) { level = "genus"}
#   if (all(!is.na(tax_table(input_rps)[,"tax_id"]))) { level = "rsv"}
#   
#   # determine taxa present above threshold in at least one group
#   ps1list = sample_names(input_rps)[sample_data(input_rps)$mode_delivery=="vaginal"]
#   ps1 = prune_samples(ps1list, input_rps)
#   ps1 = filter_taxa(ps1, function(x) { sum(x > 0) > prev_threshold*nsamples(ps1) }, TRUE)
#   ps2list = sample_names(input_rps)[sample_data(input_rps)$mode_delivery=="caesarean"]
#   ps2 = prune_samples(ps2list, input_rps)
#   ps2 = filter_taxa(ps2, function(x) { sum(x > 0) > prev_threshold*nsamples(ps2) }, TRUE)
#   taxlist = unique(c(rownames(tax_table(ps1)),rownames(tax_table(ps2))))
#   
#   # if <50 responders or <50 non-responders, amend n_per_group to match minimum group size
#   mingroup = min(c(length(ps1list), length(ps2list)))
#   if (mingroup < 50) { n_per_group = mingroup }
#   
#   # create RF input
#   rf_reference_ps = prune_taxa(taxlist, input_rps)
#   rf_input = data.frame(t(otu_table(rf_reference_ps))) # Create matrix of OTU abundances 
#   #all(names(rf_input) == rownames(tax_table(rf_reference_ps)))
#   
#   # append rsv or genus assignment
#   if (level == "genus") { names(rf_input) = tax_table(rf_reference_ps)[,"Genus"] }
#   if (level == "rsv") { names(rf_input) = tax_table(rf_reference_ps)[,"tax_id"] }
#   
#   # all(rownames(rf_input) == rownames(sample_data(rf_reference_ps)))
#   rf_input$mode_delivery = factor(sample_data(rf_reference_ps)$mode_delivery)
#   
#   # create mode-specific subsets
#   rf_input_1 = subset(rf_input, mode_delivery=="vaginal")
#   rf_input_2 = subset(rf_input, mode_delivery=="caesarean")
#   
#   # perform crossval on whole dataset
#   crossval_stats = list()
#   crossval_importance = list()
#   
#   for (i in 1:iter) {
#     # random draw for each iteration
#     rf_s = rbind(rf_input_1[sample(nrow(rf_input_1),n_per_group),],rf_input_2[sample(nrow(rf_input_2),n_per_group),])
#     # perform crossval for iteration above
#     cv.rf = crossval_update(predfun, X=rf_s[,!grepl("mode_delivery", names(rf_s))], Y=factor(rf_s$mode_delivery), K=folds, B=1, verbose=FALSE)
#     # write outputs
#     crossval_stats[[i]] = cv.rf$stat.cv
#     crossval_importance[[i]] = cv.rf$Gini_importance_output
#   }
#   
#   # collate crossvalidation statistics
#   stat.cv = as.data.frame(crossval_stats[[1]])
#   for (i in 2:iter) { stat.cv = rbind(stat.cv, as.data.frame(crossval_stats[[i]])) }
#   stat.cv$accuracy = (stat.cv$TP+stat.cv$TN)/(stat.cv$TP+stat.cv$TN+stat.cv$FP+stat.cv$FN)
#   
#   # collate crossvalidation statistics
#   gini.cv = as.data.frame(crossval_importance)
#   gini.cv$mean = rowMeans(gini.cv[,1:(iter*folds)])
#   gini.cv$sd = rowSds(as.matrix(gini.cv[,1:(iter*folds)]))
#   gini.cv = gini.cv[,(ncol(gini.cv)-1):ncol(gini.cv)]
#   
#   if (level == "genus") { gini.cv$taxonomy = gini.cv$tax_id = rownames(gini.cv) }
#   if (level == "rsv") { 
#     gini.cv$tax_id = rownames(gini.cv) 
#     gini.cv = merge(gini.cv, as(tax_table(rf_reference_ps), "matrix")[,c("tax_id", "taxonomy", "taxonomy_species")], by = "tax_id", keep.x = TRUE) }
#   gini.cv <- gini.cv[order(-gini.cv$mean),]
#   gini.cv$taxonomy <- factor(gini.cv$taxonomy, levels = gini.cv$taxonomy)
#   
#   # append prevalence, abundance data and Fisher p
#   for (i in 1:nrow(gini.cv)) {
#     taxon = gini.cv$tax_id[i]
#     gini.cv$vaginal_present[i] = sum(rf_input_1[,taxon]>0)
#     gini.cv$vaginal_present_strict[i] = sum(rf_input_1[,taxon]>=0.001)
#     gini.cv$vaginal_n[i] = nrow(rf_input_1)
#     gini.cv$vaginal_prev[i] = sum(rf_input_1[,taxon]>0)/nrow(rf_input_1)
#     gini.cv$vaginal_prev_strict[i] = sum(rf_input_1[,taxon]>=0.001)/nrow(rf_input_1)
#     gini.cv$vaginal_mean[i] = mean(rf_input_1[,taxon])
#     gini.cv$vaginal_sd[i] = sd(rf_input_1[,taxon])
#     
#     gini.cv$caesarean_mode[i] = "caesarean"
#     gini.cv$caesarean_present[i] = sum(rf_input_2[,taxon]>0)
#     gini.cv$caesarean_present_strict[i] = sum(rf_input_2[,taxon]>=0.001)
#     gini.cv$caesarean_n[i] = nrow(rf_input_2)
#     gini.cv$caesarean_prev[i] = sum(rf_input_2[,taxon]>0)/nrow(rf_input_2)
#     gini.cv$caesarean_prev_strict[i] = sum(rf_input_2[,taxon]>=0.001)/nrow(rf_input_2)
#     gini.cv$caesarean_mean[i] = mean(rf_input_2[,taxon])
#     gini.cv$caesarean_sd[i] = sd(rf_input_2[,taxon])
#     gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>0), sum(rf_input_1[,taxon]==0),
#                                                sum(rf_input_2[,taxon]>0), sum(rf_input_2[,taxon]==0)),nrow=2))$p.value
#     gini.cv$fisher_p_strict[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>=0.001), sum(rf_input_1[,taxon]<0.001),
#                                                       sum(rf_input_2[,taxon]>=0.001), sum(rf_input_2[,taxon]<0.001)),nrow=2))$p.value
#   }
#   
#   gini.cv$fisher_padj = NA
#   gini.cv$fisher_padj[!(gini.cv$vaginal_prev<0.05 & gini.cv$caesarean_prev<0.05) & !(gini.cv$vaginal_prev>0.95 & gini.cv$caesarean_prev>0.95)] = 
#     p.adjust(gini.cv$fisher_p[!(gini.cv$vaginal_prev<0.05 & gini.cv$caesarean_prev<0.05) & !(gini.cv$vaginal_prev>0.95 & gini.cv$caesarean_prev>0.95)], method="fdr")
#   gini.cv$fisher_signif = as.numeric(gini.cv$fisher_padj<=0.1)
#   
#   gini.cv$fisher_padj_strict = NA
#   gini.cv$fisher_padj_strict[!(gini.cv$vaginal_prev_strict<0.05 & gini.cv$caesarean_prev_strict<0.05) & !(gini.cv$vaginal_prev_strict>0.95 & gini.cv$caesarean_prev_strict>0.95)] = 
#     p.adjust(gini.cv$fisher_p_strict[!(gini.cv$vaginal_prev_strict<0.05 & gini.cv$caesarean_prev_strict<0.05) & !(gini.cv$vaginal_prev_strict>0.95 & gini.cv$caesarean_prev_strict>0.95)], method="fdr")
#   gini.cv$fisher_signif_strict = as.numeric(gini.cv$fisher_padj_strict<=0.1)
#   
#   gini.cv$prev_diff = gini.cv$vaginal_prev-gini.cv$caesarean_prev
#   gini.cv$prev_diff_strict = gini.cv$vaginal_prev_strict-gini.cv$caesarean_prev_strict
#   
#   gini.cv$enriched_vaginal = gini.cv$prev_diff>0
#   gini.cv$enriched_vaginal_strict = gini.cv$prev_diff_strict>0
#   
#   gini.cv$rf_rank = 1:nrow(gini.cv)
#   gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
#   gini.cv$rf_n = n_per_group
#   
#   # add plot colour scheme
#   gini.cv$fisher_col = "ns"
#   gini.cv$fisher_col[gini.cv$prev_diff>0 & gini.cv$fisher_padj<0.1] = "vaginal"
#   gini.cv$fisher_col[gini.cv$prev_diff<0 & gini.cv$fisher_padj<0.1] = "caesarean"
#   
#   gini.cv$fisher_col_strict = "ns"
#   gini.cv$fisher_col_strict[gini.cv$prev_diff_strict>0 & gini.cv$fisher_padj_strict<0.1] = "vaginal"
#   gini.cv$fisher_col_strict[gini.cv$prev_diff_strict<0 & gini.cv$fisher_padj_strict<0.1] = "caesarean"
#   
#   # add cross-val mean + sd
#   gini.cv$crossval_mean = mean(stat.cv$accuracy)
#   gini.cv$crossval_sd = sd(stat.cv$accuracy)
#   
#   # save RF results
#   write.csv(stat.cv,file=paste0("outputs/RF_crossval/delivery/crossval_stats_deliverymode_",subset,"_",level,".csv"))
#   write.csv(gini.cv, file=paste0("outputs/RF_crossval/delivery/crossval_gini_deliverymode_",subset,"_",level,".csv"))
#   return(list(gini.cv=gini.cv))
# } 
# 
# ### function to create rfcv summary plot
# rfcv_plot = function(input_df) {
#   ggplot(input_df, aes(x=nvar, y=mean, color=level)) + geom_point(size=3, alpha=0.8) +
#     geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0, alpha=0.5) + # facet_grid(level~.) +
#     scale_x_continuous(trans='log10', breaks=subset(input_df, level=="rsv")$nvar[1:round(length(input_df$nvar)/2,0)*2-1]) + ylim(0.5,1) +
#     geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") + scale_color_manual(values=c("#034e7b", "#74a9cf")) +
#     ylab("cross-validation accuracy (%)") + xlab("No. variables") + theme( legend.position="top")
# }
# 
# 


###################################
### seroconversion RF functions ###
###################################

### function to create gini summary plot
# gini_plot_outcome = function(input_df) {
#   input_df = input_df[order(input_df$rf_rank),]
#   input_top20 = input_df[1:20,]
#   input_top20$taxonomy <- factor(input_top20$taxonomy, levels = input_top20$taxonomy)
#   ggplot(input_top20, aes(x = taxonomy, y = mean, fill=enriched_responder)) +
#     scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
#     xlab("") + ylab("importance (Gini)") +
#     geom_bar(position="dodge", stat = "identity", alpha=0.8) + coord_flip() + theme_bw() +
#     theme(legend.position = "none") + scale_fill_manual(values=c("#d95f02","#045a8d"))
# }
# 
# ### function to create prevalence/abundance plot
# rf_prev_plot_outcome = function(input_df=gini.cv) {
#   input_df = input_df[order(input_df$rf_rank),]
#   input_top20 = input_df[1:20,]
#   input_top20$taxonomy <- factor(input_top20$taxonomy, levels = input_top20$taxonomy)
#   ggplot(input_top20, aes(x = taxonomy, y = nonresponder_prev)) +
#     scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
#     xlab("") + ylab("prevalence") +
#     geom_point(aes(colour = "0", size = nonresponder_mean), stat = "identity", alpha=0.7) +
#     geom_point(aes(y = responder_prev, colour = "1", size = responder_mean), stat = "identity", alpha=0.7) +
#     coord_flip() + theme_bw() + ylim(0,1) +
#     theme(axis.text.y = element_blank(), legend.position = "right") + labs(colour = "outcome", size = "mean \nabundance") +
#     scale_colour_manual(values=c("#d95f02","#045a8d"))
# }
# 
# # fisher volcano plot
# fisher_volcano_outcome = function(df) {
#   ggplot(df, aes(x = prev_diff, y = -log10(fisher_padj))) + 
#     xlab("difference in prevalence") + ylab("log10 adjusted p") +
#     geom_point(aes(colour=factor(df$fisher_col)), size=1, alpha=0.5) + theme_bw() + xlim(-0.2,0.2) +
#     theme(legend.position = "none") + scale_colour_manual(values=c(nonresponder = "#d95f02", responder = "#045a8d", ns = "grey")) + ggtitle("Fisher's test")
# }

