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





###########################
### country RF function ###
###########################

run_rf = function(input_rps, input_ps, prev_threshold=0.05, folds=5, iterations=iter, n_per_group=50, output_dir="output_module3") {

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
  rf_reference_rps = prune_taxa(taxlist, input_rps)
  rf_input = data.frame(t(otu_table(rf_reference_rps))) 

  # append taxonomy assignment
  if (level == "genus") { names(rf_input) = tax_table(rf_reference_rps)[,"Genus"] }
  if (level == "rsv") { names(rf_input) = tax_table(rf_reference_rps)[,"tax_id"] }
  
  # create country-specific subsets
  rf_input$country_arm = factor(sample_data(rf_reference_rps)$country_arm)
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
    gini.cv = merge(gini.cv, as(tax_table(rf_reference_rps), "matrix")[,c("tax_id", "taxonomy")], by = "tax_id", all = TRUE) }
  gini.cv <- gini.cv[order(-gini.cv$mean),]
  gini.cv$taxonomy <- factor(gini.cv$taxonomy, levels = gini.cv$taxonomy)
  
  # append prevalence, abundance data and Fisher's test p values for presence/absence
  for (i in 1:nrow(gini.cv)) {
    taxon = gini.cv$tax_id[i]
    gini.cv$group1_country[i] = as.character(country1)
    gini.cv$group1_present[i] = sum(rf_input_1[,taxon]>0)
    gini.cv$group1_n[i] = nrow(rf_input_1)
    gini.cv$group1_prev[i] = round(sum(rf_input_1[,taxon]>0)/nrow(rf_input_1)*100,1)
    gini.cv$group1_mean[i] = round(mean(rf_input_1[,taxon])*100,1)
    gini.cv$group1_sd[i] = round(sd(rf_input_1[,taxon])*100,1)
    
    gini.cv$group2_country[i] = as.character(country2)
    gini.cv$group2_present[i] = sum(rf_input_2[,taxon]>0)
    gini.cv$group2_n[i] = nrow(rf_input_2)
    gini.cv$group2_prev[i] = round(sum(rf_input_2[,taxon]>0)/nrow(rf_input_2)*100,1)
    gini.cv$group2_mean[i] = round(mean(rf_input_2[,taxon])*100,1)
    gini.cv$group2_sd[i] = round(sd(rf_input_2[,taxon])*100,1)
    gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>0), sum(rf_input_1[,taxon]==0),
                                               sum(rf_input_2[,taxon]>0), sum(rf_input_2[,taxon]==0)),nrow=2))$p.value
  }
  
  # calculate fdr-adjusted p values - presence or absence
  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj = p.adjust(gini.cv$fisher_p, method="BH")
  # gini.cv$fisher_padj[!(gini.cv$group1_prev<5 & gini.cv$group2_prev<5) & !(gini.cv$group1_prev>95 & gini.cv$group2_prev>95)] = 
  #   p.adjust(gini.cv$fisher_p[!(gini.cv$group1_prev<5 & gini.cv$group2_prev<5) & !(gini.cv$group1_prev>95 & gini.cv$group2_prev>95)], method="BH")
  # gini.cv$fisher_signif = as.numeric(gini.cv$fisher_padj<=0.05)
  
  # calculate prevalence differences
  gini.cv$prev_diff = gini.cv$group1_prev-gini.cv$group2_prev
  gini.cv$enriched_group1 = gini.cv$prev_diff>0
  
  # calculate RF rank and top 20 features by importance score
  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  
  # add plot colour scheme
  gini.cv$fisher_col = "ns"
  gini.cv$fisher_col[gini.cv$prev_diff>0 & gini.cv$fisher_padj<0.05] = country1
  gini.cv$fisher_col[gini.cv$prev_diff<0 & gini.cv$fisher_padj<0.05] = country2
  
  # aldex
  aldex_ps = prune_taxa(taxlist, input_ps)
  aldex_rps = transform_sample_counts(aldex_ps, function(x) {x/sum(x)})
  aldex_rps_country1 = prune_samples(ps1list, aldex_rps)
  aldex_rps_country2 = prune_samples(ps2list, aldex_rps)
  features = data.frame(otu_table(aldex_ps))
  rfeatures = data.frame(otu_table(aldex_rps))
  groups = sample_data(aldex_ps)$country_arm
  x.all <- aldex(features, groups, test="t", effect=TRUE, include.sample.summary=FALSE, denom="all", verbose=FALSE)
  
  for (i in 1:nrow(features)) { 
    x.all$mean_rabund_country1[i] = round(mean(otu_table(aldex_rps_country1)[i,]*100),3)
    x.all$mean_rabund_country2[i] =  round(mean(otu_table(aldex_rps_country2)[i,]*100),3)
    x.all$mean_rabund_diff[i] =  round(mean(otu_table(aldex_rps_country1)[i,]*100),3) - round(mean(otu_table(aldex_rps_country2)[i,]*100),3)
  }
  if (level == "genus") { x.all$tax_id = data.frame(as(tax_table(aldex_ps), "matrix"))$Genus }
  if (level == "rsv") { x.all$tax_id = data.frame(as(tax_table(aldex_ps), "matrix"))$tax_id }
  gini.cv = merge(gini.cv, x.all, by = "tax_id", all=TRUE)
  
  # convert to DESeq2 object
  #ds = phyloseq_to_deseq2(aldex_ps, ~ country_arm)
  
  # calculate geometric means prior to estimate size factors
  # gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)) }
  # geoMeans = apply(counts(ds), 1, gm_mean)
  # ds = estimateSizeFactors(ds, geoMeans = geoMeans)
  # ds = DESeq(ds, fitType="local", parallel = TRUE)
  
  # extract results
  # res = as.data.frame(results(ds))
  
  # append rsv or genus assignment
  # if (level == "genus") { 
  #   res = cbind(as(res, "data.frame"), as(tax_table(aldex_ps)[rownames(res), c("Genus")], "matrix")) 
  #   res$tax_id = res$Genus
  # }
  # if (level == "rsv")  { res = cbind(as(res, "data.frame"), as(tax_table(aldex_ps)[rownames(res), c("tax_id", "taxonomy_species")], "matrix")) }
  # 
  # 
  # # extract DESeq2 results for taxa above this threshold, then recalculate adjusted p values
  # gini.cv = merge(gini.cv, res, by = "tax_id", all=TRUE)
  
  # add cross-val median + IQR
  gini.cv$crossval_median = median(stat.cv$accuracy)
  gini.cv$crossval_quartile_1 = quantile(stat.cv$accuracy)[2]
  gini.cv$crossval_quartile_3 = quantile(stat.cv$accuracy)[4]
  
  # tax = data.frame(as(tax_table(input_rps),"matrix"))
  # gini.cv = merge(gini.cv, tax[,c("tax_id", "taxonomy_species")], by = "tax_id")
  
  # save RF results
  write.csv(stat.cv,file=paste0(output_dir,"/RF_output_stats/crossval_stats_",country1,"_",country2,"_",subset,"_",level,".csv"))
  write.csv(gini.cv, file=paste0(output_dir,"/RF_output_importance/crossval_gini_",country1,"_",country2,"_",subset,"_",level,".csv"))
  
  # add sample group
  gini.cv$sample = "week of life 1"
  if (subset=="BS2") { gini.cv$sample = "week of life 4"}
  if (subset=="BS3") { gini.cv$sample = "week of life 6"}
  if (subset=="BS5") { gini.cv$sample = "week of life 10"}
  if (subset=="MS1") { gini.cv$sample = "mother"}
  if (subset=="BM2") { gini.cv$sample = "week of life 7"}
  if (subset=="BM3") { gini.cv$sample = "week of life 11"}
  
  # select reduced variable list for supplementary table
  gini.cv.r = gini.cv[,c("taxonomy","sample","group1_n","group1_prev","mean_rabund_country1","group2_n","group2_prev","mean_rabund_country2",
                         "fisher_p", "fisher_padj", "wi.ep", "wi.eBH","rf_rank")]
  if (any(gini.cv.r$fisher_padj<0.05 | gini.cv.r$wi.eBH<0.05)) { 
    gini.cv.r = subset(gini.cv.r, fisher_padj<0.05 | wi.eBH<0.05)
    gini.cv.r$enriched = as.character(country1)
    gini.cv.r$enriched[gini.cv.r$mean_rabund_country1<gini.cv.r$mean_rabund_country2] = as.character(country2)
    gini.cv.r$pmin = with(gini.cv.r, pmin(fisher_padj, wi.eBH))
    gini.cv.r = gini.cv.r %>% arrange(pmin) %>% dplyr::select(-pmin)
    write.csv(gini.cv.r, file=paste0(output_dir,"/RF_output_importance_reduced/crossval_gini_",country1,"_",country2,"_",subset,"_",level,".csv"))
  }
  return(list(gini.cv=gini.cv))
} 





#######################################
### ORV binary outcome RF functions ###
#######################################

### function to run collection of Random Forests
collate_rf_binary = function(output_name, output_dir="outputs", filtered_flag) {
  # create dataframe for collated accuracy results
  collated_rf = data.frame(sample_subset = rep(NA,length(subset_list)*length(country_list)))
  collated_rf$responders = collated_rf$nonresponders = collated_rf$country = NA
  
  # run RF in loop across sample types and country-country comparisons
  for (s in 1:length(subset_list)) {
    for (c in 1:length(country_list)) {
      input_m$outcome = input_m[,outcome]
      sample_subset = subset_list[s]
      country_subset = country_list[c]
      
      if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
        samplelist = subset(input_m, country_exposure==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID
      } else {
        samplelist = subset(input_m, country==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID 
      }
      
      input_rps = prune_samples(samplelist, rps)
      sample_data(input_rps)$outcome = as(sample_data(input_rps), "matrix")[,outcome]
      rf = run_rf_binary_outcome(input_rps, country_arm=country_subset, sample_sub=sample_subset, filtered=filtered_flag, prev_threshold=0.05, folds=5, 
                                 iterations=20, n_per_group=50)
      
      collated_rf$sample_subset[length(country_list)*(s-1)+c] = sample_subset
      collated_rf$country[length(country_list)*(s-1)+c] = country_subset
      collated_rf$nonresponders[length(country_list)*(s-1)+c] = as.character("nonresponders")
      collated_rf$responders[length(country_list)*(s-1)+c] = as.character("responders")
      collated_rf$crossval_median[length(country_list)*(s-1)+c] = rf$gini.cv$crossval_median[1]
      collated_rf$crossval_q1[length(country_list)*(s-1)+c] = rf$gini.cv$crossval_quartile_1[1]
      collated_rf$crossval_q3[length(country_list)*(s-1)+c] = rf$gini.cv$crossval_quartile_3[1]
    }
  }
  write.csv(collated_rf,file=paste0("outputs/",output_name))
}
### function to run Random Forests
run_rf_binary_outcome = function(input_rps, country_arm, sample_sub, prev_threshold=0.05, folds=5, filtered, iterations=iter, n_per_group=50) {
  
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

  # append species assignment
  #names(rf_input) = tax_table(rf_reference_rps)[,"Species"]
  
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
  gini.cv <- gini.cv[order(-gini.cv$mean),]

  # append prevalence, abundance data and Fisher's test p values for presence/absence
  for (i in 1:nrow(gini.cv)) {
    feature = gini.cv$feature_id[i]
    gini.cv$nonresponder_present[i] = sum(rf_input_1[,feature]>0)
    gini.cv$nonresponder_n[i] = nrow(rf_input_1)
    gini.cv$nonresponder_prev[i] = round(sum(rf_input_1[,feature]>0)/nrow(rf_input_1)*100,1)
    gini.cv$nonresponder_mean[i] = round(mean(rf_input_1[,feature])*100,1)
    gini.cv$nonresponder_sd[i] = round(sd(rf_input_1[,feature])*100,1)
    
    gini.cv$responder_present[i] = sum(rf_input_2[,feature]>0)
    gini.cv$responder_n[i] = nrow(rf_input_2)
    gini.cv$responder_prev[i] = round(sum(rf_input_2[,feature]>0)/nrow(rf_input_2)*100,1)
    gini.cv$responder_mean[i] = round(mean(rf_input_2[,feature])*100,1)
    gini.cv$responder_sd[i] = round(sd(rf_input_2[,feature])*100,1)
    gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,feature]>0), sum(rf_input_1[,feature]==0),
                                               sum(rf_input_2[,feature]>0), sum(rf_input_2[,feature]==0)),nrow=2))$p.value
    gini.cv$wilcox_p[i] = wilcox.test(rf_input[,feature] ~ rf_input$outcome)$p.value
    
  }
  
  # calculate fdr-adjusted p values - presence or absence
  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj = p.adjust(gini.cv$fisher_p, method="BH")
  gini.cv$wilcox_padj = NA
  gini.cv$wilcox_padj = p.adjust(gini.cv$wilcox_p, method="BH")
  
  # calculate prevalence differences
  gini.cv$mean_diff = gini.cv$nonresponder_mean-gini.cv$responder_mean
  gini.cv$prev_diff = gini.cv$nonresponder_prev-gini.cv$responder_prev
  gini.cv$enriched_nonresponder = gini.cv$mean_diff>0

  # calculate RF rank and top 20 features by importance score
  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  gini.cv$rf_n = n_per_group
  
  # add cross-val median + IQR
  gini.cv$crossval_median = median(stat.cv$accuracy)
  gini.cv$crossval_quartile_1 = quantile(stat.cv$accuracy)[2]
  gini.cv$crossval_quartile_3 = quantile(stat.cv$accuracy)[4]
  
  features = data.frame(tax_table(rf_reference_rps))
  features$feature_id = rownames(features)
  gini.cv = merge(gini.cv, features, by = "feature_id")
  
  # save RF results
  write.csv(stat.cv,file=paste0("outputs/RF_outputs/crossval_stats_",outcome,"_",country_arm,"_",sample_sub,"_",filtered,".csv"))
  write.csv(gini.cv, file=paste0("outputs/RF_outputs/crossval_gini_",outcome,"_",country_arm,"_",sample_sub,"_",filtered,".csv"))
  return(list(gini.cv=gini.cv))
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
  rf.fit = randomForest(x = Xtrain, y = Ytrain, ntree=500, importance=TRUE)
  
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
collate_rf_IgA = function(output_name, filtered_flag) {
  # create dataframe for collated accuracy results
  collated_rf = data.frame(sample_subset = rep(NA,length(subset_list)*length(country_list)))
  collated_rf$country = collated_rf$q1_tst_lm_R2 = collated_rf$q3_tst_lm_R2  = collated_rf$median_tst_lm_R2 = NA
  
  # run RF in loop across sample types and country-country comparisons
  for (s in 1:length(subset_list)) {
    for (c in 1:length(country_list)) {
      input_m$outcome = log(input_m[,"BB2_IgA"])
      sample_subset = subset_list[s]
      country_subset = country_list[c]
      
      if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
        samplelist = subset(input_m, country_exposure==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID
      } else {
        samplelist = subset(input_m, country==country_subset & sample_type==sample_subset & !is.na(outcome))$sample_ID
      }
      
      input_rps = prune_samples(samplelist, rps)
      input_rps = filter_taxa(input_rps, function(x) { sum(x > 0) > 0.05*nsamples(input_rps) }, TRUE)
      rf_input = data.frame(t(otu_table(input_rps))) # Create matrix of OTU abundances 

      # run RF using updated crossval_regression function
      log_IgA = log(as.numeric(as(sample_data(input_rps), "matrix")[,"BB2_IgA"]))
      K = 5 # number of folds
      B = 20 # number of repititions
      cv.out = crossval_regression(predfun = predfun_regression, X = rf_input, Y = log_IgA, K=K, B=B, verbose = FALSE)
      
      # collate cross-validation statistics - stats
      stats = data.frame(cv.out$stat.cv[,1:(ncol(cv.out$stat.cv)-1)])
      stats_df = data.frame(tr_lm_R2 = unlist(stats$tr_lm_R2), tr_spearman = unlist(stats$tr_lm_R2), tr_n = unlist(stats$tr_n),
                            tst_lm_R2 = unlist(stats$tst_lm_R2), tst_spearman = unlist(stats$tst_spearman), tst_n = unlist(stats$tst_n))
      write.csv(stats_df,file=paste0("outputs/RF_outputs/crossval_stats_IgA_",country_subset,"_",sample_subset,"_",filtered_flag,".csv"))
      
      # calculate importance score mean and SD across 1000 iterations
      IncMSE_output <- data.frame(cv.out$IncMSE_importance_output)
      for (i in 1:(length(rf_input))) {
        IncMSE_output$IncMSE_means[i] <- mean(as.numeric(IncMSE_output[i,1:B*K]))
        IncMSE_output$IncMSE_sds[i] <- sd(as.numeric(IncMSE_output[i,1:B*K]))
      }
      
      # append prevalence, abundance data and Fisher's test p values for presence/absence
      for (i in 1:nrow(IncMSE_output)) {
        feature = rownames(IncMSE_output)[i]
        IncMSE_output$present[i] = sum(rf_input[,feature]>0)
        IncMSE_output$n[i] = nrow(rf_input)
        IncMSE_output$prev[i] = round(sum(rf_input[,feature]>0)/nrow(rf_input)*100,1)
        IncMSE_output$mean[i] = round(mean(rf_input[,feature])*100,1)
        IncMSE_output$rho[i] = cor.test(log_IgA, rf_input[,feature], method="spearman")$estimate
        IncMSE_output$spearman_p[i] = cor.test(log_IgA, rf_input[,feature], method="spearman")$p.value
      }
      
      # calculate fdr-adjusted p values - presence or absence
      IncMSE_output$spearman_p_adj = p.adjust(IncMSE_output$spearman_p, method="BH")
      IncMSE_output$feature_id = rownames(IncMSE_output)
      
      # add pathway details
      features = data.frame(tax_table(input_rps))
      features$feature_id = rownames(features)
      IncMSE_output = merge(IncMSE_output, features, by = "feature_id")
      IncMSE_output_sorted <- IncMSE_output[order(-IncMSE_output$IncMSE_means),] 
      write.csv(IncMSE_output_sorted,file=paste0("outputs/RF_outputs/crossval_regression_IgA_",country_subset,"_",sample_subset,"_",filtered_flag,".csv"))
      
      # add RF statistics to collated_rf dataframe
      collated_rf$sample_subset[length(country_list)*(s-1)+c] = sample_subset
      collated_rf$country[length(country_list)*(s-1)+c] = country_subset
      collated_rf$median_tst_lm_R2[length(country_list)*(s-1)+c] = median(unlist(stats$tst_lm_R2))
      collated_rf$q1_tst_lm_R2[length(country_list)*(s-1)+c] = quantile(unlist(stats$tst_lm_R2))[2]
      collated_rf$q3_tst_lm_R2[length(country_list)*(s-1)+c] = quantile(unlist(stats$tst_lm_R2))[4]
    }
  }
  write.csv(collated_rf,file=paste0("outputs/",output_name))
}






# function to create gini summary plot
gini_plot = function(input_df, colour1=col1, colour2=col2, delivery=FALSE, species=TRUE) {
  if (delivery==TRUE) {
    input_df$enriched_group1 = input_df$enriched_vaginal
  }
  
  input_df = input_df[order(input_df$rf_rank),]
  input_top20 = input_df[1:10,]
  
  if (species==TRUE) {
    input_top20$taxonomy = input_top20$taxonomy_species
  }
  
  input_top20$taxonomy = gsub("_", " ", input_top20$taxonomy) 
  input_top20$taxonomy <- factor(input_top20$taxonomy, levels = input_top20$taxonomy)
  
  if (all(input_top20$enriched_group1==TRUE)) {
    ggplot(input_top20, aes(x = taxonomy, y = mean, fill=enriched_group1)) +
      scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
      xlab("") + ylab("importance (Gini)") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + scale_fill_manual(values=c(colour1)) +
      theme(axis.text = element_text(size=22), axis.title = element_text(size=30), plot.title = element_text(size=25), strip.text = element_text(size=20),
            axis.text.x = element_text(size=25))
  } else {
    ggplot(input_top20, aes(x = taxonomy, y = mean, fill=enriched_group1)) +
      scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
      xlab("") + ylab("importance (Gini)") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + scale_fill_manual(values=c(colour2,colour1)) +
      theme(axis.text.y = element_text(size=22), axis.title = element_text(size=30), plot.title = element_text(size=25), strip.text = element_text(size=20), 
            axis.text.x = element_text(size=25))
  }
}

# function to create prevalence/abundance plot
rf_prev_plot = function(input_df=gini.cv, label1=country1, colour1=col1, label2=country2, colour2=col2, species=TRUE) {
  input_df = input_df[order(input_df$rf_rank),]
  input_top20 = input_df[1:10,]
  
  if (species==TRUE) {
    input_top20$taxonomy = input_top20$taxonomy_species
  }
  input_top20$taxonomy = factor(input_top20$taxonomy, levels = input_top20$taxonomy)
  
  ggplot(input_top20, aes(x = taxonomy, y = group1_prev)) +
    scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
    xlab("") + ylab("prevalence (%)") +
    geom_point(aes(colour = label1, size=group1_mean/100), stat = "identity", alpha=0.9) +
    geom_point(aes(y = group2_prev, colour = label2, size=group2_mean/100), stat = "identity", alpha=0.9) +
    coord_flip() + theme_bw() + scale_y_continuous(limits = c(0,105), breaks=c(0,50,100)) +
    theme(axis.text.y = element_blank(), legend.position = "right") + labs(colour = "enriched group", size = "mean \nabundance") +
    scale_colour_manual(values=c(colour1,colour2)) + scale_size(limits = c(0, 0.5),  range=c(3,15)) +
    guides(colour = guide_legend(order=1), size = guide_legend(order=2)) +
    theme(axis.text = element_text(size=20), axis.title = element_text(size=30), plot.title = element_text(size=25), strip.text = element_text(size=20),
          legend.text = element_text(size=20), legend.title = element_text(size=20), axis.text.x = element_text(size=25))
  #guides(colour = "none")
  
}


