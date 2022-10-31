########################################
### Functions for RoVI analyses work ###
########################################
library(VennDiagram)

### Module 6

mother_infant_prev = function(country_select, prev_threshold=0.1, maternal="MS1") {
  # create separate dataframes for maternal and infant samples 
  t_M = subset(t, country==country_select & sample_type==maternal & mode_delivery=="vaginal")
  t_B = subset(t, country==country_select & sample_type=="BS1" & mode_delivery=="vaginal")
  
  # select families present in both datasets, then pick out corresponding subsets
  family_list = unique(t_M$family_ID)[unique(t_M$family_ID) %in% unique(t_B$family_ID)]
  sample_list_M = t_M$sample_ID_full[t_M$family_ID %in% family_list]
  sample_list_B = t_B$sample_ID_full[t_B$family_ID %in% family_list]
  
  # prune rps8 to relevant subsets
  rps8_M = prune_samples(sample_list_M, rps8) %>% prune_taxa(taxa_sums(.) > 0, .)
  rps8_B = prune_samples(sample_list_B, rps8) 
  
  ### maternal samples ###
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps8_M_no0 = rps8_M
  otu_table(rps8_M_no0)[otu_table(rps8_M_no0)==0] = NA
  
  # calculate taxon means
  tax_abund = data.frame(round(rowMeans(otu_table(rps8_M_no0), na.rm=T),6), 
                         round(rowSds(otu_table(rps8_M_no0), na.rm=T),6), rowSums(otu_table(rps8_M)>0))
  names(tax_abund) = c("mean_abundance_M", "SD_abundance_M", "present_M")
  tax_abund$prevalence_M = round(tax_abund$present_M/ncol(otu_table(rps8_M)),3)
  
  # append taxon names
  tax_abund$taxon = data.frame(as(tax_table(rps8_M),"matrix"))$taxonomy_species
  rownames(tax_abund) = tax_abund$taxon
  
  # order by decreasing abundance and filter to 20% prevalence
  tax_abund = tax_abund[order(-tax_abund$prevalence_M, -tax_abund$mean_abundance_M),]
  tax_abund = subset(tax_abund, prevalence_M>=prev_threshold)
  
  # create dataframe of feature table for corresponding taxa
  otu_M = data.frame(t(otu_table(rps8_M))) 
  if (all(colnames(otu_M)==rownames(tax_table(rps8_M)))) { 
    colnames(otu_M) = as(tax_table(rps8_M),"matrix")[,"taxonomy_species"]} else { print("error: rownames do not match")}
  otu_M = otu_M[,rownames(tax_abund)]
  otu_M$family_ID = substr(substrRight(rownames(otu_M),7),1,4)
  otu_M = otu_M[order(otu_M$family_ID),]
  
  ### infant samples ###
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps8_B_no0 = rps8_B
  otu_table(rps8_B_no0)[otu_table(rps8_B_no0)==0] = NA
  
  # calculate taxon means
  tax_abund_B = data.frame(round(rowMeans(otu_table(rps8_B_no0), na.rm=T),6),
                           round(rowSds(otu_table(rps8_B_no0), na.rm=T),6), rowSums(otu_table(rps8_B)>0))
  names(tax_abund_B) = c("mean_abundance_BS1", "SD_abundance_BS1", "present_BS1")
  tax_abund_B$prevalence_BS1 = round(tax_abund_B$present_BS1/ncol(otu_table(rps8_B)),3)
  
  # append taxon names
  tax_abund_B$taxon = data.frame(as(tax_table(rps8_B),"matrix"))$taxonomy_species
  rownames(tax_abund_B) = tax_abund_B$taxon
  tax_abund_B = tax_abund_B[rownames(tax_abund),]
  
  # create dataframe of feature table for corresponding taxa
  otu_B = data.frame(t(otu_table(rps8_B))) 
  if (all(colnames(otu_B)==rownames(tax_table(rps8_B)))) { 
    colnames(otu_B) = as(tax_table(rps8_B),"matrix")[,"taxonomy_species"] 
  } else { print("error: rownames do not match")}
  otu_B = otu_B[,rownames(tax_abund_B)]
  otu_B$family_ID = substr(substrRight(rownames(otu_B),7),1,4)
  otu_B = otu_B[order(otu_B$family_ID),]
  
  # collate results
  collated = data.frame(cbind(tax_abund, tax_abund_B[,1:4]))
  
  # add wilcox tests
  for (i in 1:nrow(collated)) {
    if (sum(as.numeric(otu_B[,i]>0)) >= 10 & sum(as.numeric(otu_M[,i]>0)) >= 10) {
      collated$wilcox_p[i] = wilcox.test(otu_M[,i]~as.numeric(otu_B[,i]>0))$p.value
      collated$spearman_p[i] = cor.test(otu_M[,i],otu_B[,i], method="spearman")$p.value
      collated$spearman_rho[i] = cor.test(otu_M[,i],otu_B[,i], method="spearman")$estimate
      
    } else { collated$wilcox_p[i] = collated$spearman_p[i] = collated$spearman_rho[i] = NA }
  }
  
  collated$fdr_wilcox_p = NA
  collated$fdr_wilcox_p[!is.na(collated$wilcox_p)] = p.adjust(collated$wilcox_p[!is.na(collated$wilcox_p)], method="fdr")
  
  collated$mean_abundance_BS1[is.na(collated$mean_abundance_BS1)]=0
  collated$signif = factor(collated$fdr_wilcox_p<0.05)
  collated$signif[is.na(collated$signif)] = FALSE
  collated
}

mother_infant_plot = function(collated, country_select, ntaxa=25) {
  collated$taxon = factor(collated$taxon, levels = collated$taxon)
  g1 = ggplot(collated[1:ntaxa,], aes(x = taxon, y = prevalence_M)) +
    scale_x_discrete(limits = rev(levels(collated$taxon)[1:ntaxa])) +
    xlab("") + ylab("prevalence") +
    geom_point(aes(size = mean_abundance_M, colour = "maternal"), stat = "identity", alpha=0.5) +
    coord_flip() + theme_bw() + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
    labs(colour = "sample", size = "mean \nabundance") +
    scale_colour_manual(values=c(my_palette[5])) + scale_size(limits = c(0, 0.25)) +
    guides(colour = guide_legend(order=1), size = guide_legend(order=2)) +
    theme(legend.position = "none", plot.title = element_text(size=12)) + 
    ggtitle(paste0(country_select, " - mother"))
  g2 = ggplot(collated[1:ntaxa,], aes(x = taxon, y = prevalence_BS1, colour = signif)) +
    scale_x_discrete(limits = rev(levels(collated$taxon)[1:ntaxa])) +
    xlab("") + ylab("prevalence") + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
    geom_point(aes(size = mean_abundance_BS1), stat = "identity", alpha=0.5) +
    coord_flip() + theme_bw() + scale_size(limits = c(0, 0.25)) +
    scale_colour_manual(values=c(my_palette[1], my_palette[4])) +
    theme(axis.text.y = element_blank(), legend.position = "right", plot.title = element_text(size=12)) + 
    labs(colour = "significant", size = "mean \nabundance")  + ggtitle("infant")
  return(list(g1,g2))
}


infant_mother_prev = function(country_select, prev_threshold=0.1, delivery_mode="vaginal", maternal="MS1") {
  # create separate dataframes for maternal and infant samples 
  t_B = subset(t, country==country_select & sample_type=="BS1" & mode_delivery==delivery_mode)
  t_M = subset(t, country==country_select & sample_type==maternal & mode_delivery==delivery_mode)
  
  # select families present in both datasets, then pick out corresponding subsets
  family_list = unique(t_M$family_ID)[unique(t_M$family_ID) %in% unique(t_B$family_ID)]
  sample_list_B = t_B$sample_ID_full[t_B$family_ID %in% family_list]
  sample_list_M = t_M$sample_ID_full[t_M$family_ID %in% family_list]
  
  # prune rps8 to relevant subsets
  rps8_B = prune_samples(sample_list_B, rps8) %>% prune_taxa(taxa_sums(.) > 0, .)
  rps8_M = prune_samples(sample_list_M, rps8)
  
  ### maternal samples ###
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps8_B_no0 = rps8_B
  otu_table(rps8_B_no0)[otu_table(rps8_B_no0)==0] = NA
  
  # calculate taxon means
  tax_abund = data.frame(round(rowMeans(otu_table(rps8_B_no0), na.rm=T),6), 
                         round(rowSds(otu_table(rps8_B_no0), na.rm=T),6), rowSums(otu_table(rps8_B)>0))
  names(tax_abund) = c("mean_abundance_BS1", "SD_abundance_BS1", "present_BS1")
  tax_abund$prevalence_BS1 = round(tax_abund$present_BS1/ncol(otu_table(rps8_B)),3)
  
  # append taxon names
  tax_abund$taxon = data.frame(as(tax_table(rps8_B),"matrix"))$taxonomy_species
  rownames(tax_abund) = tax_abund$taxon
  
  # order by decreasing abundance and filter to 20% prevalence
  tax_abund = tax_abund[order(-tax_abund$prevalence_BS1, -tax_abund$mean_abundance_BS1),]
  tax_abund = subset(tax_abund, prevalence_BS1>=prev_threshold)
  
  # create dataframe of feature table for corresponding taxa
  otu_B = data.frame(t(otu_table(rps8_B))) 
  if (all(colnames(otu_B)==rownames(tax_table(rps8_B)))) { 
    colnames(otu_B) = as(tax_table(rps8_B),"matrix")[,"taxonomy_species"]} else { print("error: rownames do not match")}
  otu_B = otu_B[,rownames(tax_abund)]
  otu_B$family_ID = substr(substrRight(rownames(otu_B),7),1,4)
  otu_B = otu_B[order(otu_B$family_ID),]
  
  ### infant samples ###
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps8_M_no0 = rps8_M
  otu_table(rps8_M_no0)[otu_table(rps8_M_no0)==0] = NA
  
  # calculate taxon means
  tax_abund_M = data.frame(round(rowMeans(otu_table(rps8_M_no0), na.rm=T),6),
                           round(rowSds(otu_table(rps8_M_no0), na.rm=T),6), rowSums(otu_table(rps8_M)>0))
  names(tax_abund_M) = c("mean_abundance_M", "SD_abundance_M", "present_M")
  tax_abund_M$prevalence_M = round(tax_abund_M$present_M/ncol(otu_table(rps8_M)),3)
  
  # append taxon names
  tax_abund_M$taxon = data.frame(as(tax_table(rps8_M),"matrix"))$taxonomy_species
  rownames(tax_abund_M) = tax_abund_M$taxon
  tax_abund_M = tax_abund_M[rownames(tax_abund),]
  
  # create dataframe of feature table for corresponding taxa
  otu_M = data.frame(t(otu_table(rps8_M))) 
  if (all(colnames(otu_M)==rownames(tax_table(rps8_M)))) { 
    colnames(otu_M) = as(tax_table(rps8_M),"matrix")[,"taxonomy_species"] 
  } else { print("error: rownames do not match")}
  otu_M = otu_M[,rownames(tax_abund_M)]
  otu_M$family_ID = substr(substrRight(rownames(otu_M),7),1,4)
  otu_M = otu_M[order(otu_M$family_ID),]
  
  # collate results
  collated = data.frame(cbind(tax_abund, tax_abund_M[,1:4]))
  
  # add wilcox tests
  for (i in 1:nrow(collated)) {
    if (sum(as.numeric(otu_B[,i]>0)) >= 10 & sum(as.numeric(otu_M[,i]>0)) >= 10) {
      collated$wilcox_p[i] = wilcox.test(otu_M[,i]~as.numeric(otu_B[,i]>0))$p.value
      collated$spearman_p[i] = cor.test(otu_B[,i],otu_M[,i], method="spearman")$p.value
      collated$spearman_rho[i] = cor.test(otu_B[,i],otu_M[,i], method="spearman")$estimate
    } else { collated$wilcox_p[i] = collated$spearman_p[i] = collated$spearman_rho[i] = NA }
  }
  
  collated$fdr_wilcox_p = NA
  collated$fdr_wilcox_p[!is.na(collated$wilcox_p)] = p.adjust(collated$wilcox_p[!is.na(collated$wilcox_p)], method="fdr")
  
  collated$mean_abundance_M[is.na(collated$mean_abundance_M)]=0
  collated$signif = factor(collated$fdr_wilcox_p<0.05)
  collated$signif[is.na(collated$signif)] = FALSE
  collated
}

infant_mother_plot = function(collated, country_select, ntaxa=25) {
  collated = collated[1:ntaxa,]
  collated$taxon = factor(collated$taxon, levels = collated$taxon)
  
  g1 = ggplot(collated, aes(x = taxon, y = prevalence_BS1)) +
    scale_x_discrete(limits = rev(levels(collated$taxon))) +
    xlab("") + ylab("prevalence") +
    geom_point(aes(size = mean_abundance_BS1, colour = "infant"), stat = "identity", alpha=0.7) +
    coord_flip() + theme_bw() + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
    labs(colour = "sample", size = "mean \nabundance") +
    scale_colour_manual(values=c(my_palette[1])) + scale_size(limits = c(0, 0.35)) +
    guides(colour = guide_legend(order=1), size = guide_legend(order=2)) +
    theme(legend.position = "none", plot.title = element_text(size=12)) + 
    ggtitle(paste0(country_select, " - infant"))
  g2 = ggplot(collated, aes(x = taxon, y = prevalence_M, colour = signif)) +
    scale_x_discrete(limits = rev(levels(collated$taxon))) +
    xlab("") + ylab("prevalence") + scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1)) +
    geom_point(aes(size = mean_abundance_M), stat = "identity", alpha=0.7) +
    coord_flip() + theme_bw() + scale_size(limits = c(0, 0.35)) +
    scale_colour_manual(values=c(my_palette[5], my_palette[4])) +
    theme(axis.text.y = element_blank(), legend.position = "right", plot.title = element_text(size=12)) + 
    labs(colour = "significant", size = "mean \nabundance")  + ggtitle("mother")
  return(list(g1,g2))
}

# longitudinal plot of alpha diversity
alpha_longitudinal = function(t1, alpha_measure, alpha_title, outcome, outcome_title) {
  t1$alpha_measure = t1[,alpha_measure]
  t1$outcome = t1[,outcome]
  ggplot(t1, aes(x = factor(week), y = alpha_measure, group=family_ID, colour=factor(outcome))) + 
    geom_line(alpha=0.04) +
    scale_color_manual(outcome_title, values = c(cs_col, vd_col)) + 
    stat_summary(aes(group = outcome), geom = "line", fun.y = mean, size = 1, alpha = 0.8) +# ylim(0,100) +
    stat_summary(aes(group = outcome), geom = "point", fun.data = "mean_cl_boot", size = 2, alpha = 0.8) + ggtitle(" ") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(breaks=c("1","4","6","10"), labels=c("1", "4", "6", "10")) + xlab("week") + ylab(alpha_title) + theme(legend.position = "top")
}

# Statistics
lme_alpha = function(df, alpha_measure, outcome) {
  lm_base <- lmer(formula(paste0(alpha_measure," ~ (1 | family_ID)")), df, REML=TRUE)
  lm_week <- lmer(formula(paste0(alpha_measure," ~  week + (1|family_ID)")), df, REML=TRUE)
  lm_outcome <- lmer(formula(paste0(alpha_measure," ~  week + ", outcome, " + (1|family_ID)")), df, REML=TRUE)
  round(anova(lm_outcome, lm_week)$Pr[2],4)
}

cross_sectional_alpha = function(df, sample_subset, alpha_measure, alpha_label) {
  df$alpha = df[,alpha_measure]
  ggplot(subset(df, sample_type==sample_subset), aes(factor(mode_delivery), y=alpha, colour=factor(mode_delivery))) + 
    geom_jitter(size = 1, alpha = 0.9, width=0.2) +
    geom_boxplot(alpha = 0.5)  + xlab("") + ylab(alpha_label) +
    theme(legend.position = "none", strip.background = element_blank()) + 
    scale_colour_manual(values = c(cs_col, vd_col)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Cross-sectional comparisons
delivery_summary = function(input_df, sample_subset, alpha_variable, anova=TRUE) {
  df_sub = subset(input_df, sample_type==sample_subset) 
  df_sub$alpha_variable = df_sub[,alpha_variable]
  df_sub = subset(df_sub, !is.na(alpha_variable))
  
  cs = round(mean(df_sub$alpha_variable[df_sub$mode_delivery=="caesarean"]),2)
  sd_cs = round(sd(df_sub$alpha_variable[df_sub$mode_delivery=="caesarean"]),2)
  cs1 = paste0(cs," (",sd_cs,")")
  n_cs = sum(df_sub$mode_delivery=="caesarean")
  
  vd = round(mean(df_sub$alpha_variable[df_sub$mode_delivery=="vaginal"]),2)
  sd_vd = round(sd(df_sub$alpha_variable[df_sub$mode_delivery=="vaginal"]),2)
  vd1 = paste0(vd," (",sd_vd,")")
  n_vd = sum(df_sub$mode_delivery=="vaginal")
  
  total = paste0((n_cs+n_vd)," (",n_cs,")")
  
  if (anova==TRUE) {
    if (alpha_variable=="Observed" | alpha_variable=="Observed_rsv") {
      p = round(summary(aov(log10(df_sub$alpha_variable) ~ df_sub$mode_delivery))[[1]][["Pr(>F)"]],4)
    } else {
      p = round(summary(aov(df_sub$alpha_variable ~ df_sub$mode_delivery))[[1]][["Pr(>F)"]],4)
    }      
  } else {
    p = round(wilcox.test(df_sub$alpha_variable ~ df_sub$mode_delivery)$p.value,4)
  }
  
  c(sample_subset, alpha_variable, cs1, vd1, total, p)
}

summary_cross_sectional = function(input_df, alpha1, alpha2) {
  res = rbind(
    delivery_summary(input_df, sample_subset="MS1", alpha1),
    delivery_summary(input_df, sample_subset="BS1", alpha1),
    delivery_summary(input_df, sample_subset="BS2", alpha1),
    delivery_summary(input_df, sample_subset="BS3", alpha1),
    delivery_summary(input_df, sample_subset="BS5", alpha1),
    delivery_summary(input_df, sample_subset="MS1", alpha2),
    delivery_summary(input_df, sample_subset="BS1", alpha2),
    delivery_summary(input_df, sample_subset="BS2", alpha2),
    delivery_summary(input_df, sample_subset="BS3", alpha2),
    delivery_summary(input_df, sample_subset="BS5", alpha2))
  res=res[,1:6]
  colnames(res) = c("subset", "variable" , "caesarean", "vaginal", "n_total (caesarean)", "p")
  res 
}

summary_cross_sectional_BM = function(input_df, alpha1, alpha2) {
  res = rbind(
    delivery_summary(input_df, sample_subset="BM1", alpha1),
    delivery_summary(input_df, sample_subset="BM2", alpha1),
    delivery_summary(input_df, sample_subset="BM3", alpha1),
    delivery_summary(input_df, sample_subset="BM1", alpha2),
    delivery_summary(input_df, sample_subset="BM2", alpha2),
    delivery_summary(input_df, sample_subset="BM3", alpha2))
  res=res[,1:6]
  colnames(res) = c("subset", "variable" , "caesarean", "vaginal", "n_total (caesarean)", "p")
  res 
}

summary_prev_plot= function(level, threshold=0.1, yshift=2, ymax=85, strict=TRUE) {
  BS1 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BS1_",level,".csv"), row.names=1)
  BS2 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BS2_",level,".csv"), row.names=1)
  BS3 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BS3_",level,".csv"), row.names=1)
  BS5 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BS5_",level,".csv"), row.names=1)
  
  if (strict==TRUE) {
    BS1$fisher_padj = BS1$fisher_padj_strict
    BS2$fisher_padj = BS2$fisher_padj_strict
    BS3$fisher_padj = BS3$fisher_padj_strict
    BS5$fisher_padj = BS5$fisher_padj_strict
    
    BS1$enriched_vaginal = BS1$enriched_vaginal_strict
    BS2$enriched_vaginal = BS2$enriched_vaginal_strict
    BS3$enriched_vaginal = BS3$enriched_vaginal_strict
    BS5$enriched_vaginal = BS5$enriched_vaginal_strict
  }
  
  fisher_summary = data.frame(
    count = c(sum(BS1$enriched_vaginal==TRUE & BS1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS2$enriched_vaginal==TRUE & BS2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS3$enriched_vaginal==TRUE & BS3$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS5$enriched_vaginal==TRUE & BS5$fisher_padj<threshold, na.rm=TRUE),
              
              sum(BS1$enriched_vaginal==FALSE & BS1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS2$enriched_vaginal==FALSE & BS2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS3$enriched_vaginal==FALSE & BS3$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS5$enriched_vaginal==FALSE & BS5$fisher_padj<threshold, na.rm=TRUE)),
    
    enriched = c(rep("vaginal",4), rep("caesarean",4)),
    week = factor(rep(c("1", "4", "6", "10"),2), levels = c("1", "4", "6", "10")),
    total = c(sum(!is.na(BS1$fisher_padj)), sum(!is.na(BS2$fisher_padj)), sum(!is.na(BS3$fisher_padj)), sum(!is.na(BS5$fisher_padj)))
  )
  
  ggplot(fisher_summary, aes(x=week, y=count, label=count, colour=enriched, group=enriched)) + 
    geom_point(size=3, alpha=0.8) + geom_line(alpha=0.5, show.legend = FALSE) + 
    geom_text(vjust = 0, nudge_y = yshift, show.legend = FALSE) +
    geom_text(aes(y = rep(ymax,nrow(fisher_summary)), label=total), colour="grey", show.legend = FALSE, fontface = "italic") +
    ylim(0,ymax) + xlab("week") +
    scale_colour_manual(values=c(cs_col, vd_col)) + 
    theme(legend.position="right") +
    ylab("number of enriched taxa") + labs(colour='enriched group')
}

summary_prev_plot_BM = function(level, threshold=0.1, yshift=2, ymax=85, strict=TRUE) {
  BM1 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BM1_",level,".csv"), row.names=1)
  BM2 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BM2_",level,".csv"), row.names=1)
  BM3 = read.csv(paste0("output_module6/RF_outputs/crossval_gini_deliverymode_BM3_",level,".csv"), row.names=1)

  if (strict==TRUE) {
    BM1$fisher_padj = BM1$fisher_padj_strict
    BM2$fisher_padj = BM2$fisher_padj_strict
    BM3$fisher_padj = BM3$fisher_padj_strict

    BM1$enriched_vaginal = BM1$enriched_vaginal_strict
    BM2$enriched_vaginal = BM2$enriched_vaginal_strict
    BM3$enriched_vaginal = BM3$enriched_vaginal_strict
  }
  
  fisher_summary = data.frame(
    count = c(sum(BM1$enriched_vaginal==TRUE & BM1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BM2$enriched_vaginal==TRUE & BM2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BM3$enriched_vaginal==TRUE & BM3$fisher_padj<threshold, na.rm=TRUE), 

              sum(BM1$enriched_vaginal==FALSE & BM1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BM2$enriched_vaginal==FALSE & BM2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BM3$enriched_vaginal==FALSE & BM3$fisher_padj<threshold, na.rm=TRUE)), 

    enriched = c(rep("vaginal",3), rep("caesarean",3)),
    week = factor(rep(c("1", "7", "11"),2), levels = c("1", "7", "11")),
    total = c(sum(!is.na(BM1$fisher_padj)), sum(!is.na(BM2$fisher_padj)), sum(!is.na(BM3$fisher_padj)))
  )
  
  ggplot(fisher_summary, aes(x=week, y=count, label=count, colour=enriched, group=enriched)) + 
    geom_point(size=3, alpha=0.8) + geom_line(alpha=0.5, show.legend = FALSE) + 
    geom_text(vjust = 0, nudge_y = yshift, show.legend = FALSE) +
    geom_text(aes(y = rep(ymax,nrow(fisher_summary)), label=total), colour="grey", show.legend = FALSE, fontface = "italic") +
    ylim(0,ymax) + xlab("week") +
    scale_colour_manual(values=c(cs_col, vd_col)) + 
    theme(legend.position="right") +
    ylab("number of enriched taxa") + labs(colour='enriched group')
}

# fisher volcano plot
fisher_volcano = function(df, title, strict=TRUE) {
  if (strict==TRUE) {
    df$fisher_padj = df$fisher_padj_strict
    df$fisher_col = df$fisher_col_strict
  }
  
  ggplot(df, aes(x = prev_diff, y = -log10(fisher_padj))) + 
    xlab("difference in prevalence") + ylab("log10 adjusted p") +
    geom_point(aes(colour=factor(df$fisher_col)), size=2, alpha=0.7) + theme_bw() + 
    theme(legend.position = "none") + xlim(-0.3,0.3) +
    labs(colour = "enriched group") +
    scale_colour_manual(values=c("vaginal"=vd_col, "ns"="grey", "caesarean"=cs_col)) + ggtitle(title)
}



### function to run Random Forests (run after deseq results so that list of taxa can be transferred)
run_rf_delivery = function(input_rps, prev_threshold=0.01, folds=5, iterations=iter, n_per_group=50) {

  # determine level
  if (all(is.na(tax_table(input_rps)[,"tax_id"]))) { level = "genus"}
  if (all(!is.na(tax_table(input_rps)[,"tax_id"]))) { level = "rsv"}

  # determine taxa present above threshold in at least one group
  ps1list = sample_names(input_rps)[sample_data(input_rps)$mode_delivery=="vaginal"]
  ps1 = prune_samples(ps1list, input_rps)
  ps1 = filter_taxa(ps1, function(x) { sum(x > 0) > prev_threshold*nsamples(ps1) }, TRUE)
  ps2list = sample_names(input_rps)[sample_data(input_rps)$mode_delivery=="caesarean"]
  ps2 = prune_samples(ps2list, input_rps)
  ps2 = filter_taxa(ps2, function(x) { sum(x > 0) > prev_threshold*nsamples(ps2) }, TRUE)
  taxlist = unique(c(rownames(tax_table(ps1)),rownames(tax_table(ps2))))

  # if <50 responders or <50 non-responders, amend n_per_group to match minimum group size
  mingroup = min(c(length(ps1list), length(ps2list)))
  if (mingroup < 50) { n_per_group = mingroup }

  # create RF input
  rf_reference_ps = prune_taxa(taxlist, input_rps)
  rf_input = data.frame(t(otu_table(rf_reference_ps))) # Create matrix of OTU abundances
  #all(names(rf_input) == rownames(tax_table(rf_reference_ps)))

  # append rsv or genus assignment
  if (level == "genus") { names(rf_input) = tax_table(rf_reference_ps)[,"Genus"] }
  if (level == "rsv") { names(rf_input) = tax_table(rf_reference_ps)[,"tax_id"] }

  # all(rownames(rf_input) == rownames(sample_data(rf_reference_ps)))
  rf_input$mode_delivery = factor(sample_data(rf_reference_ps)$mode_delivery)

  # create mode-specific subsets
  rf_input_1 = subset(rf_input, mode_delivery=="vaginal")
  rf_input_2 = subset(rf_input, mode_delivery=="caesarean")

  # perform crossval on whole dataset
  crossval_stats = list()
  crossval_importance = list()

  for (i in 1:iter) {
    # random draw for each iteration
    rf_s = rbind(rf_input_1[sample(nrow(rf_input_1),n_per_group),],rf_input_2[sample(nrow(rf_input_2),n_per_group),])
    # perform crossval for iteration above
    cv.rf = crossval_update(predfun, X=rf_s[,!grepl("mode_delivery", names(rf_s))], Y=factor(rf_s$mode_delivery), K=folds, B=1, verbose=FALSE)
    # write outputs
    crossval_stats[[i]] = cv.rf$stat.cv
    crossval_importance[[i]] = cv.rf$Gini_importance_output
  }

  # collate crossvalidation statistics
  stat.cv = as.data.frame(crossval_stats[[1]])
  for (i in 2:iter) { stat.cv = rbind(stat.cv, as.data.frame(crossval_stats[[i]])) }
  stat.cv$accuracy = (stat.cv$TP+stat.cv$TN)/(stat.cv$TP+stat.cv$TN+stat.cv$FP+stat.cv$FN)

  # collate crossvalidation statistics
  gini.cv = as.data.frame(crossval_importance)
  gini.cv$mean = rowMeans(gini.cv[,1:(iter*folds)])
  gini.cv$sd = rowSds(as.matrix(gini.cv[,1:(iter*folds)]))
  gini.cv = gini.cv[,(ncol(gini.cv)-1):ncol(gini.cv)]

  if (level == "genus") { gini.cv$taxonomy = gini.cv$tax_id = rownames(gini.cv) }
  if (level == "rsv") {
    gini.cv$tax_id = rownames(gini.cv)
    gini.cv = merge(gini.cv, as(tax_table(rf_reference_ps), "matrix")[,c("tax_id", "taxonomy", "taxonomy_species")], by = "tax_id", keep.x = TRUE) }
  gini.cv <- gini.cv[order(-gini.cv$mean),]
  gini.cv$taxonomy <- factor(gini.cv$taxonomy, levels = gini.cv$taxonomy)

  # append prevalence, abundance data and Fisher p
  for (i in 1:nrow(gini.cv)) {
    taxon = gini.cv$tax_id[i]
    gini.cv$vaginal_present[i] = sum(rf_input_1[,taxon]>0)
    gini.cv$vaginal_present_strict[i] = sum(rf_input_1[,taxon]>=0.001)
    gini.cv$vaginal_n[i] = nrow(rf_input_1)
    gini.cv$vaginal_prev[i] = sum(rf_input_1[,taxon]>0)/nrow(rf_input_1)
    gini.cv$vaginal_prev_strict[i] = sum(rf_input_1[,taxon]>=0.001)/nrow(rf_input_1)
    gini.cv$vaginal_mean[i] = mean(rf_input_1[,taxon])
    gini.cv$vaginal_sd[i] = sd(rf_input_1[,taxon])

    gini.cv$caesarean_mode[i] = "caesarean"
    gini.cv$caesarean_present[i] = sum(rf_input_2[,taxon]>0)
    gini.cv$caesarean_present_strict[i] = sum(rf_input_2[,taxon]>=0.001)
    gini.cv$caesarean_n[i] = nrow(rf_input_2)
    gini.cv$caesarean_prev[i] = sum(rf_input_2[,taxon]>0)/nrow(rf_input_2)
    gini.cv$caesarean_prev_strict[i] = sum(rf_input_2[,taxon]>=0.001)/nrow(rf_input_2)
    gini.cv$caesarean_mean[i] = mean(rf_input_2[,taxon])
    gini.cv$caesarean_sd[i] = sd(rf_input_2[,taxon])
    gini.cv$fisher_p[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>0), sum(rf_input_1[,taxon]==0),
                                               sum(rf_input_2[,taxon]>0), sum(rf_input_2[,taxon]==0)),nrow=2))$p.value
    gini.cv$fisher_p_strict[i] = fisher.test(matrix(c(sum(rf_input_1[,taxon]>=0.001), sum(rf_input_1[,taxon]<0.001),
                                                      sum(rf_input_2[,taxon]>=0.001), sum(rf_input_2[,taxon]<0.001)),nrow=2))$p.value
  }

  gini.cv$fisher_padj = NA
  gini.cv$fisher_padj[!(gini.cv$vaginal_prev<0.05 & gini.cv$caesarean_prev<0.05) & !(gini.cv$vaginal_prev>0.95 & gini.cv$caesarean_prev>0.95)] =
    p.adjust(gini.cv$fisher_p[!(gini.cv$vaginal_prev<0.05 & gini.cv$caesarean_prev<0.05) & !(gini.cv$vaginal_prev>0.95 & gini.cv$caesarean_prev>0.95)], method="fdr")
  gini.cv$fisher_signif = as.numeric(gini.cv$fisher_padj<=0.1)

  gini.cv$fisher_padj_strict = NA
  gini.cv$fisher_padj_strict[!(gini.cv$vaginal_prev_strict<0.05 & gini.cv$caesarean_prev_strict<0.05) & !(gini.cv$vaginal_prev_strict>0.95 & gini.cv$caesarean_prev_strict>0.95)] =
    p.adjust(gini.cv$fisher_p_strict[!(gini.cv$vaginal_prev_strict<0.05 & gini.cv$caesarean_prev_strict<0.05) & !(gini.cv$vaginal_prev_strict>0.95 & gini.cv$caesarean_prev_strict>0.95)], method="fdr")
  gini.cv$fisher_signif_strict = as.numeric(gini.cv$fisher_padj_strict<=0.1)

  gini.cv$prev_diff = gini.cv$vaginal_prev-gini.cv$caesarean_prev
  gini.cv$prev_diff_strict = gini.cv$vaginal_prev_strict-gini.cv$caesarean_prev_strict

  gini.cv$enriched_vaginal = gini.cv$prev_diff>0
  gini.cv$enriched_vaginal_strict = gini.cv$prev_diff_strict>0

  gini.cv$rf_rank = 1:nrow(gini.cv)
  gini.cv$rf_top20 = as.numeric(gini.cv$rf_rank<=20)
  gini.cv$rf_n = n_per_group

  # add plot colour scheme
  gini.cv$fisher_col = "ns"
  gini.cv$fisher_col[gini.cv$prev_diff>0 & gini.cv$fisher_padj<0.1] = "vaginal"
  gini.cv$fisher_col[gini.cv$prev_diff<0 & gini.cv$fisher_padj<0.1] = "caesarean"

  gini.cv$fisher_col_strict = "ns"
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict>0 & gini.cv$fisher_padj_strict<0.1] = "vaginal"
  gini.cv$fisher_col_strict[gini.cv$prev_diff_strict<0 & gini.cv$fisher_padj_strict<0.1] = "caesarean"

  # add cross-val mean + sd
  gini.cv$crossval_mean = mean(stat.cv$accuracy)
  gini.cv$crossval_sd = sd(stat.cv$accuracy)

  # save RF results
  write.csv(stat.cv,file=paste0("output_module6/RF_outputs/crossval_stats_deliverymode_",subset,"_",level,".csv"))
  write.csv(gini.cv, file=paste0("output_module6/RF_outputs/crossval_gini_deliverymode_",subset,"_",level,".csv"))
  return(list(gini.cv=gini.cv))
}

# function to create prevalence/abundance plot
rf_prev_plot_delivery = function(input_df=gini.cv, strict=TRUE, species=FALSE) {
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
  
  ggplot(input_top20, aes(x = taxonomy, y = vaginal_prev)) +
    scale_x_discrete(limits = rev(levels(input_top20$taxonomy))) +
    xlab("") + ylab("prevalence") +
    geom_point(aes(colour = "vaginal", size = vaginal_mean), stat = "identity", alpha=0.7) +
    geom_point(aes(y = caesarean_prev, colour = "caesarean", size = caesarean_mean), stat = "identity", alpha=0.7) +
    coord_flip() + theme_bw() + ylim(0,1) + scale_y_continuous(breaks=c(0,0.5,1)) +
    theme(axis.text.y = element_blank(), legend.position = "right") + labs(colour = "enriched group", size = "mean \nabundance") +
    scale_colour_manual(values=c(cs_col,vd_col)) + scale_size(limits = c(0, 0.5)) +
    guides(colour = guide_legend(order=1), size = guide_legend(order=2))
  #guides(colour = "none")
  
}


