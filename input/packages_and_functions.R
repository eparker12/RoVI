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
library(ZIBR)
library(corrplot)
library(data.table)
library(ggpubr)
library(labdsv)
library(UpSetR)
library(crossval)
library(DESeq2)
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
library(ggtree)
library(decontam)
library(inlmisc)
library(NBZIMM)
library(formattable)
library(ggExtra)
library(sjstats)

theme_set(theme_bw())

# set country colours
India_col = "#3f9d7f" #"#66c2a5"
India_OPV_col = "#3f9d7f"
India_IPV_col = "#045a8d"
India_col_2 = "#66c2a5" #"#045a8d"
Malawi_col = "#fc8d62"
UK_col = "#8da0cb"
r_col = "#045a8d"
nr_col = "#d95f02"

########################################
### Functions for RoVI analyses work ###
########################################

### Module 1

# function to extract n characters from right of string
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }

# function to plot % retention by xvariable
retention_plot = function(yvariable, xvariable, ylabel) {
  df = metadata
  df$yvar = metadata[,yvariable] 
  df$xvar = metadata[,xvariable] 
  ggplot(df, aes(xvar, yvar/input, colour=xvar)) + geom_jitter(size = 0.5, alpha = 0.1, width=0.2) +
    geom_violin(alpha = 0.5) + theme_bw() + ylab(ylabel) + xlab("") +
    ylim(0,1) + theme(legend.position = "none", strip.background = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# function to pick out potential contaminants
contam_freq = function(ps_t) {
  # create subsets of input data for each country
  ps_t_I = subset_samples(ps_t, country=="India")
  ps_t_M = subset_samples(ps_t, country=="Malawi")
  ps_t_U = subset_samples(ps_t, country=="UK")
  
  # frequency-based contaminant classification for each country subset
  contam_freqI = isContaminant(ps_t_I, conc="nanodrop_ngul", threshold=0.1, detailed=TRUE, normalize=TRUE, method="frequency")
  contam_freqM = isContaminant(ps_t_M, conc="nanodrop_ngul", threshold=0.1, detailed=TRUE, normalize=TRUE, method="frequency")
  contam_freqU = isContaminant(ps_t_U, conc="nanodrop_ngul", threshold=0.1, detailed=TRUE, normalize=TRUE, method="frequency")
  
  # in several instances, taxa may be identified as contaminants based a small number of observations in a single country, 
  # despite being more prevalent and showing a non-contaminant profile in another country; 
  # to account for this, set contaminant p value to NA if <10 observations
  contam_freqI$p[contam_freqI$prev<10] = NA
  contam_freqM$p[contam_freqM$prev<10] = NA
  contam_freqU$p[contam_freqU$prev<10] = NA
  
  # combine results
  combined_res = data.frame(taxonomy = data.frame(as(tax_table(ps_t), "matrix"))$taxonomy,
                            prev_I = contam_freqI$prev, p_I = contam_freqI$p,
                            prev_M = contam_freqM$prev, p_M = contam_freqM$p,
                            prev_U = contam_freqU$prev, p_U = contam_freqU$p)
  combined_res$pmin = pmin(combined_res$p_I, combined_res$p_M, combined_res$p_U, na.rm=T)
  combined_res$contaminant = combined_res$pmin<0.1
  combined_res$contaminant[is.na(combined_res$pmin)] = FALSE
  return(combined_res)
}

# add paired and random distances to phyloseq sample data
add_paired_distances = function(input_df, input_dist) {
  input_df$sample_ID_full = as.character(input_df$sample_ID_full)
  
  for (i in 1:nrow(input_df)) {
    ID = input_df$sample_ID_full[i]
    row = input_df[i,]
    
    # if MS1, BM1, or ctrl, turnover = NA
    if (row$sample_type == "MS1" | row$sample_type == "BM1" | row$control == "ctrl") {
      is.na(input_df$turnover[i])
      is.na(input_df$turnover_randompair[i])
    }
    
    # if BS1, turnover relative to MS1 (paired and random)
    if (row$sample_type == "BS1") {
      family = subset(input_df, family_ID==row$family_ID)
      if ("MS1" %in% family$sample_type) {
        prev_ID = family$sample_ID_full[family$sample_type=="MS1"]
        input_df$turnover[i] = input_dist[ID, prev_ID]
        random = subset(input_df, family_ID!=row$family_ID & sample_type=="MS1" & country==row$country)
        prev_ID_random = sample(random$sample_ID_full,1)
        input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
      } else {
        is.na(input_df$turnover[i])
        is.na(input_df$turnover_randompair[i])
      }
    }
    
    # if BS2, turnover relative to BS1 (paired and random)
    if (row$sample_type == "BS2") {
      family = subset(input_df, family_ID==row$family_ID)
      if ("BS1" %in% family$sample_type) {
        prev_ID = family$sample_ID_full[family$sample_type=="BS1"]
        input_df$turnover[i] = input_dist[ID, prev_ID]
        random = subset(input_df, family_ID!=row$family_ID & sample_type=="BS1" & country==row$country)
        prev_ID_random = sample(random$sample_ID_full,1)
        input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
      } else {
        is.na(input_df$turnover[i])
        is.na(input_df$turnover_randompair[i])
      }
    }
    
    # if BS3, turnover relative to BS2 (paired and random)
    if (row$sample_type == "BS3") {
      family = subset(input_df, family_ID==row$family_ID)
      if ("BS2" %in% family$sample_type) {
        prev_ID = family$sample_ID_full[family$sample_type=="BS2"]
        input_df$turnover[i] = input_dist[ID, prev_ID]
        random = subset(input_df, family_ID!=row$family_ID & sample_type=="BS2" & country==row$country)
        prev_ID_random = sample(random$sample_ID_full,1)
        input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
      } else {
        is.na(input_df$turnover[i])
        is.na(input_df$turnover_randompair[i])
      }
    }
    
    # if BS5, turnover relative to BS3 (paired and random)
    if (row$sample_type == "BS5") {
      family = subset(input_df, family_ID==row$family_ID)
      if ("BS3" %in% family$sample_type) {
        prev_ID = family$sample_ID_full[family$sample_type=="BS3"]
        input_df$turnover[i] = input_dist[ID, prev_ID]
        random = subset(input_df, family_ID!=row$family_ID & sample_type=="BS3" & country==row$country)
        prev_ID_random = sample(random$sample_ID_full,1)
        input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
      } else {
        is.na(input_df$turnover[i])
        is.na(input_df$turnover_randompair[i])
      }
    }
  }
  
  # if BM2, turnover relative to BM1 (paired and random)
  if (row$sample_type == "BM2") {
    family = subset(input_df, family_ID==row$family_ID)
    if ("BM1" %in% family$sample_type) {
      prev_ID = family$sample_ID_full[family$sample_type=="BM1"]
      input_df$turnover[i] = input_dist[ID, prev_ID]
      random = subset(input_df, family_ID!=row$family_ID & sample_type=="BM1" & country==row$country)
      prev_ID_random = sample(random$sample_ID_full,1)
      input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
    } else {
      is.na(input_df$turnover[i])
      is.na(input_df$turnover_randompair[i])
    }
  }
  
  # if BM3, turnover relative to BM2 (paired and random)
  if (row$sample_type == "BM3") {
    family = subset(input_df, family_ID==row$family_ID)
    if ("BM3" %in% family$sample_type) {
      prev_ID = family$sample_ID_full[family$sample_type=="BM2"]
      input_df$turnover[i] = input_dist[ID, prev_ID]
      random = subset(input_df, family_ID!=row$family_ID & sample_type=="BM2" & country==row$country)
      prev_ID_random = sample(random$sample_ID_full,1)
      input_df$turnover_randompair[i] = input_dist[ID, prev_ID_random]
    } else {
      is.na(input_df$turnover[i])
      is.na(input_df$turnover_randompair[i])
    }
  }
  
  return(input_df)
}

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

# function to plot taxonomic distribution pre/post filtering
taxa_filter_summary = function(ps_unfiltered, ps_filtered, level) {
  # taxa list pre/post filtering
  taxa_list_u = data.frame(as(tax_table(ps_unfiltered),"matrix"))
  taxa_list_f = data.frame(as(tax_table(ps_filtered),"matrix"))
  
  # calculate distribution of taxa by class - unfiltered
  u = data.frame(table(taxa_list_u[,level]), row.names="Var1")
  u$proportion = round(u$Freq/sum(u$Freq),4)
  u = u[order(-u$Freq),][1:10,]
  u["other",] = c(nrow(taxa_list_u)-sum(u$Freq), 1-sum(u$proportion))
  u$type = rep("unfiltered_ps1", nrow(u))
  u$taxon = rownames(u)
  
  # calculate distribution of taxa by class - filtered
  f = data.frame(table(taxa_list_f[,level]), row.names="Var1")
  f$proportion = round(f$Freq/sum(f$Freq),4)
  f = f[order(-f$Freq),][1:10,]
  f["other",] = c(nrow(taxa_list_f)-sum(f$Freq), 1-sum(f$proportion))
  f$type = rep("filtered_ps7", nrow(f))
  f$taxon = rownames(f)
  
  # plot
  t = c(brewer.pal(n=9, "Set3"), brewer.pal(n=9, "Set1"), brewer.pal(n=8, "Set2"))
  plot_u = ggplot(u, aes(y=Freq, x=type, fill=taxon)) + geom_bar(stat="identity") + xlab("") + ylab("frequency") +
    scale_fill_manual(values = c(Actinobacteria = t[1], Alphaproteobacteria = t[2], Bacilli = t[3], 
                                 Bacteroidia = t[4],  Clostridia = t[5], Coriobacteriia = t[6], Erysipelotrichia = t[7], 
                                 Gammaproteobacteria = t[8], Mollicutes = t[11], Negativicutes = t[10], other = t[9],
                                 unassigned_Bacteria = t[12], Oxyphotobacteria = t[13], unassigned = t[14], 
                                 unassigned_Eukaryota = t[15], Saccharimonadia = t[16], 
                                 Deltaproteobacteria = t[17], Deinococci = t[18], Fusobacteriia = t[19])) 
  plot_f = ggplot(f, aes(y=Freq, x=type, fill=taxon)) + geom_bar(stat="identity") + xlab("") + ylab(" ") +
    scale_fill_manual(values = c(Actinobacteria = t[1], Alphaproteobacteria = t[2], Bacilli = t[3], 
                                 Bacteroidia = t[4],  Clostridia = t[5], Coriobacteriia = t[6], Erysipelotrichia = t[7],
                                 Gammaproteobacteria = t[8], Mollicutes = t[11], Negativicutes = t[10], other = t[9],
                                 unassigned_Bacteria = t[12], Oxyphotobacteria = t[13], unassigned = t[14], 
                                 unassigned_Eukaryota = t[15], Saccharimonadia = t[16], 
                                 Deltaproteobacteria = t[17], Deinococci = t[18], Fusobacteriia = t[19])) 
  return(grid.arrange(plot_u, plot_f, ncol=2))
}

# function to generate correlation plot for 10% validation genus abundance correlations
filter_plot = function(taxon) {
  df_CGR_Imp$taxon = df_CGR_Imp[,taxon]
  df_CGR_Imp$taxon.1 = df_CGR_Imp[,paste0(taxon,".1")]
  ggplot(df_CGR_Imp, aes(x=taxon, y=taxon.1)) + geom_abline(color = "black", linetype = "dotted") + 
    geom_point(aes(colour=sample_type), size = 2, alpha = 0.6) +
    scale_colour_manual(values = c("BS1" = "#fe9929", "BS5" = "#cc4c02", "MS1" = "#662506", "BM1" = "#2c7fb8")) +
    xlab("") +  ylab("") + theme_bw() + xlim(0,1) + ylim(0,1)+
    theme(legend.position = "none") + ggtitle(taxon) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), plot.title = element_text(size=10)) 
}

# function to create multipanel plot for 10% validation genus abundance correlations
multi_plot = function(df_CGR_Imp) {
  p1 = filter_plot(names(df_CGR_Imp)[1])
  p2 = filter_plot(names(df_CGR_Imp)[2])
  p3 = filter_plot(names(df_CGR_Imp)[3])
  p4 = filter_plot(names(df_CGR_Imp)[4])
  p5 = filter_plot(names(df_CGR_Imp)[5])
  p6 = filter_plot(names(df_CGR_Imp)[6])
  p7 = filter_plot(names(df_CGR_Imp)[7])
  p8 = filter_plot(names(df_CGR_Imp)[8])
  p9 = filter_plot(names(df_CGR_Imp)[9])
  p10 = filter_plot(names(df_CGR_Imp)[10])
  p11 = filter_plot(names(df_CGR_Imp)[11])
  p12 = filter_plot(names(df_CGR_Imp)[12])
  p13 = filter_plot(names(df_CGR_Imp)[13])
  p14 = filter_plot(names(df_CGR_Imp)[14])
  p15 = filter_plot(names(df_CGR_Imp)[15])
  p16 = filter_plot(names(df_CGR_Imp)[16])
  p17 = filter_plot(names(df_CGR_Imp)[17])
  p18 = filter_plot(names(df_CGR_Imp)[18])
  p19 = filter_plot(names(df_CGR_Imp)[19])
  p20 = filter_plot(names(df_CGR_Imp)[20])
  if (names(df_CGR_Imp)[2]=="Escherichia.Shigella") { p2 = p2 + ggtitle("Escherichia/Shigella") } else { p11 = p11 + ggtitle("Escherichia/Shigella") }
  return(grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,ncol=4))
}





### Module 2

# plot categorical outcomes by country or arm
fig_cat = function(df = s, outcome, title) {
  # collate data by country unless input data consists only of Indian samples
  if (all(df$country=="India")) {
    # if all observations absent, manually create figure dataframe
    if (sum(df[,outcome], na.rm=TRUE)==0) { 
      t = data.frame(nonresponders = as.numeric(table(df$country_arm)), responders = c(0,0))
      rownames(t) = c("India (IPV)", "India (OPV)")
    } else {
      t = data.frame(nonresponders = table(df$country_arm, df[,outcome])[,1], responders = table(df$country_arm, df[,outcome])[,2])
    }
  } else {
    t = data.frame(nonresponders = table(df$country, df[,outcome])[,1], responders = table(df$country, df[,outcome])[,2])
  }
  t$covariate = rownames(t)
  t$n = t$nonresponders + t$responders
  t$proportion = t$responders/t$n
  t$lowerCI = binom.confint(t$responders, t$n, method="exact")$lower
  t$upperCI = binom.confint(t$responders, t$n, method="exact")$upper
  ggplot(t, aes(x = factor(covariate), y = proportion, fill=factor(covariate))) + geom_bar(stat="identity", alpha=0.8) + 
    scale_fill_manual(values = c("India" = India_col, "India (OPV)" = India_col, "India (IPV)" = India_IPV_col, "Malawi" = Malawi_col, "UK" = UK_col)) +
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.3) + 
    xlab("") + theme(legend.position = "none") + ggtitle(title) + scale_y_continuous(limits=c(0,1.125), breaks=c(0,0.25,0.5,0.75,1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(size=10, hjust = 0.5))
}

# summarise binary outcome by country or arm
outcome_summary = function(df, outcome) { 
  d1 = df
  d1$outcome = d1[,outcome]
  d1 = subset(d1, !is.na(outcome))
  # collate data by country unless input data consists only of Indian samples
  if (all(df$country=="India")) {
    data.frame(proportions = levels(factor(d1$country_arm)),
               positive = as.vector(table(d1$outcome, d1$country_arm)[2,]),
               n = as.vector(table(d1$country_arm)), 
               percent = as.vector(round(table(d1$outcome, d1$country_arm)[2,]/table(d1$country_arm)*100,1)))
  } else {
    data.frame(proportions = levels(factor(d1$country)),
               positive = as.vector(table(d1$outcome, d1$country)[2,]),
               n = as.vector(table(d1$country)), 
               percent = as.vector(round(table(d1$outcome, d1$country)[2,]/table(d1$country)*100,1)))
  }
}

# FDR-corrected p values for binary outcomes (Fisher's exact test)
p_summary = function(df, outcome) {
  # collate data by country unless input data consists only of Indian samples
  if (all(df$country=="India")) {
    data.frame( p_values = c("India (IPV) vs India (OPV)"), p = format(c(round(fisher.test(table(df[,outcome], df$country_arm))$p.value,4)),nsmall=4))
  } else {
    data.frame( p_values = c("India vs Malawi", "India vs UK", "Malawi vs UK"), 
                p = format(c(round(fisher.multcomp(table(df[,outcome], df$country))$p.value,4)),nsmall=4))
  }
}

# plot continuous outcomes on log10 scale by country or arm
fig_cont_log10 = function(df, outcome, ylabel, exposure=FALSE) {
  df$outcome = df[,outcome]
  if (all(df$country=="India")) {
    if (exposure==FALSE) { df$xvar = df$country_arm } else { df$xvar = df$country_exposure }
  } else { df$xvar = df$country }

  ggplot(df, aes(factor(xvar), y=outcome, colour=factor(xvar))) + geom_jitter(size = 0.5, alpha = 0.5, width=0.2) +
    geom_boxplot(alpha = 0.5) + ylab(ylabel) + xlab("") + theme(legend.position = "none", strip.background = element_blank()) + 
    scale_colour_manual(values = c("India" = India_col, "India (OPV)" = India_col,"India (IPV)" = India_IPV_col, 
                                   "India (exposed)" = India_col,"India (unexposed)" = India_IPV_col, "Malawi" = Malawi_col, "UK" = UK_col)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(trans="log10", breaks = c(0,1,10,100,1000)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(size=10, hjust = 0.5))
}


# plot continuous outcomes on log10 scale by delivery mode
fig_cont_log10_delivery = function(df, outcome, ylabel) {
  df$outcome = df[,outcome]
  ggplot(df, aes(factor(mode_delivery), y=outcome, colour=factor(mode_delivery))) + geom_jitter(size = 0.5, alpha = 0.5, width=0.2) +
    geom_boxplot(alpha = 0.5) + ylab(ylabel) + xlab("") + theme(legend.position = "none", strip.background = element_blank()) + 
    scale_colour_manual(values = c("#0570b0", "#cc4c02")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(plot.title = element_text(size=10, hjust = 0.5)) 
}

# summarise gmt by country or arm
gmt_summary = function(df, outcome, exposure=FALSE) { 
  df$outcome = df[,outcome]
  df = subset(df, !is.na(outcome))
  
  if (all(df$country=="India")) {
    # default is to compare by arm in Indian infants, but specifying exposure as TRUE will amend to neonatal exposure status
    if (exposure==FALSE) {
      gmts = format(round(rbind(
        Gmean(df$outcome[df$country_arm=="India (IPV)"], conf.level=0.95),
        Gmean(df$outcome[df$country_arm=="India (OPV)"], conf.level=0.95)),0),nsmall=0)
      data.frame(cbind(gmts = levels(factor(df$country_arm)), n = as.vector(table(df$country_arm)), gmts))
    } else {
      gmts = format(round(rbind(
        Gmean(df$outcome[df$country_exposure=="India (exposed)"], conf.level=0.95),
        Gmean(df$outcome[df$country_exposure=="India (unexposed)"], conf.level=0.95)),0),nsmall=0)
      data.frame(cbind(gmts = levels(factor(df$country_exposure)), n = as.vector(table(df$country_exposure)), gmts))
    }
  } else {
    gmts = format(round(rbind(
      Gmean(df$outcome[df$country=="India"], conf.level=0.95),
      Gmean(df$outcome[df$country=="Malawi"], conf.level=0.95),
      Gmean(df$outcome[df$country=="UK"], conf.level=0.95)),0),nsmall=0)
    data.frame(cbind(gmts = levels(factor(df$country)), n = as.vector(table(df$country)), gmts))
  }
}

# calculate p values for continuous outcomes on log10 scale (anova with Tukey correction)
p_summary_gmt = function(df, outcome, exposure=FALSE) {
  if (all(df$country=="India")) {
    # default is to compare by arm in Indian infants, but specifying exposure as TRUE will amend to neonatal exposure status
    if (exposure==FALSE) {
      data.frame(p_values = "India (IPV) vs India (OPV)", p = format(round(as.vector(tail(unlist(TukeyHSD(x=aov(log(df[,outcome]) ~ df$country_arm), "df$country_arm", conf.level=0.95)[1]),1)),4),nsmall=4))
    } else {
      data.frame(p_values = "India (exposed) vs India (unexposed)", p = format(round(as.vector(tail(unlist(TukeyHSD(x=aov(log(df[,outcome]) ~ df$country_exposure), "df$country_exposure", conf.level=0.95)[1]),1)),4),nsmall=4))
    }
  } else {
      data.frame(p_values = c("India vs Malawi", "India vs UK", "Malawi vs UK"), 
               p = format(round(as.vector(tail(unlist(TukeyHSD(x=aov(log(df[,outcome]) ~ df$country), "df$country", conf.level=0.95)[1]),3)),4),nsmall=4))
  }
}

# estimate correlation between covariate and specified outcome measure
corr_function = function(df = s, outcome_measure, covariate, country_exposure_subset, return_option=c("corr", "p")) {
  if (country_exposure_subset=="India") { t = subset(df, country==country_exposure_subset) } else { t = subset(df, country_exposure==country_exposure_subset) }
  
  # correlation coefficients for Shannon
  if (covariate=="Shannon" | covariate=="Shannon_rsv") {
    cor = cor.test(log(t[,outcome_measure]), t[,covariate], method = "pearson")
  
  # correlation coefficients for richness
    } else if (covariate=="Observed" | covariate=="Observed_rsv") {
    cor = cor.test(log(t[,outcome_measure]), log(t[,covariate]), method = "pearson")
  
  # correlation coefficients for other covariates
  } else {
    cor = cor.test(log(t[,outcome_measure]), log(t[,covariate]), method = "pearson")
  }
  if (return_option=="corr") { return(format(round(cor$estimate,4),nsmall=4)) }
  if (return_option=="p") { return(format(round(cor$p.value,4),nsmall=4)) }
}

# calculate correlation coefficients for country subset
corr_by_country = function(outcome_measure, country_exposure_subset, return_val, antibody=TRUE) {
  if (antibody==TRUE) {
    res = c(corr_function(s, outcome_measure, "MB1_IgA", country_exposure_subset, return_option=return_val),
            corr_function(s, outcome_measure, "BM1_IgA", country_exposure_subset, return_option=return_val))
    if (country_exposure_subset != "UK" & country_exposure_subset != "Malawi") { 
      res = c(res,
              corr_function(s, outcome_measure, "MB1_IgG", country_exposure_subset, return_option=return_val),
              corr_function(s, outcome_measure, "CB1_IgG", country_exposure_subset, return_option=return_val),
              corr_function(s, outcome_measure, "BB1_IgG", country_exposure_subset, return_option=return_val))
    } else { (res = c(res, NA, NA, NA)) } 
  } else {
    res = c(corr_function(s, outcome_measure, "A1AT_BS3_Concugml", country_exposure_subset, return_option=return_val),
            corr_function(s, outcome_measure, "A1AT_BS5_Concugml", country_exposure_subset, return_option=return_val),
            corr_function(s, outcome_measure, "MPO_BS3_Concngml", country_exposure_subset, return_option=return_val),
            corr_function(s, outcome_measure, "MPO_BS5_Concngml", country_exposure_subset, return_option=return_val))
    if (country_exposure_subset == "UK") { res= c(res,NA) } else { res = c(res,corr_function(s, outcome_measure, "AG_Concugml", country_exposure_subset, return_option=return_val)) }
  }
}

# effect of neonatal exposure on binary outcome
cofactor_summary = function(df, outcome, cofactor) { 
  d1 = df
  d1$outcome = d1[,outcome]
  d1$cofactor = d1[,cofactor]
  d1 = subset(d1, !is.na(outcome) & !is.na(cofactor))
  data.frame(proportions = levels(factor(d1$cofactor)),
             positive = as.vector(table(d1$outcome, d1$cofactor)[2,]),
             n = as.vector(table(d1$cofactor)), 
             percent = as.vector(round(table(d1$outcome, d1$cofactor)[2,]/table(d1$cofactor)*100,1)),
             p = format(c(round(fisher.multcomp(table(d1$outcome, d1$cofactor))$p.value,4), NA), nsmall=4))
}

# plot categorical outcome by binary cofactor
fig_cat_cofactor = function(df, outcome, title, cofactor) {
  df$cofactor = df[,cofactor]
  t = data.frame(nonresponders = table(df$cofactor, df[,outcome])[,1], responders = table(df$cofactor, df[,outcome])[,2])
  t$covariate = rownames(t)
  t$n = t$nonresponders + t$responders
  t$proportion = t$responders/t$n
  t$lowerCI = binom.confint(t$responders, t$n, method="exact")$lower
  t$upperCI = binom.confint(t$responders, t$n, method="exact")$upper
  ggplot(t, aes(x = factor(covariate), y = proportion, fill=factor(covariate))) + geom_bar(stat="identity", alpha=0.8) + 
    scale_fill_manual(values = c("India (exposed)" = India_col, "India (unexposed)" = India_col_2, "Malawi" = Malawi_col,
                                 "UK" = UK_col, "neonatal+" = "#e6550d", "neonatal-" = "#43a2ca", "vaginal" = "#43a2ca", "caesarean" = "#e6550d", "Sab (0/1)" = "#43a2ca", "Sab (both)" = "#e6550d")) +
    geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.3) + 
    xlab("") + theme(legend.position = "none") + ggtitle(title) + scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
    theme(plot.title = element_text(size=10, hjust = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}





### Module 3

# function to plot prevalence vs mean abundance (when present) for sample subset
taxon_abundance_function = function(input_rps, genus=TRUE, input_country="India", input_colour=India_col, MS1=FALSE) {
  
  # select subset of samples (infant or maternal)
  if (MS1==FALSE) {
    sample_list = rownames(sample_data(input_rps)[sample_data(input_rps)$sample_type!="MS1" & sample_data(input_rps)$country==input_country])
  } else {
    sample_list = rownames(sample_data(input_rps)[sample_data(input_rps)$sample_type=="MS1" & sample_data(input_rps)$country==input_country])
  }
  rps_t = prune_samples(sample_list, input_rps) %>% prune_taxa(taxa_sums(.) > 0, .)
  
  # set 0s to NA so that these do not contribute to mean abundance calculations
  rps_t_no0 = rps_t
  otu_table(rps_t_no0)[otu_table(rps_t_no0)==0] = NA
  
  # calculate taxon means
  tax_abund = data.frame(round(rowMeans(otu_table(rps_t_no0), na.rm=T),6), round(rowSds(otu_table(rps_t_no0), na.rm=T),6), rowSums(otu_table(rps_t)>0))
  names(tax_abund) = c("mean_abundance", "SD_abundance", "present")
  tax_abund$prevalence = round(tax_abund$present/ncol(otu_table(rps_t)),3)
  
  # append taxon names
  if (genus==TRUE) {
    tax_abund$taxon = data.frame(as(tax_table(rps_t),"matrix"))$Genus
  } else { 
    tax_abund$taxon = data.frame(as(tax_table(rps_t),"matrix"))$taxonomy_species
  }
  
  # order by decreasing abundance, add rank
  tax_abund = tax_abund[order(-tax_abund$mean_abundance),]
  tax_abund$rank = 1:nrow(tax_abund)
  N_over_20_prev = sum(tax_abund$prevalence>=0.20)
  N_total = nrow(tax_abund)
  
  # Calculate cumulative abundance distribution
  tax_abund$cumulative = NA
  tax_abund$cumulative[1] = tax_abund$mean_abundance[1] 
  for (i in 2:nrow(tax_abund)) { tax_abund$cumulative[i] = tax_abund$cumulative[i-1] + tax_abund$mean_abundance[i] }
  tax_abund$label = tax_abund$taxon
  tax_abund$label[tax_abund$prevalence<0.5] = NA
  ggplot(tax_abund, aes(x=prevalence, y=mean_abundance, label=label)) +
    geom_point(alpha=0.5, size=2, colour=input_colour) + 
    geom_text(size = 3, colour=input_colour, check_overlap = TRUE, nudge_y=0.01) + 
    ylab("mean relative abundance (when present)") +
    ylim(0,0.5) +
    scale_x_continuous(limits = c(0,1.1), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    annotate("text", x = 0.55, y = 0.5, label = paste0(N_over_20_prev, "/", N_total, " taxa with prevalence ≥20%"), colour=input_colour,  fontface = 'italic')
}

# linear mixed-effects models - country
lme_alpha = function(df, outcome, comparison="country") {
  lm_base <- lmer(formula(paste0(outcome," ~ (1 | family_ID)")), df, REML=FALSE)
  lm_week <- lmer(formula(paste0(outcome," ~  week + (1|family_ID)")), df, REML=FALSE)
  lm_outcome <- lmer(formula(paste0(outcome," ~  week + ",comparison," + (1|family_ID)")), df, REML=FALSE)
  anova(lm_outcome, lm_week)$Pr[2]
}

# function to calculate genus means by country and sample type
genus_mean_cal = function(df, type, country_subset) {
  if (country_subset=="India") {
    g1 = subset(df, sample_type==type & (country_arm=="India (OPV)" | country_arm=="India (IPV)"))
  } else { g1 = subset(df, sample_type==type & country_arm==country_subset) }
  colMeans(g1[,1:(ncol(g1)-2)])
}

# function to pick out subset of samples and plot PCoA via chosen method (weighted or unweighted)
beta_plot = function(ps_input, subset, plot_title, method=c("Bray-Curtis (weighted)", "Bray-Curtis (unweighted)")) {
  # pick out selected subset of samples
  t = data.frame(sample_data(ps_input))
  samplelist = subset(t, sample_type==subset)$sample_ID_full
  ps_input = prune_samples(samplelist, ps_input)
  beta = data.frame(sample_data(ps_input))
  
  # calculate beta diversity distances for selected samples
  if (method == "Bray-Curtis (weighted)") { beta_dist = phyloseq::distance(ps_input, method = "bray") }
  if (method == "Bray-Curtis (unweighted)") { 
    otu_table(ps_input)[otu_table(ps_input)<0.001] = 0
    beta_dist = phyloseq::distance(ps_input, method = "bray", binary=TRUE) 
  }
  
  # ordinate and extract PCs
  beta_ord = ordinate(ps_input, method = "PCoA", distance = beta_dist)
  beta$PC1 = beta_ord$vectors[,1]
  beta$PC2 = beta_ord$vectors[,2]
  
  # create PCoA plot
  ggplot(beta, aes(PC1, PC2, colour=country)) + geom_point(size = 1.5, alpha = 0.5) +
    theme(legend.position="none", legend.title = element_blank()) + ggtitle(plot_title) +
    scale_color_manual(values=c(India = India_col, Malawi = Malawi_col, UK = UK_col)) + 
    xlab(paste0("Bray-Curtis, PC1 (", format(round(beta_ord$values[1,1],1), nsmall=1), "%)")) + 
    ylab(paste0("Bray-Curtis, PC2 (", format(round(beta_ord$values[2,1],1), nsmall=1), "%)")) +
    theme(axis.title=element_text(size=10), plot.title=element_text(size=10))
  
}

# function to plot taxonomic enrichment between countries (based on presence/absence)
summary_prev_plot= function(country1,country2, colour1, colour2, level, threshold=0.1, yshift=2, ymax=85, strict=TRUE) {
  BS1 = read.csv(paste0("output_module3/RF_output_importance/crossval_gini_",country1,"_",country2,"_BS1_",level,".csv"), row.names=1)
  BS2 = read.csv(paste0("output_module3/RF_output_importance/crossval_gini_",country1,"_",country2,"_BS2_",level,".csv"), row.names=1)
  BS3 = read.csv(paste0("output_module3/RF_output_importance/crossval_gini_",country1,"_",country2,"_BS3_",level,".csv"), row.names=1)
  BS5 = read.csv(paste0("output_module3/RF_output_importance/crossval_gini_",country1,"_",country2,"_BS5_",level,".csv"), row.names=1)
  
  # specify whether to use standard presence/absence data or ≥0.1% threshold (strict==TRUE)
  if (strict==TRUE) {
    BS1$fisher_padj = BS1$fisher_padj_strict
    BS2$fisher_padj = BS2$fisher_padj_strict
    BS3$fisher_padj = BS3$fisher_padj_strict
    BS5$fisher_padj = BS5$fisher_padj_strict
    
    BS1$enriched_group1 = BS1$enriched_group1_strict
    BS2$enriched_group1 = BS2$enriched_group1_strict
    BS3$enriched_group1 = BS3$enriched_group1_strict
    BS5$enriched_group1 = BS5$enriched_group1_strict
  }
  
  # remove taxa with NA for adjusted p value
  BS1 = subset(BS1, !is.na(fisher_padj))
  BS2 = subset(BS2, !is.na(fisher_padj))
  BS3 = subset(BS3, !is.na(fisher_padj))
  BS5 = subset(BS5, !is.na(fisher_padj))
  
  # create dataframe summarising taxonomic enrichment results
  fisher_summary = data.frame(
    count = c(sum(BS1$enriched_group1==TRUE & BS1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS2$enriched_group1==TRUE & BS2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS3$enriched_group1==TRUE & BS3$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS5$enriched_group1==TRUE & BS5$fisher_padj<threshold, na.rm=TRUE),
              sum(BS1$enriched_group1==FALSE & BS1$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS2$enriched_group1==FALSE & BS2$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS3$enriched_group1==FALSE & BS3$fisher_padj<threshold, na.rm=TRUE), 
              sum(BS5$enriched_group1==FALSE & BS5$fisher_padj<threshold, na.rm=TRUE)),
    enriched = c(rep(country1,4), rep(country2,4)),
    week = factor(rep(c("1", "4", "6/8", "10/12"),2), levels = c("1", "4", "6/8", "10/12")),
    total = c(nrow(BS1), nrow(BS2), nrow(BS3), nrow(BS5))
  )
  
  # generate plot
  ggplot(fisher_summary, aes(x=week, y=count, label=count, colour=enriched, group=enriched)) + 
    geom_point(size=3, alpha=0.8) + geom_line(alpha=0.5, show.legend = FALSE) + 
    geom_text(vjust = 0, nudge_y = yshift, show.legend = FALSE) +
    geom_text(aes(y = rep(ymax,8), label=total), colour="grey", show.legend = FALSE, fontface = "italic") +
    ylim(0,ymax) +
    scale_colour_manual(values=c(colour1, colour2)) + 
    ggtitle(paste0(country1,"-",country2)) + theme(legend.position="none") +
    ylab("number of enriched taxa") + labs(colour='enriched country')
}

# function to create genus abundance plots by country
genus_plot = function(g, genus) {
  g$week[g$week==8] = 6
  g$week[g$week==12] = 10
  g = subset(g, sample_type!="MS1")
  g$genus = g[,genus]
  ggplot(g, aes(x=week, y=genus, color=factor(country))) + geom_jitter(size=0.5,alpha=0.07, width = 0.2) +
    xlab("") + ylab("") + ggtitle(genus) + labs(colour="") +
    scale_x_continuous(breaks=c(1,4,6,10), labels=c("1", "4", "6/8", "10/12")) +
    theme_bw() + geom_smooth(method = "loess", size=0.5, alpha=0.1, aes(fill=factor(country))) + 
    scale_color_manual(values = c(India = India_col, UK = UK_col, Malawi=Malawi_col)) + scale_fill_manual(values = c(India = India_col, UK = UK_col, Malawi=Malawi_col)) +
    theme(legend.position = "none", title=element_text(size=10))
}





### Module 4

# longitudinal plot of alpha diversity
alpha_longitudinal = function(t1, alpha_measure, alpha_title, outcome, outcome_title, exposure=FALSE) {
  # if exposure==TRUE, restrict analyses to Indian infants, stratified by neonatal exposure status
  if (exposure==TRUE) {
    t1 = subset(t1, country_exposure=="India (exposed)" | country_exposure=="India (unexposed)")
    t1$country = t1$country_exposure
    ggtitle = t1$country[1]
  } else { ggtitle = t1$country[1] }
  
  # select alpha diversity metric and outcome measure
  t1$alpha_measure = t1[,alpha_measure]
  t1$outcome = t1[,outcome]
  
  # create plot
  ggplot(t1, aes(x = factor(week_fig), y = alpha_measure, group=family_ID, colour=factor(outcome))) + 
    geom_line(alpha=0.04) + scale_color_manual(outcome_title, values = c(nr_col, r_col)) + 
    stat_summary(aes(group = outcome), geom = "line", fun.y = mean, size = 1, alpha = 0.8) +
    stat_summary(aes(group = outcome), geom = "point", fun.data = "mean_cl_boot", size = 2, alpha = 0.8) + 
    ggtitle("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(breaks=c("week 1","week 4","week 6/8","week 10/12"), labels=c("1", "4", "6/8", "10/12")) + 
    xlab("week") + ylab(alpha_title) + 
    theme_bw() + theme(legend.position = "none", plot.title = element_text(size=10, hjust = 0.5)) +
    facet_grid(.~country) + theme(strip.background = element_blank())
}

# linear mixed-effects models - country
lme_alpha_ORV = function(df, alpha_measure, outcome) {
  lm_base <- lmer(formula(paste0(alpha_measure," ~ (1 | family_ID)")), df, REML=F)
  lm_week <- lmer(formula(paste0(alpha_measure," ~  week + (1|family_ID)")), df, REML=F)
  lm_outcome <- lmer(formula(paste0(alpha_measure," ~  week + ", outcome, " + (1|family_ID)")), df, REML=F)
  anova(lm_outcome, lm_week)$Pr[2]
}

# calculate lme models for binary ORV outcomes
longitudinal_results = function(t1, outcome, alpha1, alpha2, with_UK=FALSE) {
  # calculate p values for Shannon
  p_shan = c(lme_alpha_ORV(subset(t1, country=="India"), alpha1, outcome),
             lme_alpha_ORV(subset(t1, country_exposure=="India (exposed)"), alpha1, outcome), 
             lme_alpha_ORV(subset(t1, country_exposure=="India (unexposed)"), alpha1, outcome), 
             lme_alpha_ORV(subset(t1, country=="Malawi"), alpha1, outcome))
  if (with_UK==TRUE) { p_shan = c(p_shan, lme_alpha_ORV(subset(t1, country=="UK"), alpha1, outcome)) }
  
  # calculate p values for richness
  p_obs = c(lme_alpha_ORV(subset(t1, country=="India"), alpha2, outcome), 
            lme_alpha_ORV(subset(t1, country_exposure=="India (exposed)"), alpha2, outcome), 
            lme_alpha_ORV(subset(t1, country_exposure=="India (unexposed)"), alpha2, outcome), 
            lme_alpha_ORV(subset(t1, country=="Malawi"), alpha2, outcome))
  if (with_UK==TRUE) { p_obs = c(p_obs, lme_alpha(subset(t1, country=="UK"), alpha2, outcome)) }
  
  # collate results
  results = data.frame(round(p_shan,4), round(p_obs, 4))
  names(results) = c(alpha1, alpha2)
  names = c("India", "India (exposed)", "India (unexposed)", "Malawi")
  if (with_UK==TRUE) { names = c(names, "UK") }
  rownames(results) =  names
  results
}

# calculate correlation coefficients for country subset
corr_by_country_alpha = function(df, outcome_measure, country_exposure_subset, sample_subset, return_val, rsv=FALSE) {
  s1 = subset(df, sample_type==sample_subset & !is.na(outcome_measure))
  # corr_function will automatically log richness values
  if (rsv==TRUE) {
    res = c(corr_function(s1, outcome_measure, "Shannon_rsv", country_exposure_subset, return_option=return_val),
            corr_function(s1, outcome_measure, "Observed_rsv", country_exposure_subset, return_option=return_val))
  } else {
    res = c(corr_function(s1, outcome_measure, "Shannon", country_exposure_subset, return_option=return_val),
            corr_function(s1, outcome_measure, "Observed", country_exposure_subset, return_option=return_val))
  }
}

# calculate lme models for neonatal exposure status
longitudinal_results_neo = function(t1, outcome, alpha1, alpha2) {
  p_shan = c(lme_alpha(subset(t1, country=="India"), alpha1, outcome))
  p_obs = c(lme_alpha(subset(t1, country=="India"), alpha2, outcome))
  results = data.frame(round(p_shan,4), round(p_obs, 4))
  names(results) = c(alpha1, alpha2)
  names = c("India")
  rownames(results) =  names
  results
}

# function to pick out subset of samples and plot PCoA via chosen method
beta_plot_ORV = function(ps_input, subset, plot_title, country_subset, method=c("Bray-Curtis (weighted)", "Bray-Curtis (unweighted)"), outcome=c("seroconv", "dose1_shedding")) {
  # pick out selected subset of samples
  t = data.frame(sample_data(ps_input))
  t$outcome = t[,outcome]
  if (country_subset=="India (exposed)" | country_subset=="India (unexposed)") { 
    samplelist = subset(t, country_exposure==country_subset & sample_type==subset)$sample_ID_full 
  } else {
    samplelist = subset(t, country==country_subset & sample_type==subset)$sample_ID_full
  }
  ps_input = prune_samples(samplelist, ps_input)
  beta = data.frame(sample_data(ps_input))
  beta$outcome = beta[,outcome]
  
  # calculate beta diversity distances for selected samples
  if (method == "Bray-Curtis (weighted)") { beta_dist = phyloseq::distance(ps_input, method = "bray") }
  if (method == "Bray-Curtis (unweighted)") { 
    otu_table(ps_input)[otu_table(ps_input)<0.001] = 0
    beta_dist = phyloseq::distance(ps_input, method = "bray", binary=TRUE) 
  }
  
  # ordinate and extract PCs
  beta_ord = ordinate(ps_input, method = "PCoA", distance = beta_dist)
  beta$PC1 = beta_ord$vectors[,1]
  beta$PC2 = beta_ord$vectors[,2]
  
  # create PCoA plot
  ggplot(beta, aes(PC1, PC2, colour=factor(outcome))) + geom_point(size = 1.5, alpha = 0.5) +
    theme(legend.position="none", legend.title = element_blank()) + ggtitle(plot_title) +
    scale_color_manual(values=c(nr_col, r_col)) + 
    xlab(paste0("Bray-Curtis, PC1 (", format(round(beta_ord$values[1,1],1), nsmall=1), "%)")) + 
    ylab(paste0("Bray-Curtis, PC2 (", format(round(beta_ord$values[2,1],1), nsmall=1), "%)")) +
    theme(axis.title=element_text(size=10), plot.title=element_text(size=10))
}

# plot summary of cross-validation accuracy
summary_rf_plot_binary = function(collated_df, plot_title) {
  if (any(collated_df$country=="UK")) {
   collated_df$week = revalue(as.factor(collated_df$sample_subset), c("BS1" = "week 1", "BS2" = "week 4", "BS3" = "week 6/8", 
                                                                     "BS5" = "week 10/12", "MS1" = "mother"))
  } else {
    collated_df$week = revalue(as.factor(collated_df$sample_subset), c("BS1" = "week 1", "BS2" = "week 4", "BS3" = "week 6", 
                                                                       "BS5" = "week 10", "MS1" = "mother"))
  }
  ggplot(collated_df, aes(y = crossval_mean, x = factor(week), colour=factor(country))) + 
    geom_point(size=2, alpha = 0.8) + 
    scale_x_discrete(limits = rev(levels(collated_df$week))) + 
    geom_errorbar(aes(ymin=crossval_mean-crossval_sd, ymax=crossval_mean+crossval_sd), width=0) + 
    geom_hline(yintercept=0.5, linetype="dotted") + xlab("") +
    ylab("cross-validation accuracy") +
    scale_colour_manual(values = c(India = India_col, "India (exposed)" = India_col_2, "India (unexposed)" = India_col_2, Malawi = Malawi_col, UK = UK_col)) +
    ylim(0.2,0.8) + theme_bw() + theme(legend.position = "", strip.background = element_blank()) +
    coord_flip() + ggtitle(plot_title) + theme(plot.title = element_text(size=10, hjust = 0.5)) + facet_grid(.~country)
}

# function to perform ZINB models
zinb_ORV_function = function(outcome=c("age", "seroconv", "IgA", "dose1", "either", "any_response", "mode_delivery"), include_nb=TRUE) {
  if (outcome=="age") {
    # ZINB
    f <- mms(y = features, fixed = ~ week + offset(log(N)), random = ~ 1 | subject, 
             min.p = 0, method = "zinb", verbose=TRUE)
    res = data.frame(get.fixed(f, part="dist", vr.name="week", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb <- mms(y = features, fixed = ~ week + offset(log(N)), random = ~ 1 | subject,
                min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="week", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
    
  }
  
  if (outcome=="seroconv") {
    # ZINB
    f <- mms(y = features, fixed = ~ seroconv + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="seroconv", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ seroconv + week + offset(log(N)), random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="seroconv", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
  }
  
  
  if (outcome=="dose1") {
    # ZINB
    f <- mms(y = features, fixed = ~ dose1 + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="dose1", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ dose1 + week + offset(log(N)), 
               random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="dose1", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
  }
  
  if (outcome=="either") {
    # ZINB
    f <- mms(y = features, fixed = ~ either + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="either", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ either + week + offset(log(N)), 
               random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="either", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
  }
  
  if (outcome=="any_response") {
    # ZINB
    f <- mms(y = features, fixed = ~ any_response + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="any_response", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ any_response + week + offset(log(N)), 
               random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="any_response", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
  }
  
  if (outcome=="IgA") {
    # ZINB
    f <- mms(y = features, fixed = ~ IgA + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="IgA", sort.p=F))
    res$method = "zinb"
    
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ IgA + week + offset(log(N)), random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="IgA", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
  }
  
  if (outcome=="mode_delivery") {
    # ZINB
    f <- mms(y = features, fixed = ~ mode_delivery + week + offset(log(N)), 
             random = ~ 1 | subject,
             min.p = 0, method = "zinb", verbose=TRUE)
    res <- data.frame(get.fixed(f, part="dist", vr.name="mode_deliveryvaginal", sort.p=F))
    res$method = "zinb"
    
    if (include_nb==TRUE) {
    # NB (if prevalence >95%)
    f_nb = mms(y = features, fixed = ~ mode_delivery + week + offset(log(N)), random = ~ 1 | subject,
               min.p = 0.95, method = "nb", verbose=TRUE)
    res_nb = data.frame(get.fixed(f_nb, part="dist", vr.name="mode_deliveryvaginal", sort.p=F))
    res_nb$method = "nb"
    for (i in 1:nrow(res_nb)){
      if (rownames(res_nb)[i] %in% rownames(res)) { res[rownames(res_nb)[i],] = res_nb[i,] } else { res = rbind(res, res_nb[i,]) }
    }
    }
  }
  
  res$taxon = rownames(res)
  res$padj = p.adjust(res$pvalue, method="BH")
  res$signif = "ns"
  res$signif[res$padj<0.1 & res$Estimate<0]  = "negatively correlated"
  res$signif[res$padj<0.1 & res$Estimate>0]  = "positively correlated"
  res$pscale = -log10(res$padj)
  res$n_tested = ncol(features)
  res$failed_convergence = ncol(features)-nrow(res)
  res$comparison = outcome
  
  # append prevalence and mean abundance
  res$prevalence = NA
  features1 = features[,names(features) %in% rownames(res)]
  res = res[names(features1),]
  for (i in 1:nrow(res)) { res$prevalence[i] = sum(features1[,i]>0)/nrow(features1) }
  features_rel = features1/rowSums(features1)
  res$mean = colMeans(features_rel)
  return(res)
}

# functin to summarise ZINB results
zinb_summary = function(res, ps_t) {
  res_sig = res[res$signif!="ns" & res$comparison!="age", c("taxon", "Estimate", "Std.Error", "padj", "comparison")]
  rps_t = transform_sample_counts(ps_t, function(x) {x/sum(x)})
  rps_df = data.frame(t(otu_table(rps_t)))
  names(rps_df) = data.frame(as(tax_table(rps_t),"matrix"))$Genus
  
  if (nrow(res_sig)>0) {
    for (i in 1:nrow(res_sig)) {
      res_sig$mean[i] = round(mean(rps_df[,res_sig$taxon[i]]*100),3)
      res_sig$prev[i] = round(sum(rps_df[,res_sig$taxon[i]]>0)/nrow(rps_df)*100,1)
    }
    res_sig
  } else { "No significant differences according to ORV outcome"}
}





### Module 5

# modular RF importance plot
### gini plot for RV1 IgA titre module 4
modules_importance_plot = function(path, binary=TRUE) {
  
  importance = read.csv(file=paste0(path), row.names = 1)
  importance$variable = rownames(importance)
  importance = importance[1:10,]
  
  # assign modules to important variables
  importance$module = NA
  importance$module[importance$variable %in% names(module1)] = "module1"
  importance$module[importance$variable %in% names(module2)] = "module2"
  importance$module[importance$variable %in% names(module3)] = "module3"
  
  # rename variables
  importance$variable = gsub("_", " ", importance$variable) 
  importance$variable = gsub("pre exposed", "RV pre-exposed", importance$variable) 
  importance$variable = gsub("MB1 IgG", "maternal RV-IgG", importance$variable) 
  importance$variable = gsub("CB1 IgG", "cord blood RV-IgG", importance$variable) 
  importance$variable = gsub("AG Concugml", "week 6 α1AG", importance$variable)
  importance$variable = gsub("A1AT BS3 Concugml", "week 6 α1AT", importance$variable)
  importance$variable = gsub("A1AT BS5 Concugml", "week 10 α1AT", importance$variable)
  importance$variable = gsub("MPO BS3 Concngml", "week 6 MPO", importance$variable)
  importance$variable = gsub("MPO BS5 Concngml", "week 10 MPO", importance$variable)
  importance$variable = gsub("breastfed child", "breastfeeding practice", importance$variable) 
  importance$variable = gsub("OPV IPV", "OPV/IPV schedule", importance$variable) 
  importance$variable = gsub("Escherichia.Shigella", "Escherichia", importance$variable)
  importance$variable = gsub("HAZ 6 8w", "week 6 HAZ", importance$variable) 
  importance$variable = gsub("MB1 IgA", "maternal RV-IgA", importance$variable) 
  importance$variable = gsub("BM1 IgA", "breastmilk RV-IgA", importance$variable) 
  importance$variable = gsub("birth weight", "birthweight", importance$variable) 
  importance$module[importance$variable=="Escherichia"] = "module3"
  importance$variable = factor(importance$variable, levels = unique(importance$variable))
  
  cols = brewer.pal("Spectral",n = 4)
  
  if (binary==TRUE) {
    ggplot(importance, aes(x = variable, y = mean, fill=module)) +
      scale_x_discrete(limits = rev(levels(importance$variable))) +
      xlab("") + ylab("importance (Gini)") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + 
      scale_fill_manual(values=c(module1 = cols[4], module2 = cols[3], module3 = cols[2])) +
      theme(plot.title = element_text(size=9, hjust = 0.5))
  } else {
    ggplot(importance, aes(x = variable, y = IncMSE_means, fill=module)) +
      scale_x_discrete(limits = rev(levels(importance$variable))) +
      xlab("") + ylab("increase in MSE") +
      geom_bar(position="dodge", stat = "identity") + coord_flip() + theme_bw() +
      theme(legend.position = "none") + 
      scale_fill_manual(values=c(module1 = cols[4], module2 = cols[3], module3 = cols[2])) +
      theme(plot.title = element_text(size=9, hjust = 0.5))
  }
}

# modular RF outcomes plot by country
modules_fig_main = function(df, title) {
  ggplot(df[1:4,], aes(factor(module), y=crossval_mean, colour=module)) + 
    geom_point(aes(size=nvar), alpha = 0.8) +
    geom_errorbar(aes(ymin=crossval_mean-crossval_sd, ymax=crossval_mean+crossval_sd), width=0) +
    labs(colour = "module", size = "number of variables") +
    ylab("") + xlab("")  + scale_size(limits=c(0,160)) +
    scale_color_manual(values = c("#D7191C", "#FDAE61", "#ABDDA4", "#2B83BA")) +
    theme_bw() + theme(legend.position = "") + ylim(0.3,0.8) +
    coord_flip() + ggtitle(title) + theme(plot.title = element_text(size=10, hjust = 0.5)) 
}

# modular RF ouctomes plot - India exposure subsets
modules_fig_extra = function(df) {
  ggplot(df[5:6,], aes(factor(module), y=crossval_mean, colour=module)) + 
    geom_point(aes(size=nvar), alpha = 0.8) +
    geom_errorbar(aes(ymin=crossval_mean-crossval_sd, ymax=crossval_mean+crossval_sd), width=0) + 
    labs(colour = "module", size = "number of variables") +
    ylab("cross-validation accuracy") + xlab("")  + 
    scale_color_manual(values = rep("#D7191C",2)) + scale_size(limits=c(0,160)) +
    theme_bw() + theme(legend.position = "", strip.background = element_blank()) + ylim(0.3,0.8) +
    coord_flip()
}



# function to carry out logistic regression for binary outcomes
summarise_regression = function(variable, input_data, log=FALSE) {
  input_data$variable = input_data[,variable]
  
  # create subsets of data with complete measurements
  input_sub = subset(input_data, !is.na(variable) & !is.na(outcome))
  nr = subset(input_sub, outcome==0)
  r = subset(input_sub, outcome==1)
  
  # summarise n
  n_summary = paste0(nrow(input_sub)," (",nrow(r),")")
  
  # for continuous variable, summarise means in responders and non-responders
  if (nrow(r)>=10 & nrow(nr)>=10) {
    
    # calculate mean (sd), or geometric mean (95% CI)
    if(class(input_sub$variable)=="integer" | class(input_sub$variable)=="numeric") {
      nr_summary = paste0(format(round(mean(nr$variable),1),nsmall=1), " (",format(round(sd(nr$variable),1),nsmall=1),")")
      if (log==TRUE & variable!="Observed_6w" & variable!="Observed_10w" & variable!="houshold_income") { 
        gsum = Gmean(nr$variable, conf.level =0.95)
        nr_summary = paste0(format(round(gsum[1],1),nsmall=1), 
                            " (",format(round(gsum[2],1),nsmall=1),"-", format(round(gsum[3],1),nsmall=1),")") 
      }
      r_summary = paste0(format(round(mean(r$variable),1),nsmall=1), " (",format(round(sd(r$variable),1),nsmall=1),")")
      if (log==TRUE & variable!="Observed_6w" & variable!="Observed_10w" & variable!="houshold_income") { 
        gsum = Gmean(r$variable, conf.level =0.95)
        r_summary = paste0(format(round(gsum[1],1),nsmall=1), 
                           " (",format(round(gsum[2],1),nsmall=1),"-", format(round(gsum[3],1),nsmall=1),")") 
      }
      
      # run logistic regression
      if (log == TRUE) { input_sub$variable = log(input_sub$variable) }
      lr = glm(outcome ~ variable, data = input_sub, family = binomial)
      rr = odds_to_rr(lr)[2,]
      rr_clean = paste0(format(round(rr[3],2), nsmall=2)," (",format(round(rr[4],2), nsmall=2),"-",format(round(rr[5],2), nsmall=2),")")
      p = format(round(summary(lr)$coefficients[2,4],3),nsmall=3)
      
      # collate results
      c(variable, n_summary, r_summary, nr_summary, rr_clean, p)
     
    # if categorical variable, summarise prevalence in responders and non-responders 
    } else if (class(input_sub$variable)=="factor") {
      levels = levels(input_sub$variable)
      
      if (sum(input_sub$variable==levels[1])>0 & sum(input_sub$variable==levels[2])>0) {
        
        nr_summary = paste0(sum(nr$variable==levels[2]),"/",nrow(nr)," (",format(round(sum(nr$variable==levels[2])/nrow(nr)*100,1),nsmall=1),"%)")
        r_summary = paste0(sum(r$variable==levels[2]),"/",nrow(r)," (",format(round(sum(r$variable==levels[2])/nrow(r)*100,1),nsmall=1),"%)")
        lr = glm(outcome ~ variable, data = input_sub, family = binomial)
        rr = odds_to_rr(lr)[2,]
        rr_clean = paste0(format(round(rr[3],2),nsmall=2)," (",format(round(rr[4],2),nsmall=2),"-",format(round(rr[5],2),nsmall=2),")")
        p = format(round(summary(lr)$coefficients[2,4],3),nsmall=3)
        
        if ((sum(nr$variable==levels[2])<5 & sum(r$variable==levels[2])<5) | 
            (sum(nr$variable==levels[1])<5 & sum(r$variable==levels[1])<5)) {
          c(paste0(variable," (",levels[2],")"), n_summary, r_summary, nr_summary, NA, NA)
        } else {
          c(paste0(variable," (",levels[2],")"), n_summary, r_summary, nr_summary, rr_clean, p)
        }
      } else { c(variable, n_summary, NA, NA, NA, NA) }
    }
  } else { c(variable, n_summary, NA, NA, NA, NA) }
}

# function to collate regression outputs for cofactors of interest 
baseline_summary = function(input_data) {
  summary = data.frame(rbind(
    # infant characteristics
    summarise_regression("sex_baby", input_data),
    summarise_regression("mode_delivery", input_data),
    summarise_regression("place_delivery", input_data),
    summarise_regression("breastfed_child", input_data),
    summarise_regression("age_at_first_dose", input_data),
    summarise_regression("birth_weight", input_data),
    summarise_regression("HIV_status", input_data),
    summarise_regression("HAZ_6_8w", input_data),
    summarise_regression("antibiotic_exposure", input_data),
    
    # maternal characteristics
    summarise_regression("age_mother", input_data),
    summarise_regression("mother_weight", input_data),
    summarise_regression("parity", input_data),
    summarise_regression("maternal_education", input_data),
    
    # household characteristics
    summarise_regression("type_house", input_data),
    summarise_regression("water_treat_wami", input_data),
    summarise_regression("sanitation", input_data),
    summarise_regression("kitchen", input_data),
    summarise_regression("refrigerator", input_data),
    summarise_regression("crowding_wami", input_data),
    summarise_regression("household_income", input_data, log=TRUE),
    
    # Exposure variables
    summarise_regression("A1AT_BS3_Concugml", input_data, log=TRUE),
    summarise_regression("A1AT_BS5_Concugml", input_data, log=TRUE),
    summarise_regression("MPO_BS3_Concngml", input_data, log=TRUE),
    summarise_regression("MPO_BS5_Concngml", input_data, log=TRUE),
    summarise_regression("AG_Concugml", input_data, log=TRUE),
    
    # Maternal antibodies
    summarise_regression("MB1_IgA", input_data, log=TRUE),
    summarise_regression("BM1_IgA", input_data, log=TRUE),
    summarise_regression("MB1_IgG", input_data, log=TRUE),
    summarise_regression("CB1_IgG", input_data, log=TRUE),
    
    # Infant antibodies
    summarise_regression("BB1_IgG", input_data, log=TRUE),
    
    # 16S at baseline
    summarise_regression("Shannon_6w", input_data),
    summarise_regression("Observed_6w", input_data, log=TRUE),
    summarise_regression("Shannon_10w", input_data),
    summarise_regression("Observed_10w", input_data, log=TRUE)
  ), stringsAsFactors = F)
  names(summary) = c("variable", "n total (r)", "responders", "non-responders", "RR (95% CI)", "p")
  if (class(a$outcome)=="numeric") { names(summary) = c("variable", "n", "level 1", "level 2", "pearson", "beta (95% CI)", "p")}
  summary
}
