####################################################################
### Script to generate metadata from input counts and sample IDs ###
####################################################################

# import read count data from summary csv
m = read.csv("input/qiime2_stats_14102020.csv", header=T, stringsAsFactors=F) 

# final count
m$final_count = m$nonchimeric

# sample ID
m$sample_ID = NA
m$sample_ID[grepl("tL", m$sample_ID_full)] = gsub("tL", "L", m$sample_ID_full)[grepl("tL", m$sample_ID_full)]
m$sample_ID[grepl("tp", m$sample_ID_full)] = gsub("tp", "p", m$sample_ID_full)[grepl("tp", m$sample_ID_full)]

# sequence site
m$site = NA
m$site[grepl("LIMS", m$sample_ID_full)] = "Liverpool"
m$site[grepl("Imp", m$sample_ID_full)] = "London"

# run
m$run[m$site=="Liverpool"] = substr(m$sample_ID[m$site=="Liverpool"], 0, 9)
m$run[m$site=="London"] = substr(m$sample_ID[m$site=="London"], 0, 4)
m$run[grepl("LIMS14462l1", m$sample_ID_full)] = "LIMS14462l1" # differentiate between two lanes on LIMS14462
m$run[grepl("LIMS14462l2", m$sample_ID_full)] = "LIMS14462l2" # differentiate between two lanes on LIMS14462

# control
m$control = revalue(as.factor(grepl("ctrl", m$sample_ID) | grepl("carryover", m$sample_ID)), c("FALSE" = "sample", "TRUE" = "ctrl"))

# sample code
m$sample_code = "ctrl"
m$sample_code[m$control=="sample"] = substrRight(as.character(m$sample_ID)[m$control=="sample"], 7)

# sample type
m$sample_type = substrRight(m$sample_code, 3)
m$sample_type[grepl("WCS", m$sample_ID)] = "WCS"
m$sample_type[grepl("MCctrl", m$sample_ID)] = "MCctrl"
m$sample_type[grepl("MCQia", m$sample_ID)] = "MCQia"
m$sample_type[grepl("PCRctrl", m$sample_ID)] = "PCRneg"
m$sample_type[grepl("carryover", m$sample_ID)] = "carryover"
m$sample_type[grepl("BMctrl", m$sample_ID)] = "BMctrl"
m$sample_type[grepl("WCBM", m$sample_ID)] = "WCBM"
m$sample_type[grepl("BSctrl", m$sample_ID)] = "BSctrl"
m$sample_type[grepl("MSctrl", m$sample_ID)] = "MSctrl"
m$sample_type[grepl("MCQiactrlBM", m$sample_ID)] = "MCQiaBM"
m$sample_type[grepl("MCctrlBM", m$sample_ID)] = "MCctrlBM"
m$sample_type[grepl("PCRctrlBM", m$sample_ID)] = "PCRnegBM"

# country
m$country = NA
m$country[grepl("U[0-9][0-9][0-9]", m$sample_code)] = "UK"
m$country[grepl("M[0-9][0-9][0-9]", m$sample_code)] = "Malawi"
m$country[grepl("I[0-9][0-9][0-9]", m$sample_code)] = "India"
m$country[m$control=="ctrl"] = "ctrl"

# family
m$family_ID = NA
m$family_ID[m$control=="sample"] = substr(as.character(m$sample_code)[m$control=="sample"], 0, 4)

# week as continuous variable
m$week = NA
m$week[m$sample_type=="BS1" | m$sample_type=="MS1" | m$sample_type=="BM1"] = 1
m$week[m$sample_type=="BS2"] = 4
m$week[m$sample_type=="BS3" & m$country!="UK"] = 6
m$week[m$sample_type=="BS3" & m$country=="UK"] = 8
m$week[m$sample_type=="BS4" & m$country!="UK"] = 7
m$week[m$sample_type=="BS4" & m$country=="UK"] = 9
m$week[m$sample_type=="BS5" & m$country!="UK"] = 10
m$week[m$sample_type=="BS5" & m$country=="UK"] = 12
m$week[m$sample_type=="BS6" & m$country!="UK"] = 11
m$week[m$sample_type=="BS6" & m$country=="UK"] = 13
m$week[m$sample_type=="BM2" & m$country!="UK"] = 7
m$week[m$sample_type=="BM2" & m$country=="UK"] = 9
m$week[m$sample_type=="BM3" & m$country!="UK"] = 11
m$week[m$sample_type=="BM3" & m$country=="UK"] = 13

# week for combined figure plots
m$week_fig = NA
m$week_fig[m$sample_type=="BM1"] = "week 1"
m$week_fig[m$sample_type=="BM2"] = "week 7/9"
m$week_fig[m$sample_type=="BM3"] = "week 11/13"
m$week_fig[m$sample_type=="BS1"] = "week 1"
m$week_fig[m$sample_type=="BS2"] = "week 4"
m$week_fig[m$sample_type=="BS3"] = "week 6/8"
m$week_fig[m$sample_type=="BS5"] = "week 10/12"
m$week_fig[m$sample_type=="MS1"] = "mother"

# sample group
m$sample_group = NA
m$sample_group[m$control=="ctrl"] = "ctrl"
m$sample_group[m$control=="sample" & grepl("MS", m$sample_type)] = "MS"
m$sample_group[m$control=="sample" & grepl("BM", m$sample_type)] = "BM"
m$sample_group[m$control=="sample" & grepl("BS", m$sample_type)] = "BS"

# sample group fig
m$sample_group_fig = m$sample_group
m$sample_group_fig[m$sample_group=="MS"] = "maternal stool"
m$sample_group_fig[m$sample_group=="BS"] = "infant stool"
m$sample_group_fig[m$sample_group=="BM"] = "breastmilk"

# sample number
m$sample_number = m$sample_ID
m$sample_number = gsub("^.*?s","s", m$sample_number)
m$sample_number = gsub("([A-Z]).*","", m$sample_number)
m$sample_number = gsub("([c]).*","", m$sample_number)
m$sample_number = gsub("no","", m$sample_number)

# append nanodrop data to sample data
nano = read.csv("input/nanodrop.csv", stringsAsFactors = F)
m = merge(m, nano, by="sample_ID_full", all.x=TRUE)

# for WCs in nanodrop dataframe, append country
m$country[!is.na(m$wc_country)] = m$wc_country[!is.na(m$wc_country)]

# add logical variable to specify extraction and PCR controls
m$is.neg = m$sample_type=="PCRneg" | m$sample_type=="PCRnegBM" | m$sample_type=="WCBM" | m$sample_type=="WCBM"

# append full metadata
s = read.csv("input/RoVI_metadata_by_family_03112020.csv", header=TRUE, stringsAsFactors=FALSE)
s = dplyr::select(s, -c(country))
m = merge(m, s, by = "family_ID", all.x=TRUE)

# exclude samples with 0 input sequences
metadata = m[m$input>0,]
write.csv(m, "input/RoVI_metadata_by_sample.csv")

# remove data
rm(m,s,nano)