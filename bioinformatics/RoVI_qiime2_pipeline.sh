#!/bin/sh

########################################################
###################### Preperations ####################
########################################################

# 1. Download fastq files and ensure names are in the format sampleID_barcodes_L00*_R*_001.fastq
# 2. In home directory, activate qiime2
# 3. Count processors using 'grep -c processor /proc/cpuinfo' (on linux) or 'sysctl -n hw.ncpu' (on unix), and update THREADS accordingly
THREADS=60
# 4. Initiate script: nohup bash RoVI_qiime2_pipeline.sh > RoVI_qiime2_output.txt >&1 &

########################################################
############# Trim primers with cutadapt ###############
########################################################
# count files: ls -1 *_001.fastq.gz | wc -l
# create list of sample names
ls *_R1* > sample_list.txt
sed -i '' -e 's/_R1_001.fastq.gz//' sample_list.txt

mkdir primer_trimmed
parallel -j ${THREADS} -a sample_list.txt "cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTNTCTAATCC -o primer_trimmed/t{}_R1_001.fastq.gz {}_R1_001.fastq.gz -p primer_trimmed/t{}_R2_001.fastq.gz {}_R2_001.fastq.gz --discard-untrimmed -O 10 -e 0.2"

# move input files into separate folder
mkdir input_fastq
mv *.fastq.gz input_fastq

# move files into run-specific folders
run_list=( "LIMS12416" "LIMS12651" "LIMS14462" "LIMS14801" "LIMS15089" "LIMS15168" "LIMS15350" "LIMS15914" "LIMS15990" "LIMS16518" "LIMS16519" "LIMS17407" "LIMS18668" "LIMS18669" "p1p2" "p3p4")
for run_ID in  ${run_list[@]}
    do echo $run_ID;
    mkdir ${run_ID}; mv primer_trimmed/*${run_ID}* ${run_ID}
done
rm -r primer_trimmed

########################################################
############## Run Dada2 for each run ##################
########################################################
mkdir dada2

### run dada2 for single run
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path LIMS12416 --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path dada2/LIMS12416_input.qza
#qiime demux summarize --i-data dada2/LIMS12416_input.qza --o-visualization dada2/LIMS12416_summary.qzv
#qiime tools view dada2/LIMS12416_summary.qzv
#qiime dada2 denoise-paired --i-demultiplexed-seqs dada2/LIMS12416_input.qza --p-trunc-len-f 270 --p-trunc-len-r 200 --o-table dada2/LIMS12416_table.qza --o-representative-sequences dada2/LIMS12416_repseqs.qza --o-denoising-stats dada2/LIMS12416_denoisingstats.qza --p-n-threads ${THREADS}

# set temporary directory for this session
mkdir qiime2-tmp/
export TMPDIR="$PWD/qiime2-tmp/"

# run dada2 in loop across runs
for run_ID in  ${run_list[@]}
    do echo "Processing $run_ID";
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ${run_ID} --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path dada2/${run_ID}_input.qza
    qiime demux summarize --i-data dada2/${run_ID}_input.qza --o-visualization dada2/${run_ID}_summary.qzv
    qiime dada2 denoise-paired --i-demultiplexed-seqs dada2/${run_ID}_input.qza --p-trunc-len-f 270 --p-trunc-len-r 200 --o-table dada2/${run_ID}_table.qza --o-representative-sequences dada2/${run_ID}_repseqs.qza --o-denoising-stats dada2/${run_ID}_denoisingstats.qza --p-n-threads ${THREADS}
    qiime metadata tabulate --m-input-file dada2/${run_ID}_denoisingstats.qza --o-visualization dada2/${run_ID}_denoisingstats.qzv
    qiime feature-table summarize --i-table dada2/${run_ID}_table.qza --o-visualization dada2/${run_ID}_table.qzv
done

########################################################
################ Merge and summarise ###################
########################################################
qiime feature-table merge --i-tables dada2/LIMS12416_table.qza --i-tables dada2/LIMS12651_table.qza --i-tables dada2/LIMS14462_table.qza \
--i-tables dada2/LIMS14801_table.qza --i-tables dada2/LIMS15089_table.qza --i-tables dada2/LIMS15168_table.qza \
--i-tables dada2/LIMS15350_table.qza --i-tables dada2/LIMS15914_table.qza --i-tables dada2/LIMS15990_table.qza \
--i-tables dada2/LIMS16518_table.qza --i-tables dada2/LIMS16519_table.qza --i-tables dada2/LIMS17407_table.qza \
--i-tables dada2/LIMS18668_table.qza --i-tables dada2/LIMS18669_table.qza --i-tables dada2/p1p2_table.qza \
--i-tables dada2/p3p4_table.qza --o-merged-table dada2/merged_table.qza

qiime feature-table merge-seqs --i-data dada2/LIMS12416_repseqs.qza --i-data dada2/LIMS12651_repseqs.qza --i-data dada2/LIMS14462_repseqs.qza \
--i-data dada2/LIMS14801_repseqs.qza --i-data dada2/LIMS15089_repseqs.qza --i-data dada2/LIMS15168_repseqs.qza \
--i-data dada2/LIMS15350_repseqs.qza --i-data dada2/LIMS15914_repseqs.qza --i-data dada2/LIMS15990_repseqs.qza \
--i-data dada2/LIMS16518_repseqs.qza --i-data dada2/LIMS16519_repseqs.qza --i-data dada2/LIMS17407_repseqs.qza  \
--i-data dada2/LIMS18668_repseqs.qza --i-data dada2/LIMS18669_repseqs.qza --i-data dada2/p1p2_repseqs.qza \
--i-data dada2/p3p4_repseqs.qza \
--o-merged-data dada2/merged_repseqs.qza

qiime feature-table summarize --i-table dada2/merged_table.qza  --o-visualization dada2/merged_table.qzv
qiime feature-table tabulate-seqs --i-data dada2/merged_repseqs.qza --o-visualization dada2/merged_repseqs.qzv

# export documents
qiime tools export --input-path dada2/merged_table.qza --output-path dada2
mv dada2/feature-table.biom dada2/merged_table.biom
qiime tools export --input-path dada2/merged_repseqs.qza --output-path dada2
mv dada2/dna-sequences.fasta dada2/merged_repseqs.fasta
biom convert -i dada2/merged_table.biom -o dada2/merged_table.txt --to-tsv
cp dada2/merged_table.txt dada2/merged_table_clean.txt
sed -i '1d' dada2/merged_table_clean.txt
sed -i 's/#OTU ID/OTU_ID/g' dada2/merged_table_clean.txt

# deterime run stats
mkdir dada2_stats
for run_ID in  ${run_list[@]}
    do qiime tools export --input-path dada2/${run_ID}_denoisingstats.qza --output-path dada2_stats/
    mv dada2_stats/stats.tsv dada2_stats/${run_ID}_stats.tsv
done
cat dada2_stats/*.tsv > dada2_stats/merged_stats.tsv

########################################################
########### Taxonomy assignment in dada2 ###############
########################################################
Rscript RoVI_taxonomy_assignment.R > RoVI_taxonomy_assignment_output.txt >&1 &

########################################################
################### End of pipeline ####################
########################################################
