# all path variables need to be appended with '/' at the end
# arg1=$1

# the vcf_file, ensure that the two parents are located in the column 10 and 11
vcf_file="/mnt/ibm4/cxn/RNAediting/mutation/filter_snp/vcf_file/groupHC-dim2TAAA4_dim2TAAa6-20220902.vcf"

# work directory
workdir=/mnt/ibm4/cxn/RNAediting/mutation/filter_snp/dim2TAA/
# ibm4workdir=/home/cxn/Data/RNAediting/mutation/filter_snp/dim2TAA/


# BAM files were used to obtain accurate read counts at each position.
bam_dir=/mnt/san3/usr/cxn/bam/dim2TAA-BB/

# used for IGV visualization
bam_path="/mnt/ibm4/cxn/RNAediting/mutation/filter_snp/bam_path.csv"

# the reference fasta file
reference="/mnt/san1/usr/thw/info/neurospora/NC12v42contig20MtHph/NC12v42contig20MtHph.fa"

# indicates whether samples are grouped as tetrads (1) or are random spores (0)
is_cross=1

# group name
filename=$(basename $vcf_file)
cross=${filename%%.*}

# number of parallel threads used when extracting readcounts
parallel_num=10

######################################################################################################
mkdir -p ${workdir}
mkdir -p ${workdir}readcounts
cd ${workdir}

# split mult-allelic SNV sites and split into snps and indels
bcftools norm --fasta-ref ${reference} --keep-sum AD -m - ${vcf_file} | \
    bcftools view -v snps -Oz > ${workdir}${cross}.snp.norm.vcf.gz && tabix -p vcf ${workdir}${cross}.snp.norm.vcf.gz

# extract candidate positions for snps 
vcf-annotate --fill-type ${vcf_file}| \
    perl -ne 'next if(/\#/); next unless(/snp/); my ($chrom, $pos)=(split /\s+/)[0,1]; print "$chrom\t",($pos-1),"\t$pos\n";' \
    > ${workdir}${cross}.snp.bed

# Count accurate allele depths for each locus and each sample
# for snps
mkdir -p ${workdir}readcounts/MQ20NAR/  ${workdir}readcounts/MQ0AR/  ${workdir}readcounts/MQ20AR/
find -L ${bam_dir} -name "*.realn.bam" -print | \
        sed 's/.realn.bam$//' | xargs -n 1 -P ${parallel_num} -I PREFIX \
        sh -c '
            sample=`basename PREFIX | cut -d"." -f1`
            
            echo "[`date`]: Start processing ${sample} ... "
            
            # only use proper pairs, mapping quality >= 20
            samtools mpileup -d100000 -x -q20 -f "$1" -l "$3".snp.bed PREFIX.realn.bam \
                > "$2"readcounts/MQ20NAR/"$3".snp.${sample}.MQ20.NAR.mpileup
            
            varscan readcounts "$2"readcounts/MQ20NAR/"$3".snp.${sample}.MQ20.NAR.mpileup \
                --min-base-qual 20 --min-coverage 1 --output-file "$2"/readcounts/MQ20NAR/"$3".snp.${sample}.MQ20.NAR.readcounts
            
            # count in anomalous reads, mapping quality >= 0
            samtools mpileup -Ad100000 -x -q0 -f "$1" -l "$3".snp.bed PREFIX.realn.bam \
                > "$2"readcounts/MQ0AR/"$3".snp.${sample}.MQ0.AR.mpileup
            
            varscan readcounts "$2"readcounts/MQ0AR/"$3".snp.${sample}.MQ0.AR.mpileup \
                --min-base-qual 20 --min-coverage 1 --output-file "$2"/readcounts/MQ0AR/"$3".snp.${sample}.MQ0.AR.readcounts
            
            # count in anomalous reads, mapping quality >= 20
            samtools mpileup -Ad100000 -x -q20 -f "$1" -l "$3".snp.bed PREFIX.realn.bam \
                > "$2"readcounts/MQ20AR/"$3".snp.${sample}.MQ20.AR.mpileup
            
            varscan readcounts "$2"readcounts/MQ20AR/"$3".snp.${sample}.MQ20.AR.mpileup \
                --min-base-qual 20 --min-coverage 1 --output-file "$2"/readcounts/MQ20AR/"$3".snp.${sample}.MQ20.AR.readcounts  
            echo "[`date`]: Finished processing ${sample}"
        ' _ ${reference} ${workdir} ${cross}


# get the file list of all counting results
for f in $(find ${workdir}readcounts/MQ20NAR/ -name "${cross}.snp.*.MQ20.NAR.readcounts" | sort);
 do
      library=`basename $f | cut -d"." -f3`
      sample=${library}
      echo "${sample} ${library} ${f}"
 done > ${workdir}${cross}.MQ20.NAR.snp.readcounts.list
 
for f in $(find ${workdir}readcounts/MQ0AR/ -name "${cross}.snp.*.MQ0.AR.readcounts" | sort);
 do
      library=`basename $f | cut -d"." -f3`
      sample=${library}
      echo "${sample} ${library} ${f}"
 done > ${workdir}${cross}.MQ0.AR.snp.readcounts.list
 
for f in $(find ${workdir}readcounts/MQ20AR/ -name "${cross}.snp.*.MQ20.AR.readcounts" | sort);
 do
      library=`basename $f | cut -d"." -f3`
      sample=${library}
      echo "${sample} ${library} ${f}"
 done > ${workdir}${cross}.MQ20.AR.snp.readcounts.list
 
 
# fill depth fields (RC) use fillVcfDepth.pl #--minimun-vcf
python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/fill_vcf_by_RC.py \
  --vcf  ${workdir}${cross}.snp.norm.vcf.gz --list ${workdir}${cross}.MQ20.NAR.snp.readcounts.list \
  --outfile ${workdir}${cross}.snp.MQ20NAR_RC.vcf

python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/fill_vcf_by_RC.py \
  --vcf  ${workdir}${cross}.snp.norm.vcf.gz --list ${workdir}${cross}.MQ0.AR.snp.readcounts.list \
  --outfile ${workdir}${cross}.snp.MQ0AR_RC.vcf

python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/fill_vcf_by_RC.py \
  --vcf  ${workdir}${cross}.snp.norm.vcf.gz --list ${workdir}${cross}.MQ20.AR.snp.readcounts.list \
  --outfile ${workdir}${cross}.snp.MQ20AR_RC.vcf
  
#bgzip -c > ${cross}.snp.norm_RC.vcf.gz && tabix -p vcf ${cross}.snp.norm_RC.vcf.gz


#######################################################################################################
#######################################################################################################
# for single spore
if  [ ${is_cross} == 0 ];then
    for f in ${workdir}${cross}.snp.*_RC.vcf;do
      #${workdir}${cross}.snp.MQ20AR_RC.vcf.filter.csv
      python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/filter_spore_mut_v8.py \
          --vcf  ${f} \
          --outpath ${workdir} \
          --case case1 \
          --min_mut_ratio 0.3 --max_parent_mut_ratio 0.2 --max_parent_mut_reads 5 --min_sup_reads 5  --min_sup_strand_reads 2 \
          --cross 0 
    done
    
    for f in ${workdir}${cross}.snp.*_RC.vcf.filter.csv; do
      name=$(basename "$f" | cut -d'.' -f3)
      # È¥µôÄ©Î²µÄ "_RC"
      name=${name%_RC}
      echo -e "${f}\t${name}"
    done > ${workdir}${cross}.snp.RC.filter.list
    
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/merge_RCsource.py \
        --list ${workdir}${cross}.snp.RC.filter.list \
        --outfile ${workdir}${cross}.snp.merge_RC.vcf.filter.csv
        
    mkdir -p ${workdir}igv
    
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/mut_snapshot.py \
        --mutcsv ${workdir}${cross}.snp.merge_RC.vcf.filter.csv \
        --outpath ${workdir}igv/ \
        --bampath "/mnt/ibm4/cxn/RNAediting/mutation/filter_snp/bam_path.csv" \
        --cross 0 
        
        
# can add filtering condition: remove SNPs near indels
############################################    
# for cross
if  [ ${is_cross} == 1 ];then
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/vcf_to_ABregion.py \
        --vcf ${vcf_file} --out_scvcf ${workdir}${cross}.sctag.vcf \
        --outfile ${workdir}${cross}.parentsAB.csv
        
    python vcf_to_fasta.py --vcf  ${workdir}${cross}.sctag.vcf --outpath  ${workdir}
    
    #  plot the backgroud blocks for each spore
    awk '!/\#/' ${workdir}${cross}.parentsAB.csv | \
        cut -f 1 | uniq \
        > ${workdir}${cross}.spores.list
        
    awk 'BEGIN{OFS="\t";} {print $1,$3,$2,$4+1,$5;}' \
         ${workdir}${cross}.parentsAB.csv | \
        perl /opt/nfs/scripts/block2draw.pl -b - -t 4 -d 3 \
        -s ${workdir}${cross}.spores.list \
        -p ${workdir}${cross}.parentsAB.draw.pos.csv \
        > ${workdir}${cross}.parentsAB.draw.csv
    
    Rscript /opt/nfs/scripts/plot.blocks-multi_x.R \
        ${workdir}${cross}.parentsAB.draw.csv \
        ${workdir}${cross}.parentsAB.draw.pos.csv \
        A,B 16,140 \
        ${workdir}${cross}.parentsAB.draw.pdf
    
    # sample tetrad source inferred from parentsAB.draw.pdf (and sctag.tetra), the tetra_info file format is as follows:
    # cross_name  tetrad1  tetrad2  tetrad3  tetrad4, if a tetrad consists of two samples, they are separated by a comma
    # eg. dim2TAA_c09	dim2TAA_c09_1	dim2TAA_c09_3,dim2TAA_c09_7	dim2TAA_c09_4,dim2TAA_c09_5	dim2TAA_c09_6
    tetra_info="/mnt/ibm4/cxn/RNAediting/mutation/filter_snp/tetra_info/groupHC-dim2TAAA4_dim2TAAa6-20220902.sctag.tetra"
        
    for f in ${workdir}${cross}.snp.*_RC.vcf;do
    
      #${workdir}${cross}.snp.MQ20AR_RC.vcf.filter.csv
      #${workdir}${cross}.snp.MQ20AR_RC.vcf.case4.csv
      python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/filter_spore_mut_v8.py \
          --vcf  ${f} \
          --outpath ${workdir} \
          --case case1 \
          --min_sup_ratio 0.1 --min_parent_sup_reads 6 --min_progeny_sup_reads 6 --min_sup_strand_reads 1 \
          --max_parent_mut_reads 5  \
          --cross 1 --tetrainfo ${tetra_info}
    done
    
    for f in ${workdir}${cross}.snp.*_RC.vcf.filter.csv; do
      name=$(basename "$f" | cut -d'.' -f3)
      # del "_RC"
      name=${name%_RC}
      echo -e "${f}\t${name}"
    done > ${workdir}${cross}.snp.RC.filter.list
        
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/merge_RCsource.py \
        --list ${workdir}${cross}.snp.RC.filter.list \
        --outfile ${workdir}${cross}.snp.merge_RC.vcf.filter.csv
        
    # merge case4 site
    awk 'FNR==1 && NR!=1 {next} !seen[$0]++' ${workdir}${cross}.snp.*.vcf.case4.csv > ${workdir}${cross}.snp.merge_RC.vcf.case4.csv

    # modify sctag (del case4)
    awk 'BEGIN{FS=OFS="\t"} $0 !~ /^#/ {print $1, $2-1, $2}' ${workdir}${cross}.snp.merge_RC.vcf.case4.csv | tail -n +2 > ${workdir}${cross}.snp.merge_RC.case4_sites.bed
    
    
    # del case4 sites
    bcftools view -T ^${workdir}${cross}.snp.merge_RC.case4_sites.bed  ${workdir}${cross}.sctag.vcf -Ov -o ${workdir}${cross}.sctag.delcase4.vcf 
    
    
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/vcf_to_ABregion.py \
        --scvcf ${workdir}${cross}.sctag.delcase4.vcf \
        --outfile ${workdir}${cross}.delcase4.parentsAB.csv


    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/spore_to_tetra_mut.py \
        --mutcsv ${workdir}${cross}.snp.merge_RC.vcf.filter.csv \
        --parentinfo ${workdir}${cross}.delcase4.parentsAB.csv --tetrainfo ${tetra_info} \
        --case4 ${workdir}${cross}.snp.merge_RC.vcf.case4.csv --outpath  ${workdir} \
    # ${workdir}${cross}.snp.merge_RC.vcf.filter.cross.csv
    
    
    # for IGV 
    for category in uncertain_case4 uncertain_case5_22 uncertain_case5_n22 uncertain_prob_22 uncertain_prob_n22 uncertain_improb_22 uncertain_improb_n22 reliable_22 reliable_n22 unreliable_22 unreliable_n22 ; do
    echo "${category}"
    mkdir -p ${workdir}igv/${category}
      python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/mut_snapshot.py \
          --mutcsv ${workdir}${cross}.snp.merge_RC.vcf.filter.cross.csv \
          --outpath ${workdir}igv/${category}/ \
          --snap_outpath ${ibm4workdir}igv/${category}/ \
          --bampath ${bam_path} \
          --category ${category} \
          --cross 1
    done
fi
# manually check mutations

############################################    
# for single spore
if  [ ${is_cross} == 0 ];then
    for f in ${workdir}${cross}.snp.*_RC.vcf;do
      #${workdir}${cross}.snp.MQ20AR_RC.vcf.filter.csv
      python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/filter_spore_mut_v8.py \
          --vcf  ${f} \
          --outpath ${workdir} \
          --case case1 case5 \
          --min_sup_ratio 0.1 --min_parent_sup_reads 6 --min_progeny_sup_reads 6 --min_sup_strand_reads 1 \
          --max_parent_mut_reads 5  \
          --cross 0 
    done
    
    for f in ${workdir}${cross}.snp.*_RC.vcf.filter.csv; do
      name=$(basename "$f" | cut -d'.' -f3)
      # del "_RC"
      name=${name%_RC}
      echo -e "${f}\t${name}"
    done > ${workdir}${cross}.snp.RC.filter.list
    
    python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/merge_RCsource.py \
        --list ${workdir}${cross}.snp.RC.filter.list \
        --outfile ${workdir}${cross}.snp.merge_RC.vcf.filter.csv
    
    # for IGV 
    for category in uncertain_case5 uncertain_hr uncertain_lr uncertain_ph reliable unreliable ; do
    echo "${category}"
    mkdir -p ${workdir}igv/${category}
      python /mnt/ibm4/cxn/RNAediting/mutation/filter_snp/mut_snapshot.py \
          --mutcsv ${workdir}${cross}.snp.merge_RC.vcf.filter.csv \
          --outpath ${workdir}igv/${category}/ \
          --snap_outpath ${ibm4workdir}igv/${category}/ \
          --bampath ${bam_path} \
          --category ${category} \
          --cross 0
    done
fi
# manually check mutations

