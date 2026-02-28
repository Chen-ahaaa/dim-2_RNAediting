#set -euo pipefail
fq_path=/mnt/san3/usr/cxn/rawdata/Nc_SRP096774/fq
librarylayout=paired   #single or paired
out_dir=/mnt/san3/usr/cxn/RNAediting/trans/N.crassa
reference=/mnt/san3/usr/cxn/info/neurospora_genome/thwNC12v42-contig20Mt.fa
genome_index=/mnt/san3/usr/cxn/info/neurospora_genome/thwNC12v42-contig20Mt
gff_file=/mnt/san3/usr/cxn/info/neurospora_genome/thwNC12v42.gff
species=N.crass
parallel_num=2
threads=10
gene_region=/mnt/san3/usr/cxn/RNAediting/trans/N.crassa/gene_region
DNA_dir=/mnt/san3/usr/cxn/RNAediting/trans/N.crassa/DNA
###############################################################################################################
# prepare
# build hisat2 index
#hisat2-build ${reference} ${genome_index}

# use RSeQC infer_experiment.py to determine if it is strand-specific
# infer_experiment.py -r gff.bed -i sample.bam 
# F R FR RF none
paired_mode="RF"
###############################################################################################################

mkdir -p ${fq_path}/trimming

if [ ${librarylayout} == "single" ]; then
  find ${fq_path} -maxdepth 1 -type f -name "*.fastq.gz" | \
  xargs -n 1 -P ${parallel_num} -I {} bash -c '
      sample=$(basename {})
      name=${sample%%.fastq.gz}

      # Clean FASTQ reads using Trimmomatic
      java -jar /mnt/san3/usr/cxn/soft/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads "$4" \
          -threads "$4" \
          "$2"/${name}.fastq.gz \
          "$2"/trimming/${name}.fq.gz \
          ILLUMINACLIP:/mnt/san3/usr/cxn/soft/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2:True \
          LEADING:3 TRAILING:3 MINLEN:36 && \
      echo "finished trimming ${name}"
      
      if [ "$7" == "none" ]; then
        # hisat2 alignment
        hisat2 -p  "$4" --dta -x "$5" \
            -U "$2"/trimming/${name}.fq.gz \
            -S "$6"/${name}.sam \
            --new-summary 1> "$6"/${name}_mapping.log 2>&1 
      else
        hisat2 -p  "$4" --dta -x "$5" \
            --rna-strandness "$7" \
            -U "$2"/trimming/${name}.fq.gz \
            -S "$6"/${name}.sam \
            --new-summary 1> "$6"/${name}_mapping.log 2>&1 
      fi
              
      # samtools sort + index
      samtools sort -@ "$4" -O BAM -o "$6"/${name}.sort.bam "$6"/${name}.sam
      
      # remove duplicates
      java -jar /mnt/san1/usr/thw/tansoft/picard-tools-1.124/picard.jar MarkDuplicates \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT="$6"/${name}.sort.bam \
            OUTPUT="$6"/${name}.sort.dedup.bam \
            METRICS_FILE="$6"/${name}.sort.dedup.metrics \
            >> "$6"/${name}_mapping.log 2>&1 && \
            
      samtools index "$6"/${name}.sort.dedup.bam 
      rm "$6"/${name}.sam "$6"/${name}.sort.bam "$6"/${name}.sort.dedup.metrics 
      
      # divide bam file into two separate files containing sense-strand and antisense-strand read alignments
      if [ "$7" == "R" ];then
        samtools view -b -F 16 -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.antisense.bam 
        samtools view -b -f 16 -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.sense.bam 
        samtools index "$6"/${name}.sort.dedup.sense.bam
        samtools index "$6"/${name}.sort.dedup.antisense.bam
        
      elif [ "$7" == "F" ];then
        samtools view -b -F 16 -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.sense.bam 
        samtools view -b -f 16 -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.antisense.bam 
        samtools index "$6"/${name}.sort.dedup.sense.bam
        samtools index "$6"/${name}.sort.dedup.antisense.bam
        
      elif [ "$7" == "none" ];then
        samtools view -b -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.unstrand_sense.bam
        samtools index "$6"/${name}.sort.dedup.unstrand_sense.bam
        
      fi

  ' _ ${librarylayout} ${fq_path} ${parallel_num} ${threads} ${genome_index} ${out_dir} ${paired_mode}


elif [ ${librarylayout} == "paired" ];then
  find ${fq_path} -maxdepth 1 -type f -name "*_1.fastq.gz" | \
  xargs -n 1 -P ${parallel_num} -I {} bash -c '
      sample=$(basename {})
      name=${sample%%_1.fastq.gz}
      
      # Clean FASTQ reads using Trimmomatic
      java -jar /mnt/san3/usr/cxn/soft/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads "$4" \
                -threads "$4" \
                "$2"/${name}_1.fastq.gz "$2"/${name}_2.fastq.gz \
                "$2"/trimming/${name}_1.fq.gz "$2"/trimming/${name}_1unpair.fq.gz \
                "$2"/trimming/${name}_2.fq.gz "$2"/trimming/${name}_2unpair.fq.gz \
                ILLUMINACLIP:/mnt/san3/usr/cxn/soft/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True \
                LEADING:3 TRAILING:3 MINLEN:36 && \
            echo "finished trimming ${name}"
      
      if [ "$7" == "none" ]; then
      # hisat2 alignment
        hisat2 -p  "$4" --dta -x "$5" \
            -1 "$2"/trimming/${name}_1.fq.gz \
            -2 "$2"/trimming/${name}_2.fq.gz \
            -S "$6"/${name}.sam \
            --new-summary 1> "$6"/${name}_mapping.log 2>&1
      else
        hisat2 -p  "$4" --dta -x "$5" \
            --rna-strandness "$7" \
            -1 "$2"/trimming/${name}_1.fq.gz \
            -2 "$2"/trimming/${name}_2.fq.gz \
            -S "$6"/${name}.sam \
            --new-summary 1> "$6"/${name}_mapping.log 2>&1
      fi   
          
      # samtools sort + index
      samtools sort -@ "$4" -O BAM -o "$6"/${name}.sort.bam "$6"/${name}.sam
      
      # remove duplicates
      java -jar /mnt/san1/usr/thw/tansoft/picard-tools-1.124/picard.jar MarkDuplicates \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT="$6"/${name}.sort.bam \
            OUTPUT="$6"/${name}.sort.dedup.bam \
            METRICS_FILE="$6"/${name}.sort.dedup.metrics \
            >> "$6"/${name}_mapping.log 2>&1 && \
            
      samtools index "$6"/${name}.sort.dedup.bam 
      
      rm "$6"/${name}.sam "$6"/${name}.sort.bam "$6"/${name}.sort.dedup.metrics 
      
      # divide bam file into two separate files containing sense-strand and antisense-strand read alignments
      if [ "$7" == "RF" ];then
        samtools merge -f "$6"/${name}.sort.dedup.sense.bam \
          <(samtools view -b -f 163 -F 1024 "$6"/${name}.sort.dedup.bam) \
          <(samtools view -b -f 83  -F 1024 "$6"/${name}.sort.dedup.bam)
          
        samtools merge -f "$6"/${name}.sort.dedup.antisense.bam \
          <(samtools view -b -f 147 -F 1024 "$6"/${name}.sort.dedup.bam) \
          <(samtools view -b -f 99  -F 1024 "$6"/${name}.sort.dedup.bam)
        samtools index "$6"/${name}.sort.dedup.sense.bam
        samtools index "$6"/${name}.sort.dedup.antisense.bam
      elif [ "$7" == "FR" ];then
        samtools merge -f "$6"/${name}.sort.dedup.antisense.bam \
          <(samtools view -b -f 163 -F 1024 "$6"/${name}.sort.dedup.bam) \
          <(samtools view -b -f 83  -F 1024 "$6"/${name}.sort.dedup.bam)
          
        samtools merge -f "$6"/${name}.sort.dedup.sense.bam \
          <(samtools view -b -f 147 -F 1024 "$6"/${name}.sort.dedup.bam) \
          <(samtools view -b -f 99  -F 1024 "$6"/${name}.sort.dedup.bam)
        samtools index "$6"/${name}.sort.dedup.sense.bam
        samtools index "$6"/${name}.sort.dedup.antisense.bam
      elif [ "$7" == "none" ];then
        samtools view -b -F 1024 "$6"/${name}.sort.dedup.bam -o "$6"/${name}.sort.dedup.unstrand_sense.bam
        samtools index "$6"/${name}.sort.dedup.unstrand_sense.bam
      fi

  ' _ ${librarylayout} ${fq_path} ${parallel_num} ${threads} ${genome_index} ${out_dir} ${paired_mode}

fi

################################################################################################################
# detect A-to-I RNA editing events separately from the sense-strand and antisense-strand files using REDItoolDnaRna.py
# prepare DNA_bam 
# if not add DNA_bam,need to check no mutation in parents
# focus only on specific genes, gene_region is a bed file

mkdir -p ${out_dir}/region_bam
mkdir -p ${out_dir}/output
mkdir -p ${out_dir}/region_bam/DNA

if compgen -G ${DNA_dir}/*bam > /dev/null; then
  for f in ${DNA_dir}/*bam;do
    f_name=$(basename ${f})
    samtools view -hb ${f} -L ${gene_region} > ${out_dir}/region_bam/DNA/${f_name}
  done
fi

if [ -s "${gene_region}" ]; then
    echo "[`date "+%F %T"`] Gene region exists, extracting RNA BAM by region"
    RNA_input_dir="${out_dir}/region_bam"
    DNA_input_dir="${out_dir}/region_bam/DNA"
    find ${out_dir} -maxdepth 1 -type f -name "*sense.bam" | \
    xargs -n 1 -P ${parallel_num} -I {} sh -c '
      bam={}
      bam_name=$(basename ${bam})
      samtools view -hb ${bam} -L "$2" -o "$1"/region_bam/${bam_name}
    ' _ ${out_dir} ${gene_region}
    wait

else
    echo "[`date "+%F %T"`] Gene region not found, using full BAMs as input"
    RNA_input_dir="${out_dir}/first_sample"
    DNA_input_dir="${out_dir}/DNA"
fi

################################################################################################################

#reditools needs python2
source /home/cxn/miniconda3/etc/profile.d/conda.sh
conda activate reditools_py2.7


find ${RNA_input_dir} -maxdepth 1  -name "*sense.bam" | \
xargs -n 1 -P ${parallel_num} -I {} \
sh -c '
  RNA_bam="{}"
  bam_name=$(basename "${RNA_bam}")
  echo "[`date "+%F %T"`] Processing RNA: ${bam_name}"
  if compgen -G "$4"/*.bam > /dev/null; then
    for DNA in "$4"/*.bam; do
        DNA_name=$(basename ${DNA})
        python /mnt/san3/usr/cxn/soft/REDItools/NPscripts/REDItoolDnaRnav13.py \
              -i ${RNA_bam} \
              -j ${DNA} \
              -f "$1" \
              -m 20,20 \
              -q 25,25 \
              -c 5,5 \
              -v 2 \
              -a 6-0 \
              -n 0.03 \
              -t "$3" \
              -z -e -u -U \
              -F ${bam_name}_vs_${DNA_name} \
              -o  "$2"/output
    done
  else
        python /mnt/san3/usr/cxn/soft/REDItools/NPscripts/REDItoolDnaRnav13.py \
              -i ${RNA_bam} \
              -f "$1" \
              -m 20,20 \
              -q 25,25 \
              -c 5,5 \
              -v 2 \
              -a 6-0 \
              -n 0.03 \
              -t "$3" \
              -z -e -u \
              -F ${bam_name}_vs_noDNA. \
              -o  "$2"/output   
  fi
  ' _ ${reference} ${out_dir} ${threads} ${DNA_input_dir} ${RNA_input_dir}


########################################################################################################
# filter the editing sites by 1) DNA has no mutant reads; 2) AG editing in positive genes and TC editing in negative genes
conda activate base
awk '{print "##contig=<ID="$1",length="$2">"}' ${reference}.fai > ${out_dir}/vcf_header.txt

find "${out_dir}/output" -type f -path "*.sense.*/outTable*" | while read -r file; do
  direction=sense
  foldername=$(basename "$(dirname "$file")")
  sample=${foldername#DnaRna_}
  sample=${sample%%.sort.dedup*}
  DNA_sample=$(echo ${foldername} | grep -oP '(?<=_vs_)[^.]+' )
  python "/home/cxn/Data/RNAediting/transcriptome/anno_reditool_res.py" -i ${file} -o ${out_dir}/output/ -s ${species} -r ${gene_region} -header ${out_dir}/vcf_header.txt\
        -g ${gff_file} -d 0 -n ${direction}_${sample} -D ${DNA_sample}
done

find "${out_dir}/output" -type f -path "*.antisense.*/outTable*" | while read -r file; do
  direction=antisense
  foldername=$(basename "$(dirname "$file")")
  sample=${foldername#DnaRna_}
  sample=${sample%%.sort.dedup*}
  DNA_sample=$(echo ${foldername} | grep -oP '(?<=_vs_)[^.]+' )
  python "/home/cxn/Data/RNAediting/transcriptome/anno_reditool_res.py" -i ${file} -o ${out_dir}/output/ -s ${species} -r ${gene_region} -header ${out_dir}/vcf_header.txt\
        -g ${gff_file} -d 1 -n ${direction}_${sample} -D ${DNA_sample}
done

find "${out_dir}/output" -type f -path "*.unstrand_sense.*/outTable*" | while read -r file; do
  direction=unstrand_sense
  foldername=$(basename "$(dirname "$file")")
  sample=${foldername#DnaRna_}
  sample=${sample%%.sort.dedup*}
  DNA_sample=$(echo ${foldername} | grep -oP '(?<=_vs_)[^.]+' )
  python "/home/cxn/Data/RNAediting/transcriptome/anno_reditool_res.py" -i ${file} -o ${out_dir}/output/ -s ${species} -r ${gene_region} -header ${out_dir}/vcf_header.txt\
        -g ${gff_file} -d 2 -n ${direction}_${sample} -D ${DNA_sample}
done

python "/home/cxn/Data/RNAediting/transcriptome/merge_reditools_re.py" -i ${out_dir}/output/ -o ${out_dir}/output/

files=(${out_dir}/output/*final.edsite_ann.csv)
head -n 1 "${files[0]}" > ${out_dir}/output/merged_result.txt
for f in "${files[@]}"; do
    tail -n +2 "$f"
done >> ${out_dir}/output/merged_result.txt

