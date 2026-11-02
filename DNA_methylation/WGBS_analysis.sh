#!/bin/bash
set -uo pipefail

# Load configuration from the first argument (default: config.sh).
config_file="${1:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/config.sh}"
if [[ ! -f "${config_file}" ]]; then
    echo "ERROR: Configuration not found: ${config_file}" >&2
    echo "Copy config.example.sh to config.sh and edit the paths." >&2
    exit 1
fi
# shellcheck source=config.example.sh
source "${config_file}"

# Required configuration variables
for name in methy_dir genome_dir fq_dir bismark_dir region_bed tmp_dir \
            lambda_chr lambda_min_cov min_cov max_jobs align_threads extractor_threads; do
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: Missing configuration variable: ${name}" >&2
        exit 1
    fi
done

mkdir -p "${methy_dir}" "${tmp_dir}" || exit 1

# Failed sample list
failed_list="${methy_dir}/failed_samples.txt"
: > "${failed_list}" || exit 1

if [ ! -d "${fq_dir}" ]; then
    echo "ERROR: FASTQ directory not found: ${fq_dir}"
    exit 1
fi

process_sample() {
    local r1="$1"
    local sample
    local r2
    local outdir
    local log
    local bam
    local dedup_bam
    local cx_report
    local lambda_report
    local class_report

    sample=$(basename "${r1}" _1.fq.gz)
    r2="${r1%_1.fq.gz}_2.fq.gz"
    outdir="${methy_dir}/${sample}"
    log="${outdir}/${sample}.log"

    if ! mkdir -p "${outdir}" "${tmp_dir}/${sample}"; then
        echo "ERROR: Cannot create directories for: ${sample}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    {
        echo ""
        echo "========== New run: $(date) =========="
        echo "Sample: ${sample}"
        echo "R1: ${r1}"
        echo "R2: ${r2}"
        echo "TMP: ${tmp_dir}/${sample}"
        echo "Start time: $(date)"
        echo ""
    } >> "${log}"

    echo "Processing sample: ${sample}"

    if [ ! -f "${r2}" ]; then
        echo "ERROR: R2 file not found: ${r2}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    # 1. Alignment
    bam="${outdir}/${sample}_pe.bam"
    
    if [ -s "${bam}" ] && samtools quickcheck "${bam}" 2>/dev/null; then
        echo "[1/5] Alignment already completed, skipping: ${sample}" | tee -a "${log}"
    else
        echo "[1/5] Bismark alignment: ${sample}" | tee -a "${log}"

        if ! "${bismark_dir}/bismark" \
            --genome "${genome_dir}" \
            -o "${outdir}" \
            -basename "${sample}" \
            -bowtie2 \
            -p "${align_threads}" \
            --temp_dir "${tmp_dir}/${sample}" \
            -1 "${r1}" \
            -2 "${r2}" \
            &>> "${log}"
        then
            echo "ERROR: Bismark alignment failed: ${sample}" | tee -a "${log}"
            echo "${sample}" >> "${failed_list}"
            return 1
        fi
    fi

    if [ ! -s "${bam}" ]; then
        echo "ERROR: Bismark BAM missing or empty: ${bam}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi


    # 2. Deduplication
    dedup_bam="${outdir}/${sample}_pe.deduplicated.bam"
    
    if [ -s "${dedup_bam}" ] && samtools quickcheck "${dedup_bam}" 2>/dev/null; then
        echo "[2/5] Deduplication already completed, skipping: ${sample}" | tee -a "${log}"
    else
        echo "[2/5] Deduplication: ${sample}" | tee -a "${log}"

        if ! "${bismark_dir}/deduplicate_bismark" \
            -p \
            --bam \
            --output_dir "${outdir}" \
            "${bam}" \
            &>> "${log}"
        then
            echo "ERROR: Deduplication failed: ${sample}" | tee -a "${log}"
            echo "${sample}" >> "${failed_list}"
            return 1
        fi
    fi

    if [ ! -s "${dedup_bam}" ]; then
        echo "ERROR: Deduplicated BAM missing or empty: ${dedup_bam}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi


    # 3. Methylation extraction
    cx_report="${outdir}/${sample}_pe.deduplicated.CX_report.txt.gz"
    
    if [ -s "${cx_report}" ] && gzip -t "${cx_report}" 2>/dev/null; then
        echo "[3/5] Methylation extraction already completed, skipping: ${sample}" | tee -a "${log}"
    else
        echo "[3/5] Methylation extraction: ${sample}" | tee -a "${log}"

        if ! "${bismark_dir}/bismark_methylation_extractor" \
            -p \
            --gzip \
            --buffer_size 20G \
            --comprehensive \
            --bedGraph \
            --CX \
            --cytosine_report \
            --multicore "${extractor_threads}" \
            --genome_folder "${genome_dir}" \
            -o "${outdir}" \
            "${dedup_bam}" \
            &>> "${log}"
        then
            echo "ERROR: Methylation extraction failed: ${sample}" | tee -a "${log}"
            echo "${sample}" >> "${failed_list}"
            return 1
        fi
    fi

    # 4. Lambda DNA non-conversion rate
    echo "[4/5] Lambda DNA non-conversion rate: ${sample}" | tee -a "${log}"

    cx_report="${outdir}/${sample}_pe.deduplicated.CX_report.txt.gz"
    lambda_report="${outdir}/${sample}.lambda_nonconversion.tsv"

    if [ ! -s "${cx_report}" ]; then
        echo "ERROR: CX report missing or empty: ${cx_report}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    if ! (
        gzip -cd "${cx_report}" |
        awk -v sample="${sample}" -v chr="${lambda_chr}" -v min_cov="${lambda_min_cov}" '
            BEGIN {
                OFS = "\t"
                print "Sample", "Chromosome", "Covered_sites", \
                      "Methylated_C", "Unmethylated_C", \
                      "Non_conversion_percent", "Conversion_percent"
            }

            $1 == chr && ($4 + $5) >= min_cov {
                m += $4
                u += $5
                covered++
            }

            END {
                if (m + u > 0) {
                    printf "%s\t%s\t%d\t%.0f\t%.0f\t%.6f\t%.6f\n", \
                           sample, chr, covered, m, u, \
                           100*m/(m+u), 100*u/(m+u)
                } else {
                    print "ERROR: No covered lambda cytosines found" > "/dev/stderr"
                    exit 1
                }
            }
        '
    ) > "${lambda_report}.tmp" 2>> "${log}"
    then
        rm -f "${lambda_report}.tmp"
        echo "ERROR: Lambda non-conversion calculation failed: ${sample}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    if ! mv "${lambda_report}.tmp" "${lambda_report}"; then
        echo "ERROR: Could not save lambda report: ${sample}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    cat "${lambda_report}" | tee -a "${log}"


    # 5. Methylation level by region class
    echo "[5/5] Methylation by region class: ${sample}" | tee -a "${log}"

    class_report="${outdir}/${sample}.methylation_by_class.tsv"

    if ! (
        gzip -cd "${cx_report}" |
        awk -v min_cov="${min_cov}" -v lambda="${lambda_chr}" '
            BEGIN { FS = OFS = "\t" }

            # CX report is 1-based; convert to 1-bp BED intervals.
            # Output: chr, start, end, mC, uC, context
            $1 != lambda && ($4 + $5) >= min_cov {
                print $1, $2 - 1, $2, $4, $5, $6
            }
        ' |
        bedtools intersect \
            -a stdin \
            -b "${region_bed}" \
            -wa -wb |
        awk -v sample="${sample}" '
            BEGIN {
                FS = OFS = "\t"
                n = split("ALL CG CHG CHH", contexts, " ")
            }
            
            # Read BED and calculate region lengths
            NR == FNR {
                category = $4
                len = $3 - $2

                categories[category] = 1
                categories["TOTAL"] = 1

                region_length[category] += len
                region_length["TOTAL"] += len
                next
            }
            
            # Intersect output:
            # 1–6: cytosine BED
            # 7–10: region BED; column 10 is category
            {
                category = $10
                context = $6

                # Accumulate both class-specific and TOTAL statistics
                for (j = 1; j <= 2; j++) {
                    cat = (j == 1 ? category : "TOTAL")

                    key = cat SUBSEP "ALL"
                    m[key] += $4
                    u[key] += $5
                    sites[key]++

                    key = cat SUBSEP context
                    m[key] += $4
                    u[key] += $5
                    sites[key]++
                }
            }

            END {
                print "Sample", "Category", "Context", \
                      "Region_length_bp", "Covered_sites", \
                      "Methylated_C", "Unmethylated_C", \
                      "Total_C_calls", "Methylation_percent"

                for (category in categories) {
                    for (i = 1; i <= n; i++) {
                        context = contexts[i]
                        key = category SUBSEP context
                        total = m[key] + u[key]

                        percent = total > 0 \
                            ? sprintf("%.6f", 100 * m[key] / total) \
                            : "NA"

                        printf "%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%s\n", \
                               sample, category, context, \
                               region_length[category], \
                               sites[key] + 0, \
                               m[key] + 0, u[key] + 0, \
                               total, percent
                    }
                }
            }
        ' "${region_bed}" -
    ) > "${class_report}.tmp" 2>> "${log}"
    
    then
        rm -f "${class_report}.tmp"
        echo "ERROR: Methylation by class failed: ${sample}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    if ! mv -- "${class_report}.tmp" "${class_report}"; then
        echo "ERROR: Saving class report failed: ${sample}" | tee -a "${log}"
        echo "${sample}" >> "${failed_list}"
        return 1
    fi

    tee -a "${log}" < "${class_report}"


    {
        echo ""
        echo "Finished sample: ${sample}"
        echo "End time: $(date)"
    } >> "${log}"

    rm -rf -- "${tmp_dir}/${sample}"

    echo "Finished sample: ${sample}"
}



# Run samples in parallel
sample_count=0

while IFS= read -r -d '' r1; do
    sample_count=$((sample_count + 1))

    process_sample "${r1}" &

    while [ "$(jobs -rp | wc -l)" -ge "${max_jobs}" ]; do
        wait -n 2>/dev/null || true
    done
done < <(find "${fq_dir}" -name "*_1.fq.gz" -print0 | sort -z)

wait

if [ "${sample_count}" -eq 0 ]; then
    echo "ERROR: No *_1.fq.gz files found in: ${fq_dir}"
    exit 1
fi

# Final report
if [ -s "${failed_list}" ]; then
    echo "Some samples failed. See: ${failed_list}"
    cat "${failed_list}"
    exit 1
else
    rm -f "${failed_list}"
    echo "All samples finished successfully."
fi
