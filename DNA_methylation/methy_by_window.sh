#!/usr/bin/env bash
set -euo pipefail

# Calculate per-sample methylation levels in supplied BED windows.
# Input: Bismark *_pe.deduplicated.CX_report.txt.gz (1-based C positions).
# Only C sites with methylated + unmethylated reads >= --min-cov are retained.
# Window methylation (%) = 100 * sum(methylated reads) / sum(total reads).
# Samples are NEVER pooled.

usage() {
  cat <<'USAGE'
Usage:
  bash cx_window_methylation.sh \
    --sample_dir DIR --bed FILE --outdir DIR \
    [--tmpdir DIR] [--min-cov 5] [--window-size 200] 
    
Required:
  --sample_dir     Root directory of sample folders (recursive search)
  --bed            BED windows (0-based, half-open; non-overlapping)
  --outdir         Output directory

Optional:
  --tmpdir         Temporary parent directory (default: ${TMPDIR:-/tmp})
  --min-cov        Minimum coverage per C site (default: 5)
  --window-size    Expected maximum BED interval length (default: 200);
                   shorter terminal windows are allowed; does not make windows
  --help           Show help

Per sample: <sample>.window<WINDOW_SIZE>.cov<MIN_COV>.bed
Columns: chrom, start, end, n_C_sites, methylated_reads,
         total_reads, methylation_percent
Other output: sample_manifest.tsv
USAGE
}

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

SAMPLE_DIR=""
BED=""
OUTDIR=""
TMP_PARENT="${TMPDIR:-/tmp}"
MIN_COV=5
WINDOW_SIZE=200

while (($#)); do
  case "$1" in
    --sample_dir|--bed|--outdir|--tmpdir|--min-cov|--window-size)
      (($# >= 2)) || die "Missing value for $1"
      [[ -n "$2" && "$2" != --* ]] || die "Missing value for $1"
      case "$1" in
        --sample_dir) SAMPLE_DIR=$2 ;;
        --bed) BED=$2 ;;
        --outdir) OUTDIR=$2 ;;
        --tmpdir) TMP_PARENT=$2 ;;
        --min-cov) MIN_COV=$2 ;;
        --window-size) WINDOW_SIZE=$2 ;;
      esac
      shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$SAMPLE_DIR" && -n "$BED" && -n "$OUTDIR" ]] || { usage >&2; die 'Missing required argument'; }
[[ -d "$SAMPLE_DIR" ]] || die "Sample directory not found: $SAMPLE_DIR"
[[ -f "$BED" ]] || die "BED not found: $BED"
[[ "$MIN_COV" =~ ^[1-9][0-9]*$ ]] || die '--min-cov must be a positive integer'
[[ "$WINDOW_SIZE" =~ ^[1-9][0-9]*$ ]] || die '--window-size must be a positive integer'
for cmd in bedtools awk sort gzip find mktemp; do
  command -v "$cmd" >/dev/null 2>&1 || die "Missing command: $cmd"
done
mkdir -p -- "$OUTDIR" "$TMP_PARENT"
[[ -w "$OUTDIR" && -w "$TMP_PARENT" ]] || die 'Output or temporary directory is not writable'
SAMPLE_DIR=$(cd "$SAMPLE_DIR" && pwd -P)
OUTDIR=$(cd "$OUTDIR" && pwd -P)
TMP_PARENT=$(cd "$TMP_PARENT" && pwd -P)
export LC_ALL=C
WORKDIR=$(mktemp -d "$TMP_PARENT/cx_window.XXXXXXXX")
trap 'rm -rf -- "$WORKDIR"' EXIT

printf '[INFO] Preparing windows: %s\n' "$BED"
awk -v max_size="$WINDOW_SIZE" 'BEGIN {FS="\t"; OFS="\t"}
  /^#/ || /^track([ \t]|$)/ || /^browser([ \t]|$)/ || /^[[:space:]]*$/ {next}
  NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 {
    print "Invalid BED interval on line " NR > "/dev/stderr"; exit 1
  }
  ($3-$2)>max_size {
    print "BED interval exceeds --window-size on line " NR > "/dev/stderr"; exit 1
  }
  {print $1,$2,$3}
' "$BED" | sort -T "$WORKDIR" \
    -k1,1 -k2,2n -k3,3n > "$WORKDIR/windows.bed"
[[ -s "$WORKDIR/windows.bed" ]] || die 'BED has no intervals'
awk 'BEGIN{FS="\t"} NR>1 && $1==chr && $2<prev_end {
  print "Overlapping BED intervals: " $1 ":" $2 "-" $3 > "/dev/stderr"; exit 1
} {chr=$1; prev_end=$3}' "$WORKDIR/windows.bed"

MANIFEST="$OUTDIR/sample_manifest.tsv"
printf 'sample\tinput_file\toutput_file\n' > "$MANIFEST"
count=0
declare -A SEEN_OUTPUT=()

# Use process substitution to preserve the while loop counters.
while IFS= read -r -d '' CX; do
  count=$((count+1))
  REL=${CX#"$SAMPLE_DIR"/}
  PARENT=${REL%/*}
  [[ "$PARENT" != "$REL" ]] || PARENT=root
  SAMPLE=${PARENT//\//__}
  SAMPLE=${SAMPLE//[^a-zA-Z0-9._-]/_}
  OUTPUT="$OUTDIR/${SAMPLE}.window${WINDOW_SIZE}.cov${MIN_COV}.bed"
  [[ -z "${SEEN_OUTPUT[$OUTPUT]:-}" ]] || die "Multiple inputs map to the same output: $OUTPUT"
  SEEN_OUTPUT[$OUTPUT]=1
  printf '[INFO] Sample: %s\n' "$SAMPLE"

  # CX_report fields: chr, position(1-based), strand, methylated, unmethylated, context, trinucleotide.
  # Convert position to 0-based BED interval [pos-1,pos).
  gzip -dc -- "$CX" |
    awk -v min_cov="$MIN_COV" 'BEGIN{FS=OFS="\t"}
      NF>=5 && $2~/^[1-9][0-9]*$/ && $4~/^[0-9]+$/ && $5~/^[0-9]+$/ {
        m=$4+0; total=m+$5
        if(total>=min_cov) print $1,$2-1,$2,m,total,1
      }' |
    sort -T "$WORKDIR" \
      -k1,1 -k2,2n -k3,3n > "$WORKDIR/sites.bed"

  # bedtools map -c 4,5,6 -o sum,sum,sum gives per-window read sums and site count.
  # For a completely empty sites BED, write zero/NA windows without calling bedtools.
  {
    printf '#chrom\tstart\tend\tn_C_sites\tmethylated_reads\ttotal_reads\tmethylation_percent\n'
    if [[ -s "$WORKDIR/sites.bed" ]]; then
      bedtools map -sorted -a "$WORKDIR/windows.bed" -b "$WORKDIR/sites.bed" \
        -c 4,5,6 -o sum,sum,sum -null 0 |
      awk 'BEGIN {FS=OFS="\t"}
        {
          m=$4+0; total=$5+0; n=$6+0
          if(total>0) printf "%s\t%s\t%s\t%d\t%.0f\t%.0f\t%.6f\n", $1,$2,$3,n,m,total,100*m/total
          else print $1,$2,$3,0,0,0,"NA"
        }'
    else
      awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,0,0,0,"NA"}' "$WORKDIR/windows.bed"
    fi
  } > "$OUTPUT"

  printf '%s\t%s\t%s\n' "$SAMPLE" "$CX" "$OUTPUT" >> "$MANIFEST"
  printf '[INFO] Finished: %s\n' "$SAMPLE"
done < <(find "$SAMPLE_DIR" -type f -name '*_pe.deduplicated.CX_report.txt.gz' -print0 | sort -z)

((count>0)) || die "No matching CX_report files found under $SAMPLE_DIR"
printf '[INFO] Finished %d samples. Output directory: %s\n' "$count" "$OUTDIR"
