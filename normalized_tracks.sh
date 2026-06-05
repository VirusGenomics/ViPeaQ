#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  normalized_tracks.sh -i INPUT_DIR -o OUTPUT_DIR [-m peaks|bckg]
  normalized_tracks.sh PAIRS_FILE
  normalized_tracks.sh -h

Description:
  Generate normalized bedGraph/bigWig tracks from a ViPeaQ result directory,
  or run the same process for multiple INPUT_DIR / OUTPUT_DIR / MODE entries
  listed in a 3-column TSV or CSV file.

Arguments:
  -i, --input-dir   Input directory
  -o, --output-dir  Output directory
  -m, --mode        Median source: peaks (default) or bckg
  -h, --help        Show this help and exit

Batch mode:
  If a single file is provided instead of -i/-o, it must contain 3 columns:
    column 1 = INPUT_DIR
    column 2 = OUTPUT_DIR
    column 3 = MODE   (peaks or bckg)

Examples:
  normalized_tracks.sh -i ./vipeaq_results -o ./scaled_signal -m peaks
  normalized_tracks.sh samples.tsv
EOF
}

info()  { printf '[INFO] %s\n' "$*"; }
warn()  { printf '[WARN] %s\n' "$*" >&2; }
error() { printf '[ERROR] %s\n' "$*" >&2; }
die()   { error "$*"; exit 1; }

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

trim() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    printf '%s' "$s"
}

detect_delimiter() {
    local file="$1"
    local line
    line=$(awk 'NF{print; exit}' "$file")
    if [[ "$line" == *$'\t'* ]]; then
        printf '\t'
    elif [[ "$line" == *,* ]]; then
        printf ','
    else
        return 1
    fi
}

is_three_column_table() {
    local file="$1"
    local delim
    delim=$(detect_delimiter "$file") || return 1

    awk -F "$delim" '
        NF == 0 { next }
        NF != 3 { exit 1 }
        END { if (NR == 0) exit 1 }
    ' "$file"
}

split_suffix_variants() {
    local s="${1:-}"
    suffix_nounder="${s#_}"
    if [[ -n "$suffix_nounder" ]]; then
        suffix_withunder="_$suffix_nounder"
    else
        suffix_withunder=""
    fi
}

centroid_grid_preserving_bed() {
    local input_bed="$1"
    local output_bed="$2"

    awk '
    BEGIN { OFS="\t" }

    function die(msg) {
        print "ERROR: " msg > "/dev/stderr"
        exit 1
    }

    function isint(x) {
        return (x ~ /^-?[0-9]+$/)
    }

    function flush_chr(    i, step_i, width_i, full_width, step,
                                 first_end, trunc_seen, chr_end,
                                 start_last_regular, last_is_truncated,
                                 new_start, new_end) {
        if (n == 0) return

        if (n == 1) {
            print chrom, s[1], e[1], v[1]
            delete s; delete e; delete v
            n = 0
            return
        }

        step = s[2] - s[1]
        if (step <= 0)
            die("Non-positive start step in chromosome " chrom)

        full_width = e[1] - s[1]
        if (full_width <= 0)
            die("Non-positive width in chromosome " chrom)

        trunc_seen = 0
        chr_end = -1

        # Validate the chromosome block
        for (i = 1; i <= n; i++) {
            width_i = e[i] - s[i]

            if (width_i <= 0)
                die("Invalid interval in chromosome " chrom " at record " i)

            if (i < n) {
                step_i = s[i+1] - s[i]
                if (step_i != step)
                    die("Inconsistent start step in chromosome " chrom " between records " i " and " i+1)
            }

            if (width_i == full_width) {
                if (trunc_seen)
                    die("Full-width interval after truncation started in chromosome " chrom " at record " i)
            } else if (width_i < full_width) {
                if (!trunc_seen) {
                    trunc_seen = 1
                    chr_end = e[i]
                } else if (e[i] != chr_end) {
                    die("Truncated intervals do not share the same chromosome end in chromosome " chrom)
                }
            } else {
                die("Interval wider than expected in chromosome " chrom " at record " i)
            }
        }

        # First centroid boundary:
        # start1 + (full_width + step)/2
        # Example: 0 + (1000+500)/2 = 750
        first_end = s[1] + (full_width + step) / 2

        if (first_end != int(first_end))
            die("Non-integer first centroid boundary in chromosome " chrom)

        first_end = int(first_end)

        last_is_truncated = ((e[n] - s[n]) < full_width)

        # Regular centroid start of the last interval if there were no truncation
        start_last_regular = first_end + (n - 2) * step

        for (i = 1; i <= n; i++) {
            if (i == 1) {
                new_start = s[1]

                if (n == 1) {
                    new_end = e[1]
                } else if (n == 2 && last_is_truncated && start_last_regular > s[n]) {
                    new_end = s[n]
                } else {
                    new_end = first_end
                }

            } else if (i < n) {
                new_start = first_end + (i - 2) * step
                new_end   = first_end + (i - 1) * step

                # If the last interval is truncated and its regular centroid start
                # would overshoot its true start, force the previous interval to end
                # at the true start of the last input interval.
                if (i == n - 1 && last_is_truncated && start_last_regular > s[n]) {
                    new_end = s[n]
                }

            } else {  # i == n
                if (last_is_truncated) {
                    new_start = start_last_regular
                    if (new_start > s[n]) new_start = s[n]
                    new_end = e[n]
                } else {
                    new_start = start_last_regular
                    new_end   = first_end + (n - 1) * step
                }
            }

            if (new_start != int(new_start) || new_end != int(new_end))
                die("Non-integer output boundary in chromosome " chrom)

            new_start = int(new_start)
            new_end   = int(new_end)

            if (new_start >= new_end)
                die("Invalid output interval in chromosome " chrom " at record " i ": " new_start "-" new_end)

            print chrom, new_start, new_end, v[i]
        }

        delete s; delete e; delete v
        n = 0
    }

    {
        if (NF < 4)
            die("Line " NR " has fewer than 4 fields")

        if (!isint($2) || !isint($3))
            die("Non-integer coordinates at line " NR)

        if ($2 >= $3)
            die("Invalid interval at line " NR)

        if (NR == 1) {
            chrom = $1
            prev_chr = $1
            prev_s = $2 + 0
            prev_e = $3 + 0
        } else {
            if ($1 < prev_chr)
                die("Chromosomes are not sorted at line " NR)

            if ($1 == prev_chr) {
                if ($2 < prev_s || ($2 == prev_s && $3 < prev_e))
                    die("Intervals are not sorted within chromosome " $1 " at line " NR)
                if ($2 == prev_s && $3 == prev_e)
                    die("Duplicate interval at line " NR)
            }
        }

        if (NR > 1 && $1 != chrom) {
            flush_chr()
            chrom = $1
        }

        n++
        s[n] = $2 + 0
        e[n] = $3 + 0
        v[n] = $4

        prev_chr = $1
        prev_s = $2 + 0
        prev_e = $3 + 0
    }

    END {
        flush_chr()
    }
    ' "$input_bed" > "$output_bed"
}

process_one() {
    local input_dir="$1"
    local output_dir="$2"
    local mode="${3:-peaks}"
    
    local lambda target filtered_part optional_part unknown_part
    local den_col ratio_bedgraph peaks_file logfile vc basevc
    local medtar medcov factor threads scaled_bedgraph scaled_bigwig
    local mode_file median_file median_value ratio_fpkm_bedgraph ratio_fpkm_scaled_bedgraph
    local medpos="" medneg=""
    
    [[ -d "$input_dir" ]] || {
        warn "Skipping: input directory does not exist: $input_dir"
        return 1
    }

    if [[ -z "$(find "$input_dir" -maxdepth 1 -type f -print -quit)" ]]; then
        warn "Skipping: input directory contains no files: $input_dir"
        return 1
    fi

    mkdir -p "$output_dir" || {
        warn "Skipping: could not create output directory: $output_dir"
        return 1
    }

    info "--------------------------------------------------"
    info "Input directory : $input_dir"
    info "Output directory: $output_dir"

    lambda="0"
    target=""
    filtered_part=""
    optional_part=""
    unknown_part=""
    suffix_nounder=""
    suffix_withunder=""

    info "Searching for target TSV file..."

    local best_score=-999999

    while IFS= read -r path; do
        local file_name cand_filtered cand_optional cand_unknown score
        file_name=${path##*/}

        case "$file_name" in
            positives_*|negatives_*|filtered_positives_*|filtered_negatives_*|*_header.tsv)
                continue
                ;;
        esac

        if [[ $file_name =~ ^(filtered_)?[^[:space:]]+_win_count(_lambda_corrected)?(.*)\.tsv$ ]]; then
            cand_filtered="${BASH_REMATCH[1]}"
            cand_optional="${BASH_REMATCH[2]}"
            cand_unknown="${BASH_REMATCH[3]}"

            score=0

            # Prefer files with an experiment suffix, e.g.
            # _LCL_dB_8wpi_H3K9me3
            [[ -n "$cand_unknown" ]] && ((score += 100))

            # Prefer filtered target files when available.
            [[ -n "$cand_filtered" ]] && ((score += 10))

            # Slightly prefer lambda-corrected files when detected.
            [[ "$cand_optional" == "_lambda_corrected" ]] && ((score += 2))

            # Avoid genome-wide count files when target-specific files exist.
            [[ "$file_name" == genome_win_count* ]] && ((score -= 50))

            if (( score > best_score )); then
                best_score="$score"
                target="$path"
                filtered_part="$cand_filtered"
                optional_part="$cand_optional"
                unknown_part="$cand_unknown"
            fi
        fi
    done < <(find "$input_dir" -maxdepth 1 -type f -name '*.tsv' | sort)

    if [[ -n "$target" ]]; then
        split_suffix_variants "$unknown_part"
        [[ "$optional_part" == "_lambda_corrected" ]] && lambda="1"
    fi

    if [[ -z "$target" ]]; then
        warn "Skipping: no matching target TSV file found in $input_dir"
        return 1
    fi

    info "Target TSV found        : $target"
    info "Experiment suffix       : ${unknown_part:-<none>}"
    info "Suffix without '_'      : ${suffix_nounder:-<none>}"
    info "Suffix with '_'         : ${suffix_withunder:-<none>}"
    info "Lambda mode             : $lambda"
    info "Normalization           : $mode"

    median_file=""
    median_value=""
    medpos=""
    medneg=""

    if [[ "$mode" == "peaks" ]]; then
        info "Searching top positives file..."
        while IFS= read -r path; do
            file_name=${path##*/}
            if [[ $file_name =~ ^top_positives_peaks(.*)\.tsv$ ]]; then
                median_file="$path"
                break
            fi
        done < <(find "$input_dir" -maxdepth 1 -type f -name '*.tsv' | sort)

        [[ -n "$median_file" ]] || {
            warn "Skipping: no file matching ^top_positives_peaks(.*)\\.tsv$ found"
            return 1
        }

        medpos=$(
            awk -F '\t' '
                NF >= 9 && $8 != 0 { print $9 / $8 }
            ' "$median_file" | datamash median 1
        ) || {
            warn "Skipping: failed to compute medpos from $median_file"
            return 1
        }

        median_value="$medpos"
        info "Top positives file     : $median_file"
        info "medpos                 : $medpos"

    else
        info "Searching top negatives file..."
        while IFS= read -r path; do
            file_name=${path##*/}
            if [[ $file_name =~ ^top_negatives_peaks(.*)\.tsv$ ]]; then
                median_file="$path"
                break
            fi
        done < <(find "$input_dir" -maxdepth 1 -type f -name '*.tsv' | sort)

        [[ -n "$median_file" ]] || {
            warn "Skipping: no file matching ^top_negatives_peaks(.*)\\.tsv$ found"
            return 1
        }

        medneg=$(
            awk -F '\t' '
                NF >= 9 && $8 != 0 { print $9 / $8 }
            ' "$median_file" | datamash median 1
        ) || {
            warn "Skipping: failed to compute medneg from $median_file"
            return 1
        }

        median_value="$medneg"
        info "Top negatives file     : $median_file"
        info "medneg                 : $medneg"
    fi
    

    den_col=8
    [[ "$lambda" == "1" ]] && den_col=11

    ratio_fpkm_bedgraph="$output_dir/ratio_FPKM${suffix_withunder}.bedgraph"
    ratio_fpkm_overlap_scaled_bedgraph="$output_dir/ratio_overlap_targetNorm_${mode}${suffix_withunder}.bedgraph"
    ratio_fpkm_scaled_bedgraph="$output_dir/ratio_targetNorm_${mode}${suffix_withunder}.bedgraph"

    info "Writing ratio FPKM bedGraph..."
    awk -F '\t' -v c="$den_col" '
        BEGIN { OFS="\t" }

        function isint(x) {
            return (x ~ /^[0-9]+$/)
        }

        /^#/ { next }
        NF < c || NF < 9 { next }

        {
            # Case 1: BED-like table:
            # chrom  start  end  ...
            if (isint($2) && isint($3)) {
                chrom = $1
                start = $2 + 0
                end   = $3 + 0
            }

            # Case 2: featureCounts/SAF-like table:
            # Geneid  Chr  Start  End  Strand  Length  ...
            # Convert 1-based inclusive Start/End to 0-based half-open BED.
            else if (isint($3) && isint($4)) {
                chrom = $2
                start = ($3 + 0) - 1
                end   = $4 + 0
            }

            # Header or malformed line.
            else {
                next
            }

            if (start < 0) start = 0
            if (start >= end) next

            den = $(c) + 0
            val = (den == 0 ? 0 : ($9 + 0) / den)

            print chrom, start, end, val
        }
    ' "$target" > "$ratio_fpkm_bedgraph"

    if [[ ! -s "$ratio_fpkm_bedgraph" ]]; then
        warn "Skipping: ratio FPKM bedGraph is empty: $ratio_fpkm_bedgraph"
        return 1
    fi
    info "Ratio FPKM bedGraph     : $ratio_fpkm_bedgraph"

    [[ -n "$median_value" ]] || {
        warn "Skipping: median value for mode '$mode' is empty"
        return 1
    }

    if awk "BEGIN{exit !($median_value == 0)}"; then
        warn "Skipping: median value for mode '$mode' is 0"
        return 1
    fi

    info "Writing median-scaled ratio FPKM bedGraph..."
    awk -F '\t' -v m="$median_value" '
        BEGIN { OFS="\t" }
        NF >= 4 { print $1, $2, $3, $4 / m }
    ' "$ratio_fpkm_bedgraph" > "$ratio_fpkm_overlap_scaled_bedgraph"

    if [[ ! -s "$ratio_fpkm_overlap_scaled_bedgraph" ]]; then
        warn "Skipping: scaled ratio FPKM bedGraph is empty: $ratio_fpkm_overlap_scaled_bedgraph"
        return 1
    fi
    info "Scaled ratio FPKM file  : $ratio_fpkm_overlap_scaled_bedgraph"

    if ! centroid_grid_preserving_bed \
        "$ratio_fpkm_overlap_scaled_bedgraph" \
        "$ratio_fpkm_scaled_bedgraph"
    then
        warn "Skipping: failed to convert overlapped bedGraph to centroid-based bedGraph: $ratio_fpkm_overlap_scaled_bedgraph"
        return 1
    fi

    if [[ ! -s "$ratio_fpkm_scaled_bedgraph" ]]; then
        warn "Skipping: centroid-based bedGraph is empty: $ratio_fpkm_scaled_bedgraph"
        return 1
    fi
    info "Centroid ratio FPKM file: $ratio_fpkm_scaled_bedgraph"


    peaks_file=""
    other_mode="peaks"
    [[ "$mode" == "peaks" ]] && other_mode="bckg"

    for candidate in \
        "$input_dir/outfile_peaks_${mode}_${suffix_nounder}.tsv" \
        "$input_dir/outfile_peaks_${mode}_${suffix_withunder}.tsv"
    do
        if [[ -f "$candidate" ]]; then
            peaks_file="$candidate"
            break
        fi
    done

    if [[ -z "$peaks_file" ]]; then
        other_mode_file=""
        for candidate in \
            "$input_dir/outfile_peaks_${other_mode}_${suffix_nounder}.tsv" \
            "$input_dir/outfile_peaks_${other_mode}_${suffix_withunder}.tsv"
        do
            if [[ -f "$candidate" ]]; then
                other_mode_file="$candidate"
                break
            fi
        done

        if [[ -n "$other_mode_file" ]]; then
            warn "Skipping: no peaks summary file was found for requested mode '$mode'."
            warn "A corresponding file exists for mode '$other_mode' instead:"
            warn "  $other_mode_file"
            warn "This suggests that only the '$other_mode' normalization has been generated for this dataset so far."
            warn "Please run 'Plots_ViPeaQ.R -m $mode' first, then rerun normalized_tracks.sh."
            return 1
        fi

        warn "Skipping: peaks file not found. Tried:"
        warn "  $input_dir/outfile_peaks_${mode}_${suffix_nounder}.tsv"
        warn "  $input_dir/outfile_peaks_${mode}_${suffix_withunder}.tsv"
        return 1
    fi

    info "Computing median normalized coverage..."
    medtar=$(
        tail -n +2 "$peaks_file" | cut -f3 | datamash median 1
    ) || {
        warn "Skipping: failed to compute normalized median from $peaks_file"
        return 1
    }
    info "Median normalized cov   : $medtar"

        info "Searching matching log file..."

    logfile=$(
        {
            printf "%s\n" "$input_dir"/????-??-??_??-??-??_"${suffix_nounder}".log
            printf "%s\n" "$input_dir"/????-??-??_??-??-??"${suffix_withunder}".log
        } 2>/dev/null \
        | awk '!seen[$0]++' \
        | while IFS= read -r f; do [[ -f "$f" ]] && printf '%s\n' "$f"; done \
        | sort \
        | tail -n 1
    )

    if [[ -z "$logfile" || ! -f "$logfile" ]]; then
        warn "Skipping: no matching log file found for suffix variants:"
        warn "  $input_dir/????-??-??_??-??-??_${suffix_nounder}.log"
        warn "  $input_dir/????-??-??_??-??-??${suffix_withunder}.log"
        return 1
    fi
    info "Log file found          : $logfile"

    vc=$(awk -F' = ' '/^--vc[[:space:]]*=/{print $2; exit}' "$logfile")
    if [[ -z "$vc" ]]; then
        warn "Skipping: could not extract --vc from log file: $logfile"
        return 1
    fi
    if [[ ! -f "$vc" ]]; then
        warn "Skipping: extracted --vc BAM does not exist: $vc"
        return 1
    fi
    info "Raw target BAM          : $vc"

    basevc=$(basename "$vc" .bam)

    info "Computing median raw coverage..."
    medcov=$(
        samtools depth "$vc" | cut -f3 | datamash median 1
    ) || {
        warn "Skipping: failed to compute raw median from $vc"
        return 1
    }
    info "Median raw cov          : $medcov"

    if awk "BEGIN{exit !($medcov == 0)}"; then
        warn "Skipping: median raw coverage is 0"
        return 1
    fi

    factor=$(echo "scale=12; $medtar / $medcov" | bc -l)
    [[ -n "$factor" ]] || {
        warn "Skipping: failed to compute scale factor"
        return 1
    }
    info "Scale factor            : $factor"

    threads=10
    scaled_bedgraph="$output_dir/${basevc}.bedgraph"
    scaled_bigwig="$output_dir/${basevc}.bw"

    info "Running bamCoverage (bedGraph)..."
    bamCoverage \
        -b "$vc" \
        -of bedgraph \
        -bs 1 \
        -p "$threads" \
        --normalizeUsing None \
        --scaleFactor "$factor" \
        -o "$scaled_bedgraph"

    info "Running bamCoverage (bigWig)..."
    bamCoverage \
        -b "$vc" \
        -of bigwig \
        -bs 1 \
        -p "$threads" \
        --normalizeUsing None \
        --scaleFactor "$factor" \
        -o "$scaled_bigwig"

    info "Done"
    info "Generated:"
    info "  - $ratio_fpkm_scaled_bedgraph"
    info "  - $ratio_fpkm_overlap_scaled_bedgraph"
    info "  - $ratio_fpkm_bedgraph"
    info "  - $scaled_bedgraph"
    info "  - $scaled_bigwig"

    return 0
}

process_batch_file() {
    local pairs_file="$1"
    local delim
    local line_no=0
    local ok=0
    local fail=0
    local failed_entries=()

    [[ -f "$pairs_file" ]] || die "Batch file does not exist: $pairs_file"
    is_three_column_table "$pairs_file" || die "Batch file must be a non-empty 3-column TSV or CSV: $pairs_file"

    delim=$(detect_delimiter "$pairs_file") || die "Could not detect delimiter in batch file: $pairs_file"

    info "Batch mode detected"
    info "Pairs file: $pairs_file"

    while IFS= read -r line || [[ -n "$line" ]]; do
        ((line_no += 1))
        [[ -n "$(trim "$line")" ]] || continue

        local in_dir out_dir mode
        in_dir=$(printf '%s\n' "$line" | awk -F "$delim" '{print $1}')
        out_dir=$(printf '%s\n' "$line" | awk -F "$delim" '{print $2}')
        mode=$(printf '%s\n' "$line" | awk -F "$delim" '{print $3}')

        in_dir=$(trim "$in_dir")
        out_dir=$(trim "$out_dir")
        mode=$(trim "$mode")

        if [[ -z "$in_dir" || -z "$out_dir" || -z "$mode" ]]; then
            warn "Line $line_no: empty INPUT_DIR, OUTPUT_DIR, or MODE, skipping"
            failed_entries+=("line $line_no | INPUT_DIR='$in_dir' | OUTPUT_DIR='$out_dir' | MODE='$mode' | reason=empty INPUT_DIR, OUTPUT_DIR, or MODE")
            ((fail += 1))
            continue
        fi

        if [[ "$mode" != "peaks" && "$mode" != "bckg" ]]; then
            warn "Line $line_no: invalid MODE '$mode', skipping"
            failed_entries+=("line $line_no | INPUT_DIR='$in_dir' | OUTPUT_DIR='$out_dir' | MODE='$mode' | reason=invalid MODE")
            ((fail += 1))
            continue
        fi

        if process_one "$in_dir" "$out_dir" "$mode"; then
            ((ok += 1))
        else
            ((fail += 1))
            failed_entries+=("line $line_no | INPUT_DIR='$in_dir' | OUTPUT_DIR='$out_dir' | MODE='$mode'")
        fi
    done < "$pairs_file"

    info "=================================================="
    info "Batch completed"
    info "Successful runs : $ok"
    info "Failed runs     : $fail"

    if [[ ${#failed_entries[@]} -gt 0 ]]; then
        warn "Failed entries:"
        for entry in "${failed_entries[@]}"; do
            warn "  - $entry"
        done
    fi

    [[ "$ok" -gt 0 ]] || return 1
}

main() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 1
    fi

    if [[ $# -eq 1 ]]; then
        case "$1" in
            -h|--help)
                usage
                exit 0
                ;;
        esac

        require_cmd awk
        require_cmd sort
        require_cmd samtools
        require_cmd bamCoverage
        require_cmd bc
        require_cmd tail
        require_cmd cut
        require_cmd find
        require_cmd datamash

        if [[ -f "$1" ]]; then
            process_batch_file "$1"
            exit $?
        fi

        die "Single-argument mode expects a 3-column batch TSV/CSV file or -h."
    fi

    local input_dir=""
    local output_dir=""
    local mode="peaks"

    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--input-dir)
                [[ $# -ge 2 ]] || die "Missing value for $1"
                input_dir="$2"
                shift 2
                ;;
            -o|--output-dir)
                [[ $# -ge 2 ]] || die "Missing value for $1"
                output_dir="$2"
                shift 2
                ;;
            -m|--mode)
                [[ $# -ge 2 ]] || die "Missing value for $1"
                mode="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                die "Unknown argument: $1. Use -h for help."
                ;;
        esac
    done

    require_cmd awk
    require_cmd sort
    require_cmd samtools
    require_cmd bamCoverage
    require_cmd bc
    require_cmd tail
    require_cmd cut
    require_cmd find
    require_cmd datamash

    [[ -n "$input_dir" ]] || die "Input directory is required. Use -i."
    [[ -n "$output_dir" ]] || die "Output directory is required. Use -o."
    [[ "$mode" == "peaks" || "$mode" == "bckg" ]] || die "Invalid mode: $mode. Allowed values: peaks, bckg"

    process_one "$input_dir" "$output_dir" "$mode"
}

main "$@"