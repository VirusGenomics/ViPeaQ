#!/usr/bin/env bash
set -u
set -o pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <parameter_file> <output_base_path>" >&2
    exit 1
fi

PARAM_FILE=$1
BASE_PATH=$2

# Use ViPeaQ.sh located in the same directory as this wrapper
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPT="${SCRIPT_DIR}/ViPeaQ.sh"

LOG_FILE="execution_times.log"

if [[ ! -f "$PARAM_FILE" ]]; then
    echo "Error: parameter file not found: $PARAM_FILE" >&2
    exit 1
fi

if [[ ! -x "$SCRIPT" ]]; then
    echo "Error: target script is not executable or not found: $SCRIPT" >&2
    exit 1
fi

mkdir -p "$BASE_PATH" || {
    echo "Error: could not create base output directory: $BASE_PATH" >&2
    exit 1
}

BASE_PATH="$(cd "$BASE_PATH" && pwd -P)"

copy_vc_tracks() {
    local vc_path=$1
    local run_outdir=$2
    local dest_dir=$3

    local vc_base
    local vc_stem
    local copied=0
    local f

    vc_base="$(basename -- "$vc_path")"

    # Expected case: something.bam -> something
    vc_stem="${vc_base%.bam}"

    # Fallback if the vc file does not end with .bam
    if [[ "$vc_stem" == "$vc_base" ]]; then
        vc_stem="${vc_base%.*}"
    fi

    shopt -s nullglob
    local track_files=(
        "$run_outdir/${vc_stem}"*.bedgraph
        "$run_outdir/${vc_stem}"*.bw
    )
    shopt -u nullglob

    if [[ ${#track_files[@]} -eq 0 ]]; then
        echo "  warning: no track files found for vc basename: $vc_stem"
        echo "  searched in: $run_outdir"
        echo "  expected pattern:"
        echo "    ${run_outdir}/${vc_stem}*.bedgraph"
        echo "    ${run_outdir}/${vc_stem}*.bw"
        return 1
    fi

    for f in "${track_files[@]}"; do
        if [[ -f "$f" ]]; then
            cp -p -- "$f" "$dest_dir"/
            echo "  copied: $f -> $dest_dir/"
            ((copied++))
        fi
    done

    echo "  copied_track_files: $copied"
}

# Read first line, remove CR if present, split on spaces/tabs
HEADER=$(head -n 1 "$PARAM_FILE" | tr -d '\r')
IFS=$' \t' read -r -a PARAM_NAMES <<< "$HEADER"

if [[ ${#PARAM_NAMES[@]} -eq 0 ]]; then
    echo "Error: empty header in $PARAM_FILE" >&2
    exit 1
fi

# Remove accidental leading dashes from header names
# Handles vc, -vc, and --vc
for i in "${!PARAM_NAMES[@]}"; do
    PARAM_NAMES[$i]=${PARAM_NAMES[$i]#--}
    PARAM_NAMES[$i]=${PARAM_NAMES[$i]#-}
done

{
    echo "Execution times (in seconds):"
    echo "Parameter file: $PARAM_FILE"
    echo "Base output path: $BASE_PATH"
    echo "Target script: $SCRIPT"
    echo
} > "$LOG_FILE"

line_no=1

while IFS= read -r LINE || [[ -n "$LINE" ]]; do
    ((line_no++))

    LINE=${LINE//$'\r'/}

    [[ -z "${LINE//[[:space:]]/}" ]] && continue
    [[ "$LINE" =~ ^[[:space:]]*# ]] && continue

    IFS=$' \t' read -r -a PARAM_VALUES <<< "$LINE"

    if [[ ${#PARAM_VALUES[@]} -ne ${#PARAM_NAMES[@]} ]]; then
        {
            echo "Line $line_no"
            echo "  error: expected ${#PARAM_NAMES[@]} values, got ${#PARAM_VALUES[@]}"
            echo
        } >> "$LOG_FILE"
        continue
    fi

    cmd=( "$SCRIPT" )

    vc_value=""
    out_value=""

    for i in "${!PARAM_NAMES[@]}"; do
        name=${PARAM_NAMES[$i]}
        value=${PARAM_VALUES[$i]}

        if [[ "$name" == "o" ]]; then
            mkdir -p "$value" || {
                {
                    echo "Line $line_no"
                    echo "  error: failed to create output directory: $value"
                    echo
                } >> "$LOG_FILE"
                continue 2
            }

            # Store absolute output path for later copy step
            value="$(cd "$value" && pwd -P)"
            out_value="$value"
        fi

        if [[ "$name" == "vc" ]]; then
            vc_value="$value"
        fi

        # ViPeaQ.sh uses getopt long options; --vc/--o are safe here.
        cmd+=( "--$name" "$value" )
    done

    start_time=${EPOCHREALTIME:-$(date +%s.%N)}
    "${cmd[@]}"
    exit_code=$?
    end_time=${EPOCHREALTIME:-$(date +%s.%N)}

    duration=$(awk -v s="$start_time" -v e="$end_time" 'BEGIN { printf "%.6f", (e - s) }')

    {
        echo "Line $line_no"
        printf "  command: "
        printf "%q " "${cmd[@]}"
        printf "\n"
        echo "  exit_code: $exit_code"
        echo "  duration: $duration seconds"
    } >> "$LOG_FILE"

    if [[ "$exit_code" -eq 0 ]]; then
        {
            echo "  track copy:"

            if [[ -z "$vc_value" ]]; then
                echo "  warning: no vc parameter found on this line; cannot copy vc-based tracks"
            elif [[ -z "$out_value" ]]; then
                echo "  warning: no o parameter found on this line; cannot find run output directory"
            else
                copy_vc_tracks "$vc_value" "$out_value" "$BASE_PATH" || true
            fi

            echo
        } >> "$LOG_FILE" 2>&1
    else
        {
            echo "  track copy: skipped because ViPeaQ.sh failed"
            echo
        } >> "$LOG_FILE"
    fi

done < <(tail -n +2 "$PARAM_FILE")

echo "Execution completed. Times logged in $LOG_FILE"
echo "Collected track files copied to: $BASE_PATH"