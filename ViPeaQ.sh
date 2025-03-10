#!/bin/bash

# 2025
# A. Robitaille (alexis.robitaille@leibniz-liv.de)
# Leibniz Institute of Virology - Hamburg

# exit when any command fails
set -e
set +o pipefail  # Ensure pipeline failures are caught

##################
##	Functions	##
##################
to_absolute_path() {
    local input_path="$1"

    if [[ -z "$input_path" ]]; then
        echo "Error: Path is empty." >&2
        return 1
    fi

    if command -v realpath >/dev/null 2>&1; then
        realpath "$input_path"
        return
    fi

    if [[ -d "$input_path" || -f "$input_path" ]]; then
        echo "$(cd "$(dirname "$input_path")" && pwd)/$(basename "$input_path")"
    else
        echo "Error: Path does not exist or is invalid: $input_path" >&2
        return 1
    fi
}

download_file() {
    local url=$1
    local dest=$2
	echo "Attempting to download from: $url";
    # echo "Saving to: $dest";
    if [ ! -f "$dest" ]; then
        wget -N -q "$url" -O "$dest" || { echo "Failed to download $url"; return 1; }
    else
        echo "$dest already exists. Skipping download.";
    fi
}

download_genome_files() {
	local genome=$1
	case "$genome" in
		hg19)
			download_file "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes" "$out_genome/hg19.chrom.sizes"
			download_file "https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz" "$out_genome/ENCFF001TDO.bed.gz"
			gunzip -c ENCFF001TDO.bed.gz |cut -f 1,2,3 > ENCFF001TDO.bed
			#download_file "https://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz" "$out_genome/hg19_repeatmasker.bed.gz"
			#gunzip -c hg19.fa.out.gz | sed 's/ \+/\t/g' | sed 's/^\t//g' | tail -n +4 | cut -f 5,6,7 > hg19_repeatmasker
			download_file "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/gap.txt.gz" "$out_genome/gap.txt.gz"
			gunzip -c gap.txt.gz | cut -f 2,3,4 > hg19_N.bed
			;;
		mm10)
			download_file "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes" "$out_genome/mm10.chrom.sizes"
			download_file "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz" "$out_genome/ENCFF547MET.bed.gz"
			gunzip -c ENCFF547MET.bed.gz |cut -f 1,2,3 > ENCFF547MET.bed
			#download_file "https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm405-db20140131/mm10.fa.out.gz" "$out_genome/mm10_repeatmasker.bed.gz"
			#gunzip -c mm10.fa.out.gz | sed 's/ \+/\t/g' | sed 's/^\t//g' | tail -n +4 | cut -f 5,6,7 > mm10_repeatmasker
			download_file "https://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/gap.txt.gz" "$out_genome/gap.txt.gz"
			gunzip -c gap.txt.gz | cut -f 2,3,4 > mm10_N.bed
			;;
		mm9)
			download_file "https://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/mm9.chrom.sizes" "$out_genome/mm9.chrom.sizes"
			download_file "https://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz" "$out_genome/mm9_repeatmasker.bed.gz"
			gunzip -c mm9.fa.out.gz | sed 's/ \+/\t/g' | sed 's/^\t//g' | tail -n +4 | cut -f 5,6,7 > mm9_repeatmasker
			;;
		hg38)
			download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" "$out_genome/hg38.chrom.sizes"
            download_file "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz" "$out_genome/ENCFF356LFX.bed.gz"
            gunzip -c ENCFF356LFX.bed.gz |cut -f 1,2,3 > ENCFF356LFX.bed
            #download_file "https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz" "$out_genome/hg38_repeatmasker.bed.gz"
            #gunzip -c hg38.fa.out.gz | sed 's/ \+/\t/g' | sed 's/^\t//g' | tail -n +4 | cut -f 5,6,7 > hg38_repeatmasker
            download_file "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz" "$out_genome/gap.txt.gz"
            gunzip -c gap.txt.gz | cut -f 2,3,4 > hg38_N.bed
			;; 
		*)
            echo "Genome $genome is not supported.";
			echo "Try '$(cmd) -h' for more information.";
			exit 1
            ;;
		# Add other genomes similarly
	esac
}

check_chromosome_naming() {
    local bam_file="$1"

    # Extract first few chromosome names from the BAM header
    local chr_names=$(sambamba -q view -H "$bam_file" | awk '$1=="@SQ" {print $2}' | cut -d':' -f2 | head -10)

    # Check for 'chr' prefix
    if echo "$chr_names" | grep -q '^chr[0-9XYM]'; then
        echo "chr-prefixed"
    elif echo "$chr_names" | grep -q '^[0-9XYM]\+'; then
        echo "numeric"
    else
        echo "unknown"
    fi
}

##############
##	Usage	##
##############

usage()
{
	echo "$(basename "$0") -hi host_input.bam -hc host_chip.bam -vi virus_input.bam -vc virus_chip.bam -p peaks_host.narrowPeak -g hg19 -o output_dir [-v score] [-n 200] [-w 1000] [-ws 0.5] [-e exclusion.bed] [-t 2] [-c 10] [-x 10] [-s suffix]
Alpha version 1.0

Mandatory:
	-hi	Host input bam (without duplicate removal)
	-hc	Host ChiP bam (without duplicate removal)
	-vi	Virus input bam aligned against a single viral fasta genome (without duplicate removal)
	-vc	Virus Chip bam aligned against a single viral fasta genome (without duplicate removal)
	-p	Host peaks .txt (epic2) or .narrowPeak (macs2)
	-g	host genome name (eg. mm9, mm10, hg38, hg19)
	-o	output directory
	
Optional:
	-v	Value for positive peak selection (eg. score, pvalue, FDR, log2FoldChange) - default: score
	-n	Number of selected positive sites - default: 200
	-w	shifting window size of ratio (ChiP/input) calculation - default: 1000
	-ws	Shifting size - fraction of window size (w option) - float up to 1 is allowed - default: 0.5
	-e	exclusion region file bed format (chr	start	end) - Not Implemented
	-t	threads number - default: 2 
	-c	Expected FPK thershold to apply local lambda correction - default: 10
	-x	Percentile threshold for filtering the host input FPK distribution, setting two cutoff limits: the lower limit at x and the upper limit at 100 - x percentiles. Default value: 10.
	-s	Suffix for output files

	-h  show this help text"
		
	exit 2
}

#
#	Viral blacklisted bed
#

##############
##	Params	##
##############
dir=${PWD};

FULL_COMMAND="$0 $@"

#~ BASEDIR=$(dirname "$0");
BASEDIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# Gets the command name without path
cmd(){ echo $(basename $0); }

# Error message
error(){
    echo "$(cmd): invalid option -- '$1'";
    echo "Try '$(cmd) -h' for more information.";
    exit 1;
}

PARSED_ARGUMENTS=$(getopt -n $(basename $0) --alternative -o '' --longoptions hi:,hc:,vi:,vc:,h::,p:,g:,o:,v:,n:,w:,ws:,e:,t:,c:,x:,s: -- "$@")
if [[ $? -ne 0 ]]; then
    exit 1;
fi
eval set -- $PARSED_ARGUMENTS

# extract options and their arguments into variables.
no_args="true"
while true ;
do
  case "$1" in
	--hi)	hi="$2"		; shift	2	;;
	--hc)	hc="$2"		; shift	2	;;
	--vi)	vi="$2"		; shift	2	;;
	--vc)	vc="$2"		; shift	2	;;
	--p)	p="$2"		; shift	2	;;
	--g)	g="$2"		; shift	2	;;
	--o)	o="$2"		; shift	2	;;
	--v)	v="$2"		; shift	2	;;
	--n)	n="$2"		; shift	2	;;
	--w)	w="$2"		; shift	2	;;
	--ws)	ws="$2"		; shift	2	;;
	--e)	e="$2"		; shift	2	;;
	--t)	t="$2"		; shift	2	;;
	--c)	c="$2"		; shift	2	;;
	--x)	x="$2"		; shift	2	;;
	--s)	s="$2"		; shift	2	;;
    --)		shift; break ;;
   (*) usage
       exit 1
       ;;
  esac
  no_args="false"
done

[[ "$no_args" == "true" ]] && { echo "$(cmd): No options provided"; usage; exit 1; }

if [ -z "$hi" ] || [ -z "$hc" ] || [ -z "$vi" ] || [ -z "$vc" ] || [ -z "$p" ]  || [ -z "$g" ] || [ -z "$o" ]
then
   echo "$(cmd): Missing one or more mandatory options"
   echo "Try '$(cmd) -h' for more information.";
   exit 1
else
	hi=$(to_absolute_path "$hi")
	hc=$(to_absolute_path "$hc")
	vi=$(to_absolute_path "$vi")
	vc=$(to_absolute_path "$vc")
	p=$(to_absolute_path "$p")
	o=$(to_absolute_path "$o")
fi

if [ -z "$v" ]
then
	v="score"
fi

if [ -z "$n" ]
then
	n=200
else
	if ! [[ "$n" =~ ^[0-9]+$ ]]; then
		echo "\"$n\" is not a positive integer: n parameter must be a positive integer"
		echo "Try '$(cmd) -h' for more information.";
		exit 1
	fi
fi

if [ -z "$w" ]
then
	w=1000
else
	if ! [[ "$w" =~ ^[0-9]+$ ]]; then
		echo "\"$w\" is not a positive integer: w parameter must be a positive integer"
		echo "Try '$(cmd) -h' for more information.";
		exit 1
	fi
fi


if [ -z "$ws" ]
then
	ws=0.5
else
	if ! [[ $ws =~ ^0(\.[0-9]+)?$ || $ws == "1" ]]; then
		echo "$ws is not a valid floating-point number between 0 and 1";
		echo "Try '$(cmd) -h' for more information.";
		exit 1
	fi
fi

## Calculate the number of consecutive region to consider around a target region to not have any bases overlap to the target region (non-overlapped)
nonoverlap_step=$(echo "((($w / $ws) / $w) + 0.9999)" | bc -l)
nonoverlap_step_int=$(echo "$nonoverlap_step / 1" | bc)  # Truncate to an integer

if [ -z "$e" ]
then
	e=""
fi

if [ -z "$t" ]
then
	t=2
else
	if ! [[ "$t" =~ ^[0-9]+$ ]]; then
		echo "\"$t\" is not a positive integer: t parameter must be a positive integer"
		echo "Try '$(cmd) -h' for more information.";
		exit 1
	fi
fi

if [ -z "$c" ]
then
	c=10
fi

if [ -z "$x" ]
then
	x=10
fi

if [ -z "$s" ]
then
	s=""
else
	## replace all the spaces with _
	s=$(echo $s | sed 's/ /_/g')
	s="_${s}"
fi

#######################
##    LOG SETTING    ##
#######################
echo "~~~~~~~~~~~~~~~~~~~~";
echo "Log";
echo "~~~~~~~~~~~~~~~~~~~~";

outdir=${o};

mkdir -p $outdir;

DATE=$(date +"%Y-%m-%d_%H-%M-%S")

LOG_FILE="$outdir/$DATE${s}.log"

# Redirect stdout to both terminal and log
exec > >(tee -a "$LOG_FILE")

# Redirect stderr only to the log (not to the terminal)
exec 2>> "$LOG_FILE"

# Trap function: Ensure error message appears ONCE in terminal & log
trap 'status=$?; if [ $status -ne 0 ]; then 
    msg="ERROR: \"${last_command}\" command failed with exit code $status."
    echo "$msg"  # Send to terminal
fi' EXIT

# Track last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG


echo "Script started at $(date)"
echo "Command:"
echo "$FULL_COMMAND"
echo "Parameters:"

# Process the string using regex matching
while [[ $PARSED_ARGUMENTS =~ (--[a-zA-Z]+)[[:space:]]+\'([^\']*)\' ]]; do
    key="${BASH_REMATCH[1]}"       # Capture the key (e.g., --hi)
    value="${BASH_REMATCH[2]}"     # Capture the value (e.g., path inside quotes)
    echo "$key = $value"       # Format the output as required
    # Trim the processed part from the string
    PARSED_ARGUMENTS="${PARSED_ARGUMENTS#*"${BASH_REMATCH[0]}"}"
done

RSCRIPT=$(which Rscript 2>/dev/null)

if [[ -z "$RSCRIPT" ]]; then
    echo "Error: Rscript not found in PATH" >&2
    exit 1
fi

#####################
##    DOWNLOAD     ##
#####################

echo "~~~~~~~~~~~~~~~~~~~~";
echo "Download";
echo "~~~~~~~~~~~~~~~~~~~~";

out_genome=$outdir"/genome"

mkdir -p $out_genome;

cd $out_genome

download_genome_files "$g"

##
##	ADD A TEST IF EACH OF THE INPUT FILE EXISTS AND NOT EMPTY
##

cd ${dir}

##########
##	QC	##
##########
## Some statistics and associated plots on peaks file
echo "~~~~~~~~~~~~~~~~~~~~";
echo "QC";
echo "~~~~~~~~~~~~~~~~~~~~";

name=$(basename "${p%.*}");
ext="${p#*.}"
if [ "$ext" == "txt" ]	## epic2
then
	${RSCRIPT} ${BASEDIR}/plot_epic2_qc.r -i ${p} -s ${name} -o ${outdir}
elif [ "$ext" == "narrowPeak" ]	## macs2
then
	${RSCRIPT} ${BASEDIR}/plot_macs2_qc.r -i ${p} -s ${name} -o ${outdir}
	sed $'1i #chr\tstart\tend\tname\tscore\tstrand\tsignal\tpvalue\tqvalue\tpeak' ${p} > ${outdir}/macs2_peaks.bed
	p="${outdir}/macs2_peaks.bed"
else
	echo "\"$p\" does not have a correct file extention (\".txt\" or \".narrowPeak\"). The file must be unmodified epic2 or macs2 output."
	echo "Try '$(cmd) -h' for more information.";
	exit 1
fi

###################################################
## Some statistics and associated plots on bam file

#	Blacklisted
cd ${out_genome};
cat $(ls *.bed) | bedtools sort | bedtools merge > blacklisted.bed

cd ${dir}

if [[ "$ws" < 1 ]]
then
	#~ shift_size=$($w*$ws);
	# shift_size=$(printf '%.*f\n' 0 $(bc <<< "$w * $ws"));
	shift_size=$(awk "BEGIN {print int($w * $ws)}")
else
	shift_size=$ws;
fi

#################################
##	Check if paired or single end
##
paired=0
paired=$(sambamba view -h ${hi} | head -n 10000 | samtools view -c -f 1)
if [[ "$paired" < 1 ]]
then
	echo "Single end reads detected";
	mode="";
else
	echo "Paired end reads detected";
	mode="-p"
fi


##
##	Change chrom.size so it contain only chr1-22, X, Y
##
naming_convention=$(check_chromosome_naming "${hi}")

if [ "$naming_convention" == "chr-prefixed" ]
then
	echo "Chromosome naming convention: chr-prefixed"
	awk '$1 ~ /^chr[0-9XY]+$/' ${out_genome}/${g}.chrom.sizes > ${out_genome}/${g}.chrom.sizes.1
elif [ "$naming_convention" == "numeric" ]
then
	echo "Chromosome naming convention: numeric"
	awk '$1 ~ /^chr[0-9XY]+$/' ${out_genome}/${g}.chrom.sizes | awk '{sub(/^chr/, ""); print}' > ${out_genome}/${g}.chrom.sizes.1
else
	echo "Chromosome naming convention not recognized: $naming_convention"
	echo "Try '$(cmd) -h' for more information.";
	exit 1
fi

rm ${out_genome}/${g}.chrom.sizes
mv ${out_genome}/${g}.chrom.sizes.1 ${out_genome}/${g}.chrom.sizes

#########################################################################################
##	Split host genome in same size windows (no blacklist) and count cov input and chipsig

bedtools makewindows -g ${out_genome}/${g}.chrom.sizes -w $w -s $shift_size | bedtools intersect -v -a /dev/stdin -b ${out_genome}/blacklisted.bed | awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' /dev/stdin > ${out_genome}/genome_win.saf

featureCounts ${mode} -T ${t} -O -F SAF --fracOverlap 0.5 -a ${out_genome}/genome_win.saf -o ${outdir}/genome_win_count.tsv ${hi} ${hc}

tail -n +3 ${outdir}/genome_win_count.tsv | awk -v OFS='\t' '{$9 = sprintf("%.3f", $7 / ( $6 / 1000 ) )}1' > ${outdir}/genome_win_count2.tsv

awk -v OFS='\t' '{$10 = sprintf("%.3f", $8 / ( $6 / 1000 ) )}1' ${outdir}/genome_win_count2.tsv > ${outdir}/genome_win_count.tsv

rm ${outdir}/genome_win_count2.tsv

host_win_count=$(wc -l ${outdir}/genome_win_count.tsv | cut -f 1 -d ' ')

median_cov_host_input=$(cut -f 7 ${outdir}/genome_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_host_input=$(cut -f 9 ${outdir}/genome_win_count.tsv| sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median input cov ${median_cov_host_input} - Median input fpk ${median_fpk_host_input}";

median_cov_host_chip=$(cut -f 8 ${outdir}/genome_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_host_chip=$(cut -f 10 ${outdir}/genome_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median chip cov ${median_cov_host_chip} - Median chip fpk ${median_fpk_host_chip}";

##########################################################################################
##	Split viral genome in same size windows (no blacklist) and count cov input and chipsig

## Add control if empty --> BAM file header is not correctly formatted, refers to help, exit 1
genome_name=$(sambamba -q view -H ${vi} | grep "SQ" | cut -d ":" -f 2 | cut -f 1);
genome_size=$(sambamba -q view -H ${vi} | grep "SQ" | cut -d ":" -f 3)

# echo -e "${genome_name}\t${genome_size}"

echo -e "${genome_name}\t${genome_size}" > "${out_genome}/${genome_name}.chrom.sizes"

bedtools makewindows -g ${out_genome}/${genome_name}.chrom.sizes -w $w -s $shift_size | awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' /dev/stdin > "${out_genome}/${genome_name}_win.saf"

featureCounts ${mode} -O -T ${t} -F SAF --fracOverlap 0.5 -a ${out_genome}/${genome_name}_win.saf -o ${outdir}/${genome_name}_win_count.tsv ${vi} ${vc}

tail -n +3 ${outdir}/${genome_name}_win_count.tsv | awk -v OFS='\t' '{$9 = sprintf("%.3f", $7 / ( $6 / 1000 ) )}1' > "${outdir}/${genome_name}_win_count2.tsv"

awk -v OFS='\t' '{$10 = sprintf("%.3f", $8 / ( $6 / 1000 ) )}1' ${outdir}/${genome_name}_win_count2.tsv > "${outdir}/${genome_name}_win_count.tsv"

rm ${outdir}/${genome_name}_win_count2.tsv

virus_win_count=$(wc -l ${outdir}/${genome_name}_win_count.tsv | cut -f 1 -d ' ')

median_cov_virus_input=$(cut -f 7 ${outdir}/${genome_name}_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_virus_input=$(cut -f 9 ${outdir}/${genome_name}_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median input cov virus ${median_cov_virus_input} - Median input fpk virus ${median_fpk_virus_input}";

median_cov_virus_chip=$(cut -f 8 ${outdir}/${genome_name}_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_virus_chip=$(cut -f 10 ${outdir}/${genome_name}_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median chip cov virus ${median_cov_virus_chip} - Median chip fpk virus ${median_fpk_virus_chip}";

##########################################################################################
##	Calculate cov and FPK on input ChiP peaks

## Add an extra step that remove chipseq peaks falling in blacklisted region of the genome (at least 10% of the peaks length)
echo -e "Checking for potential input host peaks in blacklisted region...";

bedtools intersect -v -a ${p} -b ${out_genome}/blacklisted.bed -wa > ${outdir}/peaks_blacklist_exc_tmp.bed

bedtools intersect -a ${p} -b ${out_genome}/blacklisted.bed -wa -u > ${outdir}/peaks_excluded.bed

if [ "$naming_convention" == "chr-prefixed" ]
then
	awk '$1 ~ /^chr([0-9]+|X|Y)$/' ${outdir}/peaks_blacklist_exc_tmp.bed > ${outdir}/peaks_blacklist_exc.bed
	awk '$1 !~ /^chr([0-9]+|X|Y)$/' ${outdir}/peaks_blacklist_exc_tmp.bed >> ${outdir}/peaks_excluded.bed
else
	awk '$1 ~ /^([0-9]+|X|Y)$/' ${outdir}/peaks_blacklist_exc_tmp.bed > ${outdir}/peaks_blacklist_exc.bed
	awk '$1 !~ /^([0-9]+|X|Y)$/' ${outdir}/peaks_blacklist_exc_tmp.bed >> ${outdir}/peaks_excluded.bed
fi

rm -f ${outdir}/peaks_blacklist_exc_tmp.bed

excluded_peaks=$(wc -l ${outdir}/peaks_excluded.bed | cut -f 1 -d ' ')

echo -e "${excluded_peaks} overlapping with blacklisted regions of the genome or present in non-canonical chromosomes; the peaks excluded can be retrieved in the file peaks_excluded.bed";

echo -e " ";

#~ echo -e "${column_number}";

#~ cut -f 1-3,${column_number} ${outdir}/peaks_blacklist_exc.bed | awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' /dev/stdin > ${outdir}/host_peaks.saf

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${outdir}/peaks_blacklist_exc.bed > ${outdir}/host_peaks.saf

featureCounts ${mode} -O -T ${t} -F SAF --fracOverlap 0.5 -a ${outdir}/host_peaks.saf -o ${outdir}/host_peaks_count.tsv ${hi} ${hc}

tail -n +3 ${outdir}/host_peaks_count.tsv | awk -v OFS='\t' '{$9 = sprintf("%.3f", $7 / ( $6 / 1000 ) )}1' > ${outdir}/host_peaks_count2.tsv

awk -v OFS='\t' '{$10 = sprintf("%.3f", $8 / ( $6 / 1000 ) )}1' ${outdir}/host_peaks_count2.tsv > ${outdir}/host_peaks_count.tsv

rm ${outdir}/host_peaks_count2.tsv

peak_count=$(wc -l ${outdir}/host_peaks_count.tsv | cut -f 1 -d ' ')

median_cov_peaks_input=$(cut -f 7 ${outdir}/host_peaks_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_peaks_input=$(cut -f 9 ${outdir}/host_peaks_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median input cov peaks ${median_cov_peaks_input} - Median input fpk peaks ${median_fpk_peaks_input}";

median_cov_peaks_chip=$(cut -f 8 ${outdir}/host_peaks_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
median_fpk_peaks_chip=$(cut -f 10 ${outdir}/host_peaks_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')
#~ echo "Median chip cov peaks ${median_cov_peaks_chip} - Median chip fpk peaks ${median_fpk_peaks_chip}";

## Calculate target genome copy number per host equivalent as ${median_fpk_virus_input}/${median_fpk_host_input}
estimated_target_genome_copy_number=$(echo "scale=2; ${median_fpk_virus_input}/${median_fpk_host_input}" | bc)

#############################
##	Save statistics on a file
stats_file="${outdir}/QC_stats${s}.tsv"
touch ${stats_file}
echo -e "\tHost - ${host_win_count} regions (${w} bp)\t\tTarget - ${virus_win_count} regions (${w} bp)\t\tPeaks - ${peak_count} peaks\t" > ${stats_file}
echo -e "Median\tCoverage\tFPK\tCoverage\tFPK\tCoverage\tFPK" >> ${stats_file}
echo -e "INPUT\t${median_cov_host_input}\t${median_fpk_host_input}\t${median_cov_virus_input}\t${median_fpk_virus_input}\t${median_cov_peaks_input}\t${median_fpk_peaks_input}" >> ${stats_file}
echo -e "ChiP\t${median_cov_host_chip}\t${median_fpk_host_chip}\t${median_cov_virus_chip}\t${median_fpk_virus_chip}\t${median_cov_peaks_chip}\t${median_fpk_peaks_chip}" >> ${stats_file}

column -t ${stats_file} -o " | " -s $'\t'

echo -e " ";

echo -e "Estimated target genome copy number per host genome equivalent ${estimated_target_genome_copy_number}";

echo -e " ";
##
##	As an extra QC, add the input FPK of the peaks and relate it to the overall FPK value (in both cases without the bins/peaks with FPK equal to 0)
##

echo "~~~~~~~~~~~~~~~~~~~~";
echo "Genome Bins";
echo "~~~~~~~~~~~~~~~~~~~~";

###############################################################
##	Alerte if median Coverage/FPK too low in input for bam file
lambda_input=0;

if (( $(echo "$median_fpk_host_input < $c" |bc -l) )); then
	echo "The median Fragment Per Kilobase (FPK) of the input bam file from the host is less than ${c} (${median_fpk_host_input})."
	echo "Applying local lambda calculation for host input FPK calculation.";
	lambda_input=1;
fi

echo -e " ";

if (( $(echo "$median_fpk_virus_input < $c" |bc -l) )); then
	echo "The median Fragment Per Kilobase (FPK) of the input bam file from the target is less than ${c} (${median_fpk_virus_input})."
	echo "Applying local lambda calculation for target input FPK calculation.";
	lambda_input=1;
fi

echo -e " ";
##########################################################################
##	Apply local lambda to host/target input coverage and FPK - Perl script

if (( $(echo "$lambda_input > 0" |bc -l) )); then
	## Host - genome_win_count.tsv must be sorted by chromosome and start position
	perl ${BASEDIR}/apply_local_lambda.pl -i ${outdir}/genome_win_count.tsv -l $median_fpk_host_input -s ${nonoverlap_step_int} -c $c -o ${outdir}/genome_win_count_lambda_corrected.tsv

	## Target - ${genome_name}_win_count.tsv must be sorted by chromosome and start position
	perl ${BASEDIR}/apply_local_lambda.pl -i ${outdir}/${genome_name}_win_count.tsv -l $median_fpk_virus_input -s ${nonoverlap_step_int} -c $c -o "${outdir}/${genome_name}_win_count_lambda_corrected${s}.tsv"


	median_fpk_host_input_lambda=$(cut -f 12 ${outdir}/genome_win_count_lambda_corrected.tsv| sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')

	median_fpk_virus_input_lambda=$(cut -f 12 "${outdir}/${genome_name}_win_count_lambda_corrected${s}.tsv" | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')

	# rounded_median_fpk_host_input_lambda=$(printf "%.3f" "$median_fpk_host_input_lambda")
	rounded_median_fpk_host_input_lambda=$(awk "BEGIN {printf \"%.3f\", $median_fpk_host_input_lambda}")

	# rounded_median_fpk_virus_input_lambda=$(printf "%.3f" "$median_fpk_virus_input_lambda")
	rounded_median_fpk_virus_input_lambda=$(awk "BEGIN {printf \"%.3f\", $median_fpk_virus_input_lambda}")

	echo "The newly calculated median Fragment Per Kilobase (FPK) of the input bam file from the host is ${rounded_median_fpk_host_input_lambda} (previous value ${median_fpk_host_input}).";
	echo "The newly calculated median Fragment Per Kilobase (FPK) of the input bam file from the target is ${rounded_median_fpk_virus_input_lambda} (previous value ${median_fpk_virus_input}).";

	echo -e " ";

	${RSCRIPT} ${BASEDIR}/ChiP_Input_FPK_QC.R ${outdir}/genome_win_count_lambda_corrected.tsv 1 ${outdir} "genome_FPK${s}.pdf"

	cut -f 2- ${outdir}/host_peaks_count.tsv > ${outdir}/host_peaks_count.bed
	cut -f 2- ${outdir}/genome_win_count_lambda_corrected.tsv | sed '1d' > ${outdir}/genome_win_count_lambda_corrected.bed

	## First column are entry in A and last column are overlapping entry in B
	bedtools intersect -a ${outdir}/host_peaks_count.bed -b ${outdir}/genome_win_count_lambda_corrected.bed -f 0.75 -F 0.75 -e -wa -wb > ${outdir}/intersect_peaks_lambda_corrected.bed
	
	cut -f 10,11,12,13,14,15,16,17,18,19,20 ${outdir}/intersect_peaks_lambda_corrected.bed > "${outdir}/positives_win_count_lambda_corrected.tsv"
	
	# Only report those entries in A that have _no overlaps_ with B
	bedtools intersect -v -a ${outdir}/genome_win_count_lambda_corrected.bed -b ${outdir}/host_peaks_count.bed > "${outdir}/negatives_win_count_lambda_corrected.tsv"
	
	temp_file=$(mktemp)
	cut -f 2- "${outdir}/${genome_name}_win_count_lambda_corrected${s}.tsv" | sed '1d' > "${temp_file}"
 	mv "$temp_file" "${outdir}/${genome_name}_win_count_lambda_corrected${s}.tsv"


	# The three files are formatted as follows:
	# Chromosome	Start	End	Strand	Length	CovInput	CovChiP	FPKInput	FPKChiP	CovInputCorrected	FPKInputCorrected
	## Here the 3 distribution to give to R are:
	## 1) positives_win_count_lambda_corrected.tsv
	## 2) negatives_win_count_lambda_corrected.tsv
	## 3) ${genome_name}_win_count_lambda_corrected${s}.tsv

else
	## Need to add header to genome_win_count.tsv
	# SAF	Chromosome	Start	End	Strand	Length	CovInput	CovChiP	FPKInput	FPKChiP

	sed $'1i SAF\tChromosome\tStart\tEnd\tStrand\tLength\tCovInput\tCovChiP\tFPKInput\tFPKChiP' ${outdir}/genome_win_count.tsv > ${outdir}/genome_win_count_header.tsv

	cut -f 2- ${outdir}/${genome_name}_win_count.tsv > "${outdir}/filtered_${genome_name}_win_count${s}.tsv"
	
	${RSCRIPT} ${BASEDIR}/ChiP_Input_FPK_QC.R ${outdir}/genome_win_count_header.tsv 0 ${outdir} "genome_FPK${s}.pdf"
	
	cut -f 2- ${outdir}/host_peaks_count.tsv > ${outdir}/host_peaks_count.bed
	cut -f 2- ${outdir}/genome_win_count.tsv > ${outdir}/genome_win_count.bed

	## First column are entry in A and last column are overlapping entry in B
	bedtools intersect -a ${outdir}/host_peaks_count.bed -b ${outdir}/genome_win_count.bed -wa -wb > ${outdir}/intersect_peaks.bed

	cut -f 10,11,12,13,14,15,16,17,18 ${outdir}/intersect_peaks.bed > "${outdir}/positives_win_count.tsv"

	bedtools intersect -v -a ${outdir}/genome_win_count.bed -b ${outdir}/host_peaks_count.bed > "${outdir}/negatives_win_count.tsv"
	
	## Here the 3 distribution to give to R are:
	## 1) positives_win_count.tsv
	## 2) negatives_win_count.tsv
	## 3) filtered_${genome_name}_win_count${s}.tsv

fi

##########################################################################
##	Host sites selection
##  select_positives_and_negatives_regions_host.pl --> readapted with shell commands (at the end to be removed from ViPeaQ)

echo "~~~~~~~~~~~~~~~~~~~~";
echo "Top positives and negatives";
echo "~~~~~~~~~~~~~~~~~~~~";

## here get the column number matching the v parameter
# Enable case-insensitive matching
shopt -s nocasematch

# Extract the header (first line of the file)
header=$(head -n 1 "$p")

# Convert the header to an array of column names
IFS=$'\t' read -r -a columns <<< "$header"

# Initialize a variable to store the column number
column_number=-1

# Find the column number by comparing each column name to the target regex
for i in "${!columns[@]}"; do
	if [[ "${columns[$i]}" =~ $v ]]; then
		column_number=$((i + 1))
		break
	fi
done

# Disable case-insensitive matching
shopt -u nocasematch

# Check if the column was found
if [[ $column_number -eq -1 ]]; then
	echo "The file $p does not have a header string matching '$v'.";
	echo "Try '$(cmd) -h' for more information.";
	exit 1
fi

sorted_counts=$(mktemp)
sorted_peaks=$(mktemp)
 
sort -k1,1 "${outdir}/host_peaks_count.tsv" > "$sorted_counts"
awk 'OFS="\t" {print $1"."$2"."$3, $0}' "${outdir}/peaks_blacklist_exc.bed" | sort -k1,1 > "$sorted_peaks"

column_number=$(($column_number+1))

# Determine the number of columns in sorted_counts
num_columns_sorted_counts=$(awk '{print NF; exit}' "$sorted_counts")

# Construct the -o option string
o_option="1.2"
for ((i=3; i<=num_columns_sorted_counts; i++)); do
  o_option+=",1.$i"
done
o_option+=",2.$column_number"

join -1 1 -2 1 -o $o_option -t $'\t' "$sorted_counts" "$sorted_peaks" > "${outdir}/host_peaks_count.tsv"

#~ # Clean up temporary files
rm "$sorted_counts" "$sorted_peaks"


######################################################
##	No filter on fpk input of the peaks (for now)	##
######################################################
lambda_peak=0

if (( $(echo "$lambda_peak > 0" |bc -l) )); then
	## remove peaks with corrected fpk (num_columns_sorted_counts-1) below $rounded_median_fpk_peaks_input_lambda
	columns_corrected_fpk=$(($num_columns_sorted_counts-1))
	#~ threshold_fpk=$rounded_median_fpk_peaks_input_lambda
	threshold_fpk=1
else
	## remove peaks with fpk (num_columns_sorted_counts-2) below 10
	columns_corrected_fpk=$(($num_columns_sorted_counts-2))
	#~ threshold_fpk=10	
	threshold_fpk=0
fi

##
##	Here threshold FPK must be:
##	calculate the fraction of bins over the all genome that have a FPK value in the input = 0 --> flag if high
##	calculate the median FPK value of the bin inc. the ones equal to 0, and excluding the one equal to 0
##	From the distribution of the FPK value in the input when excluding the one equal to 0, calculate the threshold values as eg. the lowest and highest 10th percentile
##	Exclude all positives peaks that have a median FPK value below or above the lowest and highest 10h percentile thershold
##	Take the top n peaks based on the eg. score (user-defined)
##	report the median FPK value of the top n peaks and compare it to the median FPK value of all the bins in the genome (excluding the one equal to 0)
##

zero_lines=$(awk -F'\t' '$9 == 0 {count++} END {print count}' ${outdir}/genome_win_count.tsv)
fraction=$(awk "BEGIN {print $zero_lines / $host_win_count}")
echo "Fraction of genome bins with a FPK value equal to 0: $fraction"


median_fpk_host_input_no_zero=$(awk -F'\t' '$9 != 0 {print $9}' ${outdir}/genome_win_count.tsv | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')

echo "Median input FPK value of all genome bins: ${median_fpk_host_input}"
echo "Median input FPK value of all genome bins with a FPK value different from 0: ${median_fpk_host_input_no_zero}"

# Filter non-zero values, sort them numerically
values=$(awk -F'\t' '$9 != 0 {print $9}' ${outdir}/genome_win_count.tsv | sort -n)

# Total number of non-zero values
count=$(echo "$values" | wc -l)

# Calculate positions for the lowest and highest percentiles --> I should do the same for the bin
low_index=$(awk "BEGIN {printf(\"%.0f\", ($x / 100) * $count)}")
high_index=$(awk "BEGIN {printf(\"%.0f\", $count - (($x / 100) * $count) + 1)}")

# Extract the values
low_percentile=$(echo "$values" | sed -n "${low_index}p")
high_percentile=$(echo "$values" | sed -n "${high_index}p")

echo "The $x-th lowest percentile FPK value is: $low_percentile"
echo "The $x-th highest percentile FPK value is: $high_percentile"

echo "The host peaks will be filtered to excluded peaks with input FPK below $low_percentile and above $high_percentile."

initial_peaks=$(wc -l ${outdir}/host_peaks_count.tsv | cut -f 1 -d ' ')

awk -F $'\t' -v col="$columns_corrected_fpk" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/host_peaks_count.tsv > ${outdir}/host_peaks_count_filtered.tsv

sort -t $'\t' -k"$num_columns_sorted_counts","$num_columns_sorted_counts"nr ${outdir}/host_peaks_count_filtered.tsv > ${outdir}/host_peaks_count_filtered_sorted.tsv

remain_peaks=$(wc -l ${outdir}/host_peaks_count_filtered_sorted.tsv | cut -f 1 -d ' ')

excluded_peaks=$(( $initial_peaks - $remain_peaks ))

excluded_percentage=$(echo "scale=2; $excluded_peaks / $initial_peaks * 100" | bc)

echo "The number of excluded peaks is $excluded_peaks out of $initial_peaks ($excluded_percentage%)."

if (( $n > $remain_peaks )); then
	if (( $remain_peaks > 0 )); then
		head -n "$remain_peaks" ${outdir}/host_peaks_count_filtered_sorted.tsv > "${outdir}/top_positives_peaks${s}.tsv"
		echo "Following the filtering of the peaks, less than $n peaks are remaining. The calculation will carry on with $remain_peaks top peaks."
	else
		echo "No top peaks can be selected. Please check the input peaks file."
		echo "Try '$(cmd) -h' for more information.";
		exit 1;
	fi
else
	head -n "$n" ${outdir}/host_peaks_count_filtered_sorted.tsv > "${outdir}/top_positives_peaks${s}.tsv"
fi

median_fpk_host_input_top_positives=$(awk -F'\t' -v col="$columns_corrected_fpk" '$col != 0 {print $col}' "${outdir}/top_positives_peaks${s}.tsv" | sort -n | awk 'NF {a[NR] = $1} END {print (NR % 2 ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2)}')

negatives_regions=0

if (( $n > $remain_peaks )); then
	echo "Median input FPK value of the top $remain_peaks positives peaks is: ${median_fpk_host_input_top_positives}"
	negatives_regions=$remain_peaks
else
	echo "Median input FPK value of the top $n positives peaks is: ${median_fpk_host_input_top_positives}"
	negatives_regions=$n
fi

cut -f1,2,3 ${out_genome}/blacklisted.bed > ${out_genome}/exclusion_file.bed

cut -f1,2,3 "${outdir}/top_positives_peaks${s}.tsv" >> ${out_genome}/exclusion_file.bed

cut -f1,2,3 "${outdir}/top_positives_peaks${s}.tsv" > ${outdir}/top_positives_peaks.bed

exclusion_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+= $3 - $2 }END{print SUM}' ${out_genome}/exclusion_file.bed)

genome_size=$(samtools view -H ${hi} | grep '^@SQ' | awk '{sum += substr($3, 4)} END {print sum}')

available_size=$(($genome_size-$exclusion_size))

i=0;
fpk_column=9

while (( $i < $negatives_regions )); do
	bedtools shuffle -maxTries 1000 -chrom -noOverlapping -excl ${out_genome}/exclusion_file.bed -i ${outdir}/top_positives_peaks.bed -g ${out_genome}/${g}.chrom.sizes > ${outdir}/negatives_peaks_tmp.bed
	
	cat ${outdir}/negatives_peaks_tmp.bed >> ${outdir}/negatives_peaks.bed
	cat ${outdir}/negatives_peaks_tmp.bed >> ${out_genome}/exclusion_file.bed

	awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${outdir}/negatives_peaks_tmp.bed > ${outdir}/negatives_peaks.saf

	featureCounts ${mode} -O -T ${t} -F SAF --fracOverlap 0.5 -a ${outdir}/negatives_peaks.saf -o ${outdir}/negative_sites_host_count.tsv ${hi} ${hc}

	tail -n +3 ${outdir}/negative_sites_host_count.tsv | awk -v OFS='\t' '{$9 = sprintf("%.3f", $7 / ( $6 / 1000 ) )}1' > ${outdir}/negative_sites_host_count2.tsv

	awk -v OFS='\t' '{$10 = sprintf("%.3f", $8 / ( $6 / 1000 ) )}1' ${outdir}/negative_sites_host_count2.tsv > ${outdir}/negative_sites_host_count.tsv

	rm ${outdir}/negative_sites_host_count2.tsv

	awk -F $'\t' -v col="$fpk_column" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/negative_sites_host_count.tsv >> ${outdir}/negative_sites_host_count_filtered.tsv

	nb=$(wc -l ${outdir}/negative_sites_host_count_filtered.tsv | awk '{print $1}')

	rm -f ${outdir}/negative_sites_host_count.tsv
	rm -f ${outdir}/negative_sites_host_count.tsv.summary

	i=$nb
done

## Negatives
remain_peaks_neg=$(wc -l ${outdir}/negative_sites_host_count_filtered.tsv | cut -f 1 -d ' ')

if (( $n > $remain_peaks_neg )); then
	if (( $remain_peaks_neg > 0 )); then
		head -n "$remain_peaks_neg" ${outdir}/negative_sites_host_count_filtered.tsv | cut -f 2- > "${outdir}/top_negatives_peaks${s}.tsv"
		echo "Following the searching of negatives peaks of same chromosomes and size distribution of the positives peaks, only $n peaks could be found."
	else
		echo "No negatives peaks can be selected. Please check the blacklisted file."
		echo "Try '$(cmd) -h' for more information.";
		exit 1;
	fi
else
	head -n "$n" ${outdir}/negative_sites_host_count_filtered.tsv | cut -f 2- > "${outdir}/top_negatives_peaks${s}.tsv"
fi



## Positives
# Chromosome	Start	End	Strand	Length	CovInput	CovChiP	FPKInput	FPKChiP	Score
# ${outdir}/top_positives_peaks.tsv
## Negatives
# Chromosome	Start	End	Strand	Length	CovInput	CovChiP	FPKInput	FPKChiP
# ${outdir}/top_negatives_peaks.tsv

echo "~~~~~~~~~~~~~~~~~~~~";
echo "Output writing";
echo "~~~~~~~~~~~~~~~~~~~~";

if (( $(echo "$lambda_input > 0" |bc -l) )); then
	## Here the 3 distribution to give to R are:
	# Chromosome	Start	End	Strand	Length	CovInput	CovChiP	FPKInput	FPKChiP	CovInputCorrected	FPKInputCorrected
	## 1) positives_win_count_lambda_corrected.tsv
	## 2) negatives_win_count_lambda_corrected.tsv
	## 3) ${genome_name}_win_count_lambda_corrected.tsv

	input_fpk_col=11
	awk -F $'\t' -v col="$input_fpk_col" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/positives_win_count_lambda_corrected.tsv > "${outdir}/positives_win_count_lambda_corrected_filtered${s}.tsv"
	awk -F $'\t' -v col="$input_fpk_col" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/negatives_win_count_lambda_corrected.tsv > "${outdir}/negatives_win_count_lambda_corrected_filtered${s}.tsv"

	${RSCRIPT} ${BASEDIR}/ChiP_statistics_all.R \
	"${outdir}/top_positives_peaks${s}.tsv" \
	"${outdir}/top_negatives_peaks${s}.tsv" \
	"${outdir}/${genome_name}_win_count_lambda_corrected${s}.tsv" \
	"${outdir}/positives_win_count_lambda_corrected_filtered${s}.tsv" \
	"${outdir}/negatives_win_count_lambda_corrected_filtered${s}.tsv" \
	"${outdir}" \
	"${lambda_input}" \
	"${s}"

else

	input_fpk_col=9
	awk -F $'\t' -v col="$input_fpk_col" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/positives_win_count.tsv > "${outdir}/positives_win_count_filtered${s}.tsv"
	awk -F $'\t' -v col="$input_fpk_col" -v low="$low_percentile" -v high="$high_percentile" '$col > low && $col < high' ${outdir}/negatives_win_count.tsv > "${outdir}/negatives_win_count_filtered${s}.tsv"

	## Here the 3 distribution to give to R are:
	## 1) positives_win_count.tsv
	## 2) negatives_win_count.tsv
	## 3) filtered_${genome_name}_win_count.tsv --> for now it is unfiltered on FPK value input

	${RSCRIPT} ${BASEDIR}/ChiP_statistics_all.R \
	"${outdir}/top_positives_peaks${s}.tsv" \
	"${outdir}/top_negatives_peaks${s}.tsv" \
	"${outdir}/filtered_${genome_name}_win_count${s}.tsv" \
	"${outdir}/positives_win_count_filtered${s}.tsv" \
	"${outdir}/negatives_win_count_filtered${s}.tsv" \
	"${outdir}" \
	"${lambda_input}" \
	"${s}"

fi

##################
##	Cleaning	##
##################
## Some statistics and associated plots on peaks file
echo "~~~~~~~~~~~~~~~~~~~~";
echo "Cleaning";
echo "~~~~~~~~~~~~~~~~~~~~";

rm -f ${outdir}/macs2_peaks.bed
rm -f ${out_genome}/blacklisted.bed
rm -f ${out_genome}/genome_win.saf
rm -f ${out_genome}/genome_win_count.tsv
rm -f ${outdir}/genome_win_count.tsv.summary
rm -f ${out_genome}/${genome_name}.chrom.sizes
rm -f ${out_genome}/${genome_name}_win.saf
rm -f ${outdir}/${genome_name}_win_count.tsv
rm -f ${outdir}/${genome_name}_win_count.tsv.summary
rm -f ${outdir}/peaks_blacklist_exc.bed
rm -f ${outdir}/host_peaks.saf
rm -f ${outdir}/host_peaks_count.tsv
rm -f ${outdir}/host_peaks_count.tsv.summary
rm -f ${outdir}/host_peaks_count.bed
rm -f ${outdir}/genome_win_count_lambda_corrected.bed
rm -f ${outdir}/genome_win_count_lambda_corrected.tsv
rm -f ${outdir}/intersect_peaks_lambda_corrected.bed
rm -f ${outdir}/genome_win_count_header.tsv
rm -f ${outdir}/genome_win_count.bed
rm -f ${outdir}/intersect_peaks.bed
rm -f ${outdir}/host_peaks_count.tsv
rm -f ${outdir}/genome_win_count.tsv
rm -f ${outdir}/host_peaks_count_filtered.tsv
rm -f ${outdir}/host_peaks_count_filtered_sorted.tsv
rm -f ${out_genome}/exclusion_file.bed
rm -f ${outdir}/negatives_peaks_tmp.bed
rm -f ${outdir}/negatives_peaks.bed
rm -f ${outdir}/negatives_peaks.saf
rm -f ${outdir}/negative_sites_host_count_filtered.tsv
rm -f ${outdir}/positives_win_count_lambda_corrected.tsv
rm -f ${outdir}/negatives_win_count_lambda_corrected.tsv
rm -f ${outdir}/positives_win_count.tsv
rm -f ${outdir}/negatives_win_count.tsv
rm -f ${outdir}/top_positives_peaks.bed
rm -f ${outdir}/macs2_peakqc.summary.txt
rm -f ${outdir}/macs2_peakqc.plots.pdf

exit 0