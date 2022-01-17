#!/bin/bash

# ChIA-PIPE
#         Step 3: Perform peak calling
# 2018
# The Jackson Laboratory for Genomic Medicine


# The help message:
function usage
{
    echo -e "usage: bash 3.call_peaks.sh --conf ${conf}
    " 
}

# Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# Source config file
source ${conf}

# Add dependency dir to path
export PATH=${dep_dir}:${PATH}

# Create name of the log file
log_file=3.${run}.call_peaks.log

# Print arguments to ensure correct parsing
echo "
Arguments:
    run=${run}
    peak_caller=${peak_caller}
    input_control=${input_control}
    bin_dir=${bin_dir}
    out_dir=${out_dir}
" >> ${log_file}


## Perform ChIP-Seq peak calling and plot density from ChIA-PET data

# Convert to bedgraph
# Sort bam for samtools counting bases and for BASIC visualization
if [ ! -f ${run}.for.BROWSER.bam ]
then
    echo -e "`date` Converting file formats..\n" >> ${log_file}
    samtools sort -@ 16 ${run}.singlelinker.paired.UU.nr.bam \
        -o ${run}.singlelinker.paired.UU.nr.sorted.bam
    samtools sort -@ 16 ${run}.singlelinker.single.UxxU.nr.bam \
        -o ${run}.singlelinker.single.UxxU.nr.sorted.bam
    samtools sort -@ 16 ${run}.none.UU.nr.bam \
        -o ${run}.none.UU.nr.sorted.bam

    samtools merge ${run}.for.BROWSER.bam \
        ${run}.singlelinker.paired.UU.nr.sorted.bam \
        ${run}.singlelinker.single.UxxU.nr.sorted.bam \
        ${run}.none.UU.nr.sorted.bam
    
    # Create BAM index
    samtools index ${run}.for.BROWSER.bam ${run}.for.BROWSER.bam.bai
    
    # Make bedgraph
    bedtools genomecov -ibam ${run}.for.BROWSER.bam \
        -bg > ${run}.for.BROWSER.bedgraph
    
    echo -e "`date` Sorting bedgraph..\n" >> ${log_file}
    # Sort bedgraph
    ${bin_dir}/util/scripts/bedSort \
        ${run}.for.BROWSER.bedgraph \
        ${run}.for.BROWSER.sorted.bedgraph
    
    # Make bigwig
    ${bin_dir}/util/scripts/bedGraphToBigWig \
        ${run}.for.BROWSER.sorted.bedgraph \
        ${chrom_sizes} \
        ${run}.for.BROWSER.bigwig    
fi


if [ ${clean} == true ]
then
    ## Remove redundant BAMs
    # BAMs of reads with no linker
    rm -f ${run}.none.*.bam

    # BAMs of reads with only one usable tag
    rm -f ${run}.singlelinker.single.*.bam

    # For reads with two tags, save only the final
    # deduplicated and sorted BAM file
    ls ${run}.singlelinker.paired.*.bam | grep -v "sorted" | xargs rm -f

    # Remove FASTQ files of linker filtering
    rm -f ${run}.*.fastq.gz

    # Remove SAM files
    rm -f ${run}.*.sam.gz
    
    # Remove unsorted bedgraph
    rm -f ${run}.for.BROWSER.bedgraph
        
fi



if [ ${peak_caller} == "spp" ] || [ ${peak_caller} == "SPP" ]
then
    # Report SPP to log file
    echo -e "`date` Peak calling using SPP..\n" >> ${log_file}
    
    # Call peaks using SPP
    R --vanilla < ${bin_dir}/util/scripts/spp.R --args ${run}.for.BROWSER.bam \
        ${input_control} ${bin_dir} ${z_thresh}
else
    # Load MACS2
    
    # Report MACS2 to log file
    echo -e "`date` Peak calling using MACS2..\n" >> ${log_file}
    
    if [ ${input_control} == 'none' ] || [ ${input_control} == 'None' ]
    then
        # Call peaks using MACS2 without input control
        macs2 callpeak --keep-dup all --nomodel -t ${run}.for.BROWSER.bam \
            -f BAM -g hs -n ${run}.no_input_all 1>> ${log_file} 2>> ${log_file}
    else
        # Call peaks using MACS2 with input control
        macs2 callpeak --keep-dup all --nomodel -t ${run}.for.BROWSER.bam \
            -c ${input_control} \
            -f BAM -g hs -n ${run}.all 1>> ${log_file} 2>> ${log_file}
    fi
fi

echo -e "`date` ENDED ${run} peak calling ..\n" >> ${log_file}
echo "$0 done" >> ${log_file}


#### Annotate loops with peak support and call CCDs
# Annotate loops
be3_file="${run}.e500.clusters.cis.BE3"
peak_file="${run}.for.BROWSER.spp.z6.broadPeak"

${dep_dir}/python ${bin_dir}/util/scripts/annotate_loops_with_peak_support.py \
    -l ${be3_file} -p ${peak_file}

# Create subsetted files by peak support
# 2 anchors
cat ${be3_file}.peak_annot | awk '{ if ( $8 == 2 ) print }' \
    > ${be3_file}.peak_annot.E2

# 1 or more anchor
cat ${be3_file}.peak_annot | awk '{ if ( $8 >= 1 ) print }' \
    > ${be3_file}.peak_annot.BE1

# Convert loops to WashU format
# 2 anchors
${dep_dir}/python ${bin_dir}/util/scripts/convert_loops_to_washu_format.py \
    -l ${be3_file}.peak_annot.E2

# 1 or more anchor
${dep_dir}/python ${bin_dir}/util/scripts/convert_loops_to_washu_format.py \
    -l ${be3_file}.peak_annot.BE1


## Annotate enhancer-promoter loops
# Split loops into anchors
${dep_dir}/python ${bin_dir}/util/scripts/split_loops_into_anchors.py \
    -l ${be3_file}.peak_annot.BE1

left_anch=${be3_file}.peak_annot.BE1.left_anchors
right_anch=${be3_file}.peak_annot.BE1.right_anchors


# Intersect anchors with enhancers
bedtools intersect -u -a ${left_anch} -b ${enhancer_bed_file} \
    > ${left_anch}.enhancers

bedtools intersect -u -a ${right_anch} -b ${enhancer_bed_file} \
    > ${right_anch}.enhancers


# Intersect anchors with promoters
bedtools intersect -u -a ${left_anch} -b ${promoter_bed_file} \
    > ${left_anch}.promoters

bedtools intersect -u -a ${right_anch} -b ${promoter_bed_file} \
    > ${right_anch}.promoters

# Annotate enhancer promoter loops
${dep_dir}/python ${bin_dir}/util/scripts/annotate_enhancer_promoter_loops.py \
    --left_enhancers ${left_anch}.enhancers \
    --right_enhancers ${right_anch}.enhancers \
    --left_promoters ${left_anch}.promoters \
    --right_promoters ${right_anch}.promoters


####
# Sort peak-supported loops by PET count
sort -k7,7n ${be3_file}.peak_annot.E2 > ${be3_file}.peak_annot.E2.freq_sorted

## Get PET count cutoff for calling CCDs
## Top one-third (67th percentile)

# Line number in sorted file
cutoff_line=$( echo -e \
    "$( cat ${be3_file}.peak_annot.E2.freq_sorted | wc -l ) * 0.67 / 1" | bc )

# PET count cutoff
min_pet_count=$( sed "${cutoff_line}q;d" \
    ${be3_file}.peak_annot.E2.freq_sorted | awk '{ print $7 }' )

### Call CCDs

# Set input file for calling CCDs
ccd_input="${be3_file}.peak_annot.E2"

# Make bed file of loop spans
#awk -v OFS='\t' '{ print $1, $3, $5, $7 }' ${ccd_input} \
#    > ${ccd_input}.inner

${dep_dir}/python ${bin_dir}/util/scripts/get_loop_anchor_midpoints.py \
    -l ${ccd_input} \
    -m ${min_pet_count}


# Sort loops
sort -k1,1 -k2,2n ${ccd_input}.anchor_mids > ${ccd_input}.anchor_mids.sorted

# Merge loops
bedtools merge -i ${ccd_input}.anchor_mids.sorted \
    > ${ccd_input}.anchor_mids.sorted.merge

# Add ccd span
awk -v OFS='\t' '{ $4 = $3 - $2} 1' ${ccd_input}.anchor_mids.sorted.merge \
    > ${ccd_input}.anchor_mids.sorted.merge.span

# Filter out CCDs with span < 10kb
awk -v OFS='\t' '{ if ( $4 > 25000 ) print $0 }' \
    ${ccd_input}.anchor_mids.sorted.merge.span \
    > ${ccd_input}.anchor_mids.sorted.merge.span.ccd

# Remove CCD intermediate files
rm -f *.anchor_mids.sorted 
rm -f *.anchor_mids.sorted.merge
rm -f *.anchor_mids.sorted.merge.span


# Create a Juicebox 2D annotation file of CCDs
awk 'BEGIN{printf("chr1\tx1\tx2\tchr2\ty1\ty2\tcolor\n");} 
    {printf("%s\t%s\t%s\t%s\t%s\t%s\t0,204,255\n",$1,$2,$3,$1,$2,$3);}' \
    ${ccd_input}.anchor_mids.sorted.merge.span.ccd \
    > ${run}_ccd_juicebox_2d_annotation.txt

