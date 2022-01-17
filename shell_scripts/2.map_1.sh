#!/bin/bash

# Set the output to be redirected to lof file
LOG_LOCATION=/Chiapipe_mod/gm19238/
exec 19>$LOG_LOCATION/logfile_2step_paired.log
export BASH_XTRACEFD=19

set -x

# ChIA-PIPE
#         Step 2: Map the reads
# 2018
# The Jackson Laboratory for Genomic Medicine


## The help message:
function usage
{
    echo -e "usage: bash 2.map.sh --conf ${conf} --tag_name ${tag_name} 
    " 
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -c | --conf )           shift
                                conf=$1
                                ;;
        -t | --tag_name )       shift
                                tag_name=$1   
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


## Source config file
source ${conf}

# Add dependency dir to path
export PATH=${dep_dir}:${PATH}

# Set mapping parameters and file suffix
if [ ${tag_name} == "singlelinker.single" ]; then
    map_qual=10
    suffix="UxxU"
else
    map_qual=30
    suffix="UU"
fi

# Update mapping parameters for hybrid data
if [ ${hybrid} != "none" ]; then
    map_qual=0
fi

## Create name of the log file
log_file=2.${run}.map_${tag_name}.log



## Cluster tags and make files for visualization with Juicebox and HiGlass
if [ ${tag_name} == "singlelinker.paired" ]
then
    	
          
    ### Make .hic file for Juicebox
    ### BAM --> pairs --> hic
    # BAM to pairs
    paired_tag_bam="${run}.singlelinker.paired.UU.nr.bam"
    ${bin_dir}/util/pairix_src/util/bam2pairs/bam2pairs -c ${chrom_sizes} \
        ${paired_tag_bam} ${run}
    
    # Pairs to .hic
    juicer="${bin_dir}/util/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar"
    hic_file="ChIA-PET_${genome}_${cell_type}_${ip_factor}_${run}_${run_type}_pairs.hic"
    
    java -Xmx2g -jar ${juicer} pre -r \
        2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 \
        ${run}.bsorted.pairs.gz ${hic_file} ${chrom_sizes}
    
    
    # Copy .hic file to URL directory
    # for automatic viewing of 2D contact map
    # url_dir="/projects/encode/to_ctencode01/dev"
    
    # if [ -d ${url_dir} ]; then
    #     cp ${hic_file} ${url_dir}
    # fi
    
    ### Multi-resolution matrix file for HiGlass
    ## Minimum resolution for HiGlass
    #resolution=1000
    
    ## Divide chromosomes into bins
    #${bin_dir}/conda/bin/cooler makebins \
    #    ${chrom_sizes} \
    #    ${resolution} \
    #    > ${out_dir}/temp_chrom_bins.bed
    
    ## Make .cool file (1000bp resolution)
    #${bin_dir}/conda/bin/cooler cload pairix -p 20\
    #    ${out_dir}/tmp_chrom_bins.bed \
    #    ${out_dir}/${run}.e500.cooler.bsorted.pairs.gz \
    #    ${out_dir}/${run}.e500.higlass.cool
    
    ## Make final, multi-resolution .cool file (with zooming ability)
    #${bin_dir}/conda/bin/cooler zoomify -p 20\
    #    --no-balance ${out_dir}/${run}.e500.higlass.cool
    
    ## Remove intermediary file
    #rm ${out_dir}/temp_chrom_bins.bed
    #rm ${out_dir}/${run}.e500.higlass.cool

fi

# Report completion of mapping
echo "$0 done" >> 2.${run}.map_${tag_name}.done
echo "$0 done" >> ${log_file}
echo "`date`" >> ${log_file}
