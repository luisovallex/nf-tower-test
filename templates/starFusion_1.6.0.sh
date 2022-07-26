#!/bin/bash
#STAR_Fusion_1.6.0

source $baseDir/templates/common.sh
script_begin ${task.process} "$sample_name"

{
    echo "runid = ${runid}"
    echo "sample_name = ${sample_name}"
    echo "chimeric junction concat file: ${chimeric_junctions_merged}, metrics file : ${metrics}"

    # retrieve totalreads from  metricsfile
    totalreads="\$(cut -f 20 ${metrics} | tail -1)"
    echo "successfully retreived total reads"

    awk '\$1="chr"\$1, \$4="chr"\$4' ${chimeric_junctions_merged} | tr " " "\t" > sfv1.6.0.txt
    echo "# Nreads \${totalreads}" >> sfv1.6.0.txt

    $params.starfusion_new_path/STAR-Fusion \
                --STAR_PATH $params.star_new_path/bin/Linux_x86_64_static/STAR \
                --genome_lib_dir ${genome_ctat_lib} \
                -J sfv1.6.0.txt -O pga_fusion_results \
                --min_novel_junction_support 1 \
                --min_alt_pct_junction 1 \
                --examine_coding_effect \
                --min_spanning_frags_only 5 \
                --max_promiscuity 10 \
                --min_pct_dom_promiscuity 20 \
                --require_LDAS 0 \
                --no_annotation_filter \
                --min_FFPM 0;
    ev=\$?
    sleep 60

    

 } 2>&1

script_exit
