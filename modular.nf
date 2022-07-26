#!/usr/bin/env nextflow

import groovy.json.*
product = "envisia"
params.product = "afirmav2"
params.runid = ""
params.uuid = ""
params.runpath = ""

params.configLocation = "rnaa-config.json"
println "RNAA config is: " + params.configLocation
def configfile = new File(params.configLocation.toString())
def config = new JsonSlurper().parseText(configfile.text)
params.config = config
runid = params.runpath.split('/')[0]

params.pkgroot = "/usr/local/pipeline/packages"
params.star_path = "${params.pkgroot}/STAR_2.4.1b"
params.star_new_path = "${params.pkgroot}/STAR-2.7.1a"
params.javatmp = "/mnt/tmp"
params.picard_path = "${params.pkgroot}/picard-tools-1.123"
params.samtools_path = "${params.pkgroot}/samtools-1.0"
params.htseq_path = "${params.pkgroot}/htseq-count"
params.rnaqc_path = "${params.pkgroot}/RNASeQC_v1.1.8"
params.gatk_path = "${params.pkgroot}/GATK-3.3-0-g37228af"
params.starfusion_path = "${params.pkgroot}/STAR-Fusion-0.5.4"
params.pysam_path = "${params.pkgroot}/pysam-0.15.3"
params.starfusion_new_path = "${params.pkgroot}/STAR-Fusion-v1.6.0"

params.ptag = "current"
params.samplebase = "all"
params.results = "s3://k8s-temporal-efs/results/"
// params.results = "/efs-output/ngs_scratch/results/"
// params.results = "s3://nextflow-k8s-output-dev/results/"
params.location_prefix = "/efs/ngs_scratch/production/fastq/"

channelspec = params.location_prefix + "/" + params.runpath + "/" + "*fastq.gz"

// Resources
params.reference = "/efs/ngs_scratch/data/ngs_clia/newref/_REFs_002/"
//params.reference = "s3://k8s-temporal-efs/work/210113_NB501870_0629_AHW3VCBGXG/stage/2e/4101237d5b11dbe94bdcd3bda89fe9/_REFs_002/"
reference = file(params.reference)

params.resources = "/efs/ngs_scratch/data/ngs_clia/resources/"
resources = file(params.resources)

params.qcreference = "/efs/ngs_scratch/data/ngs_clia/newref/"
qcreference = file(params.qcreference)

params.fusiongenome = "/efs/ngs_scratch/data/ngs_clia/newref/Hg19_CTAT_resource_lib/"
fusiongenome = file(params.fusiongenome)

params.genome_ctat_lib = "/efs/ngs_scratch/data/ngs_clia/newref/GRCh37_gencode_v19_CTAT_lib/"
genome_ctat_lib = file(params.genome_ctat_lib)


// file refactor XXX
fastqs = Channel.fromPath( channelspec )
  .toSortedList()
  .flatMap()
  .map { path -> [((path.baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][2]) + '.' + ((path.baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][3])[-1..-1], path ] }
  .groupTuple(sort: true)
  .map { k,v -> ['sample_name': ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][1]),
                  'sample_id': ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][2]),
                  'lane_id': 'L' + ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][3])[-1..-1],
                  'read1_path': v[0],
                  'read2_path': v[1]] }

/*
// get all the fastq files
def getFastqs (c) {
Channel.fromPath( c )
  .toSortedList()
  .flatMap()
  .map { path -> [((path.baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][2]) + '.' + ((path.baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][3])[-1..-1], path ] }
  .groupTuple(sort:true)
  .map { k,v -> ['sample_name': ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][1]),
                  'sample_id': ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][2]),
                  'lane_id': 'L' + ((v[0].baseName =~ /([-\w]*)_(\w*)_(\w*)_(\w*)_(\w*)/)[0][3])[-1..-1],
                  'read1_path': v[0],
                  'read2_path': v[1]] }
        .view()
}

getFastqs (channelspec).set { fastqs }
*/

// completion handler

workflow.onComplete {
        println "Pipeline completed at: $workflow.complete"
        println "Command line: $workflow.commandLine"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        println "Duration: $workflow.duration"
}

workflow.onError {
        println "Error: $workflow.errorMessage"
        println "Full report: $workflow.errorReport"
}
/*
 *  Processes
 */

process printConfig {
        echo true
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'
        // pod = [nodeSelector: "instanceType=4g-2cpu"]

        output:
                val "printConfig" into pc

        """
                echo "product: ${product}"
                echo "configuration: "
                echo ${params.config}
                echo -n "Started pipeline at: "
                date
        """
}

process clearMemory {
        echo true
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'
        pod = [nodeSelector: "${config.instanceType.clearMemory}, uuid=${params.uuid}"]
        cpus "1"
        memory "2"
        
        input:
                file reference
                val c from pc

        output:
                val "cleared" into prerunMemoryStatus

        script:
                template 'clearMemory.sh'
}

process startup {
        tag { runid + ":" + params.samplebase }

        input:
                val mem from prerunMemoryStatus

        output:
                val "started" into kickoff

        script:
                template "startup.sh"
}

process alignment {
        echo true
        tag { sample_name + ":" + lane_id }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.alignment}, uuid=${params.uuid}"]
        cpus "${config.clustercores.alignment}" 
        memory "${config.clustermem.alignment}" 

    // limit here to avoid overloading the S3 download pool
        maxForks 16
        maxRetries 5
        //errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }
        errorStrategy 'retry'

        storeDir { "${params.results}/$runid/$sample_name/align/$lane_id" }
        //publishDir { "${params.results}/$runid/$sample_name/align/$lane_id" }

        input:
                file reference
                set sample_name, sample_id, lane_id, file(read1), file(read2) from fastqs
                //val k from kickoff.first()
                val k from kickoff

        output:
                set sample_name, file ("${lane_id}.Aligned.sortedByCoord.out.bam")  into alignments
                set sample_name, file ("${lane_id}.Chimeric.out.sam")  into fusion_sams
                set sample_name, file ("${lane_id}.Chimeric.out.junction")  into fusion_jxns

                set sample_name, file ("${lane_id}.Log.final.out")  into star_logs

                set sample_name, sample_id, lane_id  into meta

        script:
                template "alignment.sh"
}

process merge_new {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.merge}, uuid=${params.uuid}"]
        cpus "${config.clustercores.merge}"
        memory "${config.clustermem.merge}"
        //clusterOptions "--mem=32768"
        //time "30m"

        storeDir { "${params.results}/$runid/$sample_name/bam" }
        //publishDir { "${params.results}/$runid/$sample_name/bam" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file('bam') from alignments.groupTuple(sort: true, size: 4)

        // output:
	// 	set sample_name, file('merged.bam') into merged1, merged2, merged3
	// 	file("cleanup") into memoryStatus
	output:
		set sample_name, file('merged.bam'), file('merged.bai') into merged1, merged2, merged3, merged4, merged5
		file("cleanup") into memoryStatus
        script:
                template "merge.sh"
}

process clearMemoryOnExit {
        echo true
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'
        pod = [nodeSelector: "${config.instanceType.clearMemoryOnExit}, uuid=${params.uuid}"]
        cpus "1"
        memory "2"

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                file reference
                val status from memoryStatus

        output:
                val "cleared" into memoryCleanupStatus

        script:
                template 'clearMemory.sh'
}

process sortByName {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.sortByName}, uuid=${params.uuid}"]
        cpus "${config.clustercores.sortByName}" 
        memory "${config.clustermem.sortByName}" 

        storeDir { "${params.results}/$runid/$sample_name/exp" }
        //publisDir { "${params.results}/$runid/$sample_name/exp" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from merged1
                file resources

        output:
                set sample_name, file('sortedByName.bam')  into sorted1, sorted2

        script:
                template "sortByName.sh"
}

process counts {
        echo true
        // counts by gene_id
        // modified threshold from 20 down to 10 (to align with Hajime's pipeline
        // modified to take sorted bamfile
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.counts}, uuid=${params.uuid}"]
        cpus "${config.clustercores.counts}" 
        memory "${config.clustermem.counts}" 
        //clusterOptions "--mem=5120"
        //time "2h"

        storeDir { "${params.results}/$runid/$sample_name/exp" }
        //publisDir { "${params.results}/$runid/$sample_name/exp" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }


        input:
                set sample_name, file(bam) from sorted1
                file resources

        output:
                set sample_name, file('counts.txt')  into counts

        script:
                template "counts.sh"
}

process countsByGene {
        echo true
        // counts by gene
        // modified threshold from 20 down to 10 (to align with Hajime's pipeline
        // modified to take sorted bamfile
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.countsByGene}, uuid=${params.uuid}"]
        cpus "${config.clustercores.countsByGene}" 
        memory "${config.clustermem.countsByGene}"

        storeDir { "${params.results}/$runid/$sample_name/exp" }
        //publisDir { "${params.results}/$runid/$sample_name/exp" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from sorted2
                file resources

        output:
                set sample_name, file('countsByGene.txt')  into countsByGene

        script:merged2
                template "countsByGene.sh"
}

// Mark duplicates and index
process markDuplicates {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.markDuplicates}, uuid=${params.uuid}"]
        cpus "${config.clustercores.markDuplicates}"
        memory "${config.clustermem.markDuplicates}"

        storeDir { "${params.results}/$runid/$sample_name/bam" }
        //publisDir { "${params.results}/$runid/$sample_name/bam" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from merged2

        output:
                set sample_name, file("dedupped.bam"), file("dedupped.bai"), file("duplicates-metrics.txt")  into dedupped1, dedupped2, dedupped3, mark_duplicates_logs
                //set sample_name, file("duplicates-metrics.txt")  into mark_duplicates_logs

        script:
                template "markDuplicates.sh"
}

// after this point is AfirmaV2-specific

// RNA QC
process rnaqc {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.rnaqc}, uuid=${params.uuid}"]
        cpus "${config.clustercores.rnaqc}"
        memory "${config.clustermem.rnaqc}"

        storeDir { "${params.results}/$runid/$sample_name/rqc" }
        //publisDir { "${params.results}/$runid/$sample_name/rqc" }

        maxForks 6

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam), file(bai), file("duplicates-metrics.txt") from dedupped1

	output:
		set sample_name, file("index.html"), file ("report.html"), file ("countMetrics.html"), file ("refGene.txt"), file ("refGene.txt.idx"), file ("gc"), file (sample_name), file ("rRNA_intervals.list"), file ("meanCoverage*"), file ("gapLength*"), file ("genes.rpkm.gct")  into rnaqc_output
                set sample_name, file ("metrics.tsv") into qcmetrics1, qcmetrics2, qcmetrics3

        when:
                params.product == "afirmav2"

        script:
                template "rnaqc.sh"
}

//find exon skipping
// process PGA_MET_E14_Skipping {
// 	echo true
// 	tag {sample_name}

//         pod = [nodeSelector: "${config.instanceType.metE14Skipping}, uuid=${params.uuid}"]
//         cpus "${config.clustercores.metE14Skipping}"
//         memory "${config.clustermem.metE14Skipping}"
       
//         //publisDir { "${params.results}/$runid/$sample_name/pga_exonskipping" }
// 	storeDir { "${params.results}/$runid/$sample_name/pga_exonskipping" }

//         maxRetries 5
//         errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }
	
// 	input:
// 		set sample_name, file(bam) from merged4
// 		set sample_name2, file(metrics) from qcmetrics1
// 		file resources
		
// 	output:	
// 		set sample_name, file('MET_E14_Skipping.csv'), file('MET_E14_Skipping_JR_detailpd.csv') into met_e14_skipping
               	
// 	//when:
// 	//	params.product == "perceptaga"

// 	script:	
// 		// run MET_E14_Skipping : pyton file_name -i sampleName -n totalreads -f bamFilePath -a annotationFilePath eventName
//                 template "metE14Skipping.sh"
		
// }

// process PGA_EGFR_Var_III_Skipping {
// 	echo true
// 	tag {sample_name}

//         pod = [nodeSelector: "${config.instanceType.egfrViiiSkipping}, uuid=${params.uuid}"]
// 	cpus "${config.clustercores.egfrViiiSkipping}"
// 	memory "${config.clustermem.egfrViiiSkipping}"

//         //publisDir { "${params.results}/$runid/$sample_name/pga_exonskipping" }
// 	storeDir { "${params.results}/$runid/$sample_name/pga_exonskipping" }

//         maxRetries 5
//         errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }
	
// 	input:
// 		set sample_name, file(bam) from merged5
// 		set sample_name2, file(metrics) from qcmetrics2
// 		file resources
		
// 	output:	
// 		set sample_name, file('EGFR_Var_III.csv'), file('EGFR_Var_III_JR_detailpd.csv') into egfr_viii_skipping
	
// 	//when:
// 	//	params.product == "perceptaga"

// 	script:	
// 		// run EGFR_Var_III 
// 	        template "egfrVarIIIskipping.sh"
// }

// downsampling
process getReadInfo {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.getReadInfo}, uuid=${params.uuid}"]
        cpus "${config.clustercores.getReadInfo}" 
        memory "${config.clustermem.getReadInfo}"

        storeDir { "${params.results}/$runid/$sample_name/bam" }
        //publisDir { "${params.results}/$runid/$sample_name/bam" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from merged3
                set sample_name2, file("dedupped.bam"), file("dedupped.bai"), file ("duplicates-metrics.txt") from mark_duplicates_logs

        output:
                set sample_name, file(bam)  into readinfo_samplename1, readinfo_samplename2
//              stdout()  into (readinfo_downsample_percentage1, readinfo_downsample_percentage2)
                file ("downsample-output")  into (readinfo_downsample_percentage1, readinfo_downsample_percentage2)

        when:
                params.product == "afirmav2"


        script:
        """
        /usr/local/setup/nextflow/bin/vcalc.pl duplicates-metrics.txt > downsample-output
        sleep 60
        """
//      /usr/local/setup/nextflow/bin/vcalc.pl ${params.results}/$runid/$sample_name/bam/duplicates-metrics.txt

}

process downsample4m {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.downsample4m}, uuid=${params.uuid}"]
        cpus "${config.clustercores.downsample4m}"
        memory "${config.clustermem.downsample4m}"

        storeDir { "${params.results}/$runid/$sample_name/bam" }
        //publisDir { "${params.results}/$runid/$sample_name/bam" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from readinfo_samplename1
                //val ds_percent from readinfo_downsample_percentage1.map { it.tokenize(":")[0] }
                val ds_percent from readinfo_downsample_percentage1.map { it.text.trim().tokenize(":")[0] }

        output:
                set file("downsampled.4m.bam"), file ("dedupped.4m.bam"), file ("duplicates-metrics.4m.txt")  into downsample4m_output


        when:
                params.product == "afirmav2"

//      exec:
//              println "4m: sample name is ${sample_name}, tuple is ${ds_percent}"
        script:
                template "downsample4m.sh"

}

process downsample8m {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.downsample8m}, uuid=${params.uuid}"]
        cpus "${config.clustercores.downsample8m}"
        memory "${config.clustermem.downsample8m}"

        storeDir { "${params.results}/$runid/$sample_name/bam" }
        //publisDir { "${params.results}/$runid/$sample_name/bam" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam) from readinfo_samplename2
                //val ds_percent from readinfo_downsample_percentage2.map { it.tokenize(":")[1] }
                val ds_percent from readinfo_downsample_percentage2.map { it.text.trim().tokenize(":")[1] }

        output:
                set file("downsampled.8m.bam"), file ("dedupped.8m.bam"), file ("duplicates-metrics.8m.txt")  into downsample8m_output


        when:
                params.product == "afirmav2"

//      exec:
//              println "8m: sample name is ${sample_name}, tuple is ${ds_percent}"
        script:
                template "downsample8m.sh"
}

// samtools stats
// XXX added retry due to race condition? refactor  into multiple steps
process samtoolsStats {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.samtoolsStats}, uuid=${params.uuid}"]
        cpus "${config.clustercores.samtoolsStats}"
        memory "${config.clustermem.samtoolsStats}"

        storeDir { "${params.results}/$runid/$sample_name/stats" }
        //publisDir { "${params.results}/$runid/$sample_name/stats" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam), file(bai), file("duplicates-metrics.txt") from dedupped3
                file resources

        output:
                set sample_name, file("nomult.bam"), file ("nomult.bai"), file ("dedupped.stats.txt"), file ("nomult.stats.txt"), file ("dedupped.stattbl.txt"), file ("nomult.stattbl.txt"), file ("dedupped.stattbl2.txt"), file ("nomult.stattbl2.txt"), file ("dedupped.flags.txt"), file ("nomult.flags.txt"), file ("nomult.header.txt")  into samtool_stats_output


        when:
                params.product == "afirmav2"

        script:
                template "samtoolsStats.sh"
}

// // fusion, finally!
process concatJunctions {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.concatJunctions}, uuid=${params.uuid}"]
        cpus "${config.clustercores.concatJuctions}"
        memory "${config.clustermem.concatJuctions}"

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        storeDir { "${params.results}/$runid/$sample_name/fus" }
        //publisDir { "${params.results}/$runid/$sample_name/fus" }

        input:
                set sample_name, file('chimeric_junctions') from fusion_jxns.groupTuple();

	output:
		set sample_name, file("chimeric_junctions_merged") into junctions1, junctions2

        when:
                params.product == "afirmav2"


        script:
                template "concatJunctions.sh"
}

process mergeChimericSams {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.mergeChimericSams}, uuid=${params.uuid}"]
        cpus "${config.clustercores.mergeChimericSams}"
        memory "${config.clustermem.mergeChimericSams}"

        storeDir { "${params.results}/$runid/$sample_name/fus" }
        //publisDir { "${params.results}/$runid/$sample_name/fus" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }



        input:
                set sample_name, file('chimeric_sam') from fusion_sams.groupTuple()
        output:
                set sample_name, file("concat.chimeric.out.sam")   into fusions


        when:
                params.product == "afirmav2"

        script:
                template "mergeChimericSams.sh"
}

process starFusion {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.starFusion}, uuid=${params.uuid}"]
        cpus "${config.clustercores.starFusion}"
        memory "${config.clustermem.starFusion}"

        storeDir { "${params.results}/$runid/$sample_name/fus" }
        //publisDir { "${params.results}/$runid/$sample_name/fus" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(chimeric_sam) from fusions
                set sample_name, file(chimeric_junctions_merged) from junctions1

        output:
                set sample_name, file ("fusion_results")  into fusion_results


        when:
                params.product == "afirmav2"

        script:
                template "starFusion.sh"
}

//Additional Star Fusion 1.6.0 work for PGA samples
// process PGA_StarFusion {
// 	echo true
// 	tag { sample_name }

//         pod = [nodeSelector: "${config.instanceType.extendStarFusion}, uuid=${params.uuid}"]
// 	cpus "${config.clustercores.extendStarFusion}"
// 	memory "${config.clustermem.extendStarFusion}"

//         //publisDir { "${params.results}/$runid/$sample_name/fus" }
// 	storeDir { "${params.results}/$runid/$sample_name/fus" }

//         maxRetries 5
//         errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

// 	input:
// 		set sample_name, file(metrics) from qcmetrics3
// 		set sample_name, file(chimeric_junctions_merged) from junctions2
              
// 	output:
// 		set sample_name, file("pga_fusion_results") into pga_fusion_results
// 		//file("star-fusion.fusion_predictions.abridged.tsv"), file("star-fusion.fusion_predictions.tsv"), file("star-fusion.fusion_predictions.abridged.coding_effect.tsv") into extendstarfusion_output

// 	//when:
// 	//	params.product == "perceptaga"


// 	script:
// 		template "starFusion_1.6.0.sh"
// }


// beyond this point corresponds to Hajime's GATK pipeline
process split {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.split}, uuid=${params.uuid}"]
        cpus "${config.clustercores.split}"
        memory "${config.clustermem.split}"

        storeDir { "${params.results}/$runid/$sample_name/split" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        ////publisDir { "${params.results}/$sample_name/split" }

        input:
                set sample_name, file(bam), file(bai), file("duplicates-metrics.txt") from dedupped2
                file resources

        output:
                set sample_name, file("split.bam"), file("split.bai")  into split

        when:
                params.product == "afirmav2"


        script:
                template "split.sh"
}


// /*
//  * 4. GATK: Indel Realignment (optional but recommended--using a list of known indels will both speed up processing and improve accuracy, although it is not required
//  * Create a target list of intervals to be realigned (RealignerTargetCreator step would need to be done just once for a single set of indels)):
//  * (this argument recommended to speed up the process *if* this is only a temporary file; otherwise, use the default value)
//  */
process indelRealignerTarget {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.indelRealignerTarget}, uuid=${params.uuid}"]
        cpus "${config.clustercores.indelRealignerTarget}"
        memory "${config.clustermem.indelRealignerTarget}"

        storeDir { "${params.results}/$runid/$sample_name/realign" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }
        //publishDir { "${params.results}/$sample_name/realign" }

        input:
                set sample_name, file(bam), file(bai) from split
                file resources

        output:
                set sample_name, file(bam), file(bai), file ("target-intervals.list")  into realignTarget

        when:
                params.product == "afirmav2"


        script:
                template "indelRealignerTarget.sh"
}

process indelRealignment {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.indelRealignment}, uuid=${params.uuid}"]
        cpus "${config.clustercores.indelRealignment}"
        memory "${config.clustermem.indelRealignment}"

        storeDir { "${params.results}/$runid/$sample_name/realign" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/realign" }

        input:
                set sample_name, file(bam), file(bai), file ("target-intervals.list") from realignTarget
                file resources

        output:
                set sample_name, file("realigned.bam"), file("realigned.bai")  into realigned

        when:
                params.product == "afirmav2"


        script:
                template "indelRealignment.sh"
}

process recalPre {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.recalPre}, uuid=${params.uuid}"]
        cpus "${config.clustercores.recalPre}"
        memory "${config.clustermem.recalPre}"

        storeDir { "${params.results}/$runid/$sample_name/recalibration" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/recalibration" }

        input:
                set sample_name, file(bam), file(bai) from realigned
                file resources

        output:
                set sample_name, file(bam), file(bai),  file("pre-recal-tbl.txt")  into pre_recal_table1, pre_recal_table2, pre_recal_table3

        when:
                params.product == "afirmav2"


        script:
                template "recalPre.sh"
}

process recalibration {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.recalibration}, uuid=${params.uuid}"]
        cpus "${config.clustercores.recalibration}"
        memory "${config.clustermem.recalibration}"

        storeDir { "${params.results}/$runid/$sample_name/recalibration" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/recalibration" }

        input:
                set sample_name, file(bam), file(bai), file("pre-recal-tbl.txt") from pre_recal_table1
                file resources

        output:
                set sample_name, file("calibrated.bam"), file("calibrated.bai")  into calibrated1, calibrated2, calibrated3, calibrated4, calibrated5

        when:
                params.product == "afirmav2"


        script:
                template "recalibration.sh"
}

process recalPost {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.recalPost}, uuid=${params.uuid}"]
        cpus "${config.clustercores.recalPost}"
        memory "${config.clustermem.recalPost}"

        storeDir { "${params.results}/$runid/$sample_name/recalibration" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/recalibration" }

        input:
                set sample_name, file(bam), file(bai), file("pre-recal-tbl.txt") from pre_recal_table2
                file resources

        output:
                set sample_name, file ("post-recal-tbl.txt")  into post_recal_table

        when:
                params.product == "afirmav2"


        script:
                template "recalPost.sh"

}

process recalAnalyze {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.recalAnalyze}, uuid=${params.uuid}"]
        cpus "${config.clustercores.recalAnalyze}"
        memory "${config.clustermem.recalAnalyze}"

        storeDir { "${params.results}/$runid/$sample_name/recalibration" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/recalibration" }

        input:
                set sample_name, file(bam), file(bai), file("pre-recal-tbl.txt") from pre_recal_table3
                set sample_name2, file ("post-recal-tbl.txt") from post_recal_table
                file resources

        output:
                file("recalplots.csv")  into recal_covariates

        when:
                params.product == "afirmav2"


        script:
                template "recalAnalyze.sh"

}

process variantSamtoolsDepths {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantSamtoolsDepths}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantSamtoolsDepths}"
        memory "${config.clustermem.variantSamtoolsDepths}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, file(bam), file(bai) from calibrated1
                file resources

        output:
                set sample_name, file ("ThyroVIPv2.depths.txt"), file ("ThyroVIPv2.nomult.depths.txt")  into variant_depths

        when:
                params.product == "afirmav2"


        script:
                template "variantSamtoolsDepths.sh"
}

process variantGenotyping {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantGenotyping}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantGenotyping}"
        memory "${config.clustermem.variantGenotyping}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file(bam), file(bai) from calibrated2
                file resources

        output:
                set sample_name, file("ThyroVIPv2.genotypes.raw.vcf")  into raw_thyro_variants

        when:
                params.product == "afirmav2"


        script:
                template "variantGenotyping.sh"
}

process variantDiscovery {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantDiscovery}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantDiscovery}"
        memory "${config.clustermem.variantDiscovery}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file(bam), file(bai) from calibrated3
                file resources
                file reference

        output:
                set sample_name, file("IllmRNAAccessCap.discover.raw.vcf"), file(bam), file(bai)  into raw_discovery_variants

        when:
                params.product == "afirmav2"


        script:
                template "variantDiscovery.sh"
}


process mitoVariantDiscovery {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.mitoVariantDiscovery}, uuid=${params.uuid}"]
        cpus "${config.clustercores.mitoVariantDiscovery}"
        memory "${config.clustermem.mitoVariantDiscovery}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file(bam), file(bai) from calibrated4
                file resources
                file reference

        output:
                set sample_name, file("MTgenes.mt.raw.vcf")  into raw_mito_variants

        when:
                params.product == "afirmav2"


        script:
                template "mitoVariantDiscovery.sh"
}

process datVariantDiscovery {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.datVariantDiscovery}, uuid=${params.uuid}"]
        cpus "${config.clustercores.datVariantDiscovery}"
        memory "${config.clustermem.datVariantDiscovery}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file(bam), file(bai) from calibrated5
                file resources
                file reference

        output:
                set sample_name, file("av2dat.genotypes.raw.vcf")  into raw_dat_variants

        when:
                params.product == "afirmav2"


        script:
                template "datVariantDiscovery.sh"
}

process variantFilterGenotypingThyroVIPv2 {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantFilterGenotypingThyroVIPv2}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantFilterGenotypingThyroVIPv2}"
        memory "${config.clustermem.variantFilterGenotypingThyroVIPv2}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file("ThyroVIPv2.genotypes.raw.vcf") from raw_thyro_variants
                file resources

        output:
                set sample_name, file("ThyroVIPv2.genotypes.filter.vcf")  into filtered_genotype_variants

        when:
                params.product == "afirmav2"


        script:
                template "variantFilterGenotypingThyroVIPv2.sh"
}

process variantFilterDiscoveryRNAAccessCap {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantFilterDiscoveryRNAAccessCap}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantFilterDiscoveryRNAAccessCap}"
        memory "${config.clustermem.variantFilterDiscoveryRNAAccessCap}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:

                set sample_name, file("IllmRNAAccessCap.discover.raw.vcf"), file(bam), file(bai) from raw_discovery_variants
                file resources

        output:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.vcf"), file(bam), file(bai)  into filtered_discovery_variants1, filtered_discovery_variants2, filtered_discovery_variants3

        when:
                params.product == "afirmav2"


        script:
                template "variantFilterDiscoveryRNAAccessCap.sh"
}

process mitoVariantFiltering {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.mitoVariantFiltering}, uuid=${params.uuid}"]
        cpus "${config.clustercores.mitoVariantFiltering}"
        memory "${config.clustermem.mitoVariantFiltering}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:

                set sample_name, file("MTgenes.mt.raw.vcf") from raw_mito_variants
                file resources

        output:
                set sample_name, file ("MTgenes.mt.filter.vcf")  into filtered_mito_variants1, filtered_mito_variants2

        when:
                params.product == "afirmav2"


        script:
                template "mitoVariantFiltering.sh"
}

process datVariantFiltering {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.datVariantFiltering}, uuid=${params.uuid}"]
        cpus "${config.clustercores.datVariantFiltering}"
        memory "${config.clustermem.datVariantFiltering}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:

                set sample_name, file("av2dat.genotypes.raw.vcf") from raw_dat_variants
                file resources

        output:
                set sample_name, file ("av2dat.genotypes.filter.vcf")  into filtered_dat_variants

        when:
                params.product == "afirmav2"


        script:
                template "datVariantFiltering.sh"
}

process variantEvaluationDiscovery {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantEvaluationDiscovery}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantEvaluationDiscovery}"
        memory "${config.clustermem.variantEvaluationDiscovery}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.vcf"), file(bam), file(bai) from filtered_discovery_variants1

        output:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.eval.txt")  into filtered_discovery_variants_eval

        when:
                params.product == "afirmav2"


        script:
                template "variantEvaluationDiscovery.sh"
}

process variantFilterDiscovery1000gHiconf {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantFilterDiscovery1000gHiconf}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantFilterDiscovery1000gHiconf}"
        memory "${config.clustermem.variantFilterDiscovery1000gHiconf}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.vcf"), file(bam), file(bai) from filtered_discovery_variants2
        output:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.1000GHiConf.mask.vcf")  into filtered_discovery_variants_1000g_hiconf

        when:
                params.product == "afirmav2"


        script:
                template "variantFilterDiscovery1000gHiconf.sh"
}

process variantFilterDiscoveryReadBackedPhasing {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.variantFilterDiscoveryReadBackedPhasing}, uuid=${params.uuid}"]
        cpus "${config.clustercores.variantFilterDiscoveryReadBackedPhasing}"
        memory "${config.clustermem.variantFilterDiscoveryReadBackedPhasing}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.vcf"), file(bam), file(bai) from filtered_discovery_variants3
        output:
                set sample_name, file ("IllmRNAAccessCap.discover.filter.phased.vcf")  into filtered_discovery_variants_phased

        when:
                params.product == "afirmav2"


        script:
                template "variantFilterDiscoveryReadBackedPhasing.sh"
}

process mitoVariantEvaluation {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.mitoVariantEvaluation}, uuid=${params.uuid}"]
        cpus "${config.clustercores.mitoVariantEvaluation}"
        memory "${config.clustermem.mitoVariantEvaluation}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file ("MTgenes.mt.filter.vcf") from filtered_mito_variants1

        output:
                set sample_name, file ("MTgenes.mt.filter.eval.txt")  into filtered_mito_variants_eval

        when:
                params.product == "afirmav2"


        script:
                template "mitoVariantEvaluation.sh"
}

process mitoVariantMasking {
        echo true
        tag { sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.mitoVariantMasking}, uuid=${params.uuid}"]
        cpus "${config.clustercores.mitoVariantMasking}"
        memory "${config.clustermem.mitoVariantMasking}"

        storeDir { "${params.results}/$runid/$sample_name/vcf" }

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        //publishDir { "${params.results}/$sample_name/vcf" }

        input:
                set sample_name, file ("MTgenes.mt.filter.vcf") from filtered_mito_variants2
        output:
                set sample_name, file ("MTgenes.mt.filter.1000GHiConf.mask.vcf")  into filtered_mito_variants_1000g_hiconf

        when:
                params.product == "afirmav2"


        script:
                template "mitoVariantMasking.sh"
}



process uploadResults {
        echo true
        tag { runid:sample_name }
        container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'

        pod = [nodeSelector: "${config.instanceType.uploadResults}, uuid=${params.uuid}"]
        memory "2g"

        maxRetries 5
        errorStrategy { sleep(Math.pow(2, task.attempt) * 15 as long); return 'retry' }

        input:
                set sample_name, sample_id, lane_id from meta
                set sample_id2, final_output from filtered_discovery_variants_1000g_hiconf
                set sample_id3, final_output2 from filtered_discovery_variants_eval

        output:
                set sample_id, sample_name  into cleanup

        when:
                params.product == "afirmav2"


        script:
                template "uploadResults.sh"
}


// process cleanupDirFiles {
//         echo true
//         tag { runid:sample_name }
//     container '183449805178.dkr.ecr.us-east-1.amazonaws.com/nextflow-worker:latest'
//     pod = [nodeSelector: "${config.instanceType.clearMemory}"]
//         input:
//                 set sample_id, sample_name from cleanup

//         when:
//                 params.product == "afirmav2" && sample_name != null


//         script:
//                 template "cleanupFiles.sh"
// }