#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
                        	Pre-process NGS reads
========================================================================================

preprocess_ngs_reads
Pipeline to trim, filter, and quality check Illumina (paired-end) or Nanopore reads.
#### Homepage / Documentation
https://github.com/mattlhowes/nextflow-pipelines/preprocess_ngs_reads
#### Author: Matt Howes
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
	Usage:
	Pipeline to trim, filter, and quality check Illumina paired-end fastq or Nanopore read fast5/fastq data.
	
	For illumina pe fastq read data: pipeline conducts Fastp and FastQC
	For nanopore fast5 or fastq read data: pipeline conducts Porechop and Nanoplot

	The command to run the pipeline is as follows:

	Illumina paired end reads (2 fastq files per sample):
	nextflow run preprocess_fastq --input "path/*_reads_{1,2}.fastq.gz" --sequencer_type illumina
    
	Nanopore fastq read data (1 fastq file per sample):
	nextflow run preprocess_fastq --input "path/*.fastq" --sequencer_type nanopore

	Mandatory arguments:
	--input		Path to input read data (must be surrounded with quotes if specifying pairs of illumina fastq files)
	--sequencer_type		'illumina' or 'nanopore' sequencer type which produced the reads
	-profile                    	Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
	Options:
	--outdir		The output directory where the results will be saved, default is "results".
	--skip_fastqc		Illumina only. Skips FastQC step.
	--skip_nanoplot		Nanopore only. Skips Nanoplot step.
      
	"""
	.stripIndent()
}


/*
 * Configuration 
 */

// Help
params.help = false
if (params.help) {
    helpMessage()
    exit 0
}


// General options
params.outdir = "results"

// Illumina options
params.skip_fastqc = false

// Nanopore options
params.skip_nanoplot = false


// Log
log.info """
	Read input: ${params.input}
	Output directory: ${params.outdir}
      """
      .stripIndent()


workflow {
   	if (params.sequencer_type == 'illumina'){
		read_pairs_ch = channel.fromFilePairs( params.input, checkIfExists: true )
		FASTP(read_pairs_ch)
		if (params.skip_fastqc) {   	
    			FASTQC(read_pairs_ch)
		}
	}
	else if (params.sequencer_type == 'nanopore'){
		nanopore_fq_read = channel.fromPath( params.input, checkIfExists: true )
		TRIM_ADAPTORS(nanopore_fq_read) 
		if (params.skip_fastqc) {   	  	
    			NANOPLOT(TRIM_ADAPTORS.out)
		}
	}
	else
		error "Error: --sequencer_type required: 'illumina' or 'nanopore'."
}
 

/*
 * Fastp trimming and filtering FOR Illumina reads
 */

process FASTP {
    	tag "FASTP on $sample_id"
	publishDir "${params.outdir}/${sample_id}/fastp/", mode: 'copy'
 
	when: (params.assembly_type == 'short')

    	input:
    	tuple val(sample_id), path(input)
 
    	output:
    	path("*")
	    

    	script:
    	"""
    	fastp --in1 "${input[0]}" --in2 "${input[1]}" \\
	    --out1 "clean_${input[0]}" --out2 "clean_${input[1]}" \\
          --cut_front \\
          --cut_front_mean_quality 20 \\
          --cut_front_window_size 4 \\
          --cut_right \\
          --trim_front1 5 \\
          --trim_front2 5 \\
          --detect_adapter_for_pe \\
          -l 75 \\
          -h "${sample_id}.html" \\
          -j "${sample_id}.json"
    	"""
}


/*
 * FastQC quality report for Illumina reads
 */

process FASTQC {
	tag "FASTQC on $sample_id"
    	publishDir "${params.outdir}/${sample_id}/fastqc/", mode: 'copy'
 
    	input:
    	tuple val(sample_id), path(input)
 
    	output:
    	path("*.{html,zip}")

    	script:
    	"""
    	fastqc "${input[0]}" "${input[1]}"
    	"""
}


/*
 * Adapter trimming for nanopore reads
 */

process TRIM_ADAPTORS {
	tag "TRIM_ADAPTORS on $sample_id"
    	publishDir "${params.outdir}/${sample_id}/trimming/", mode: 'copy'

    	input:
    	path(read)
	
    	output:
	path("*")

    	script:
	sample_id = "${read.baseName}"
    	"""
    	porechop -i "${read}" -t "${task.cpus}" -o trimmed_"${read}"
    	"""
}


/*
 * Quality check for nanopore reads and plots of quality and length
 */

process NANOPLOT {
    	tag "NANOPLOT on $sample_id"
	publishDir "${params.outdir}/${sample_id}/nanoplot_QC/", mode: 'copy'
	cpus = "$cpus"

    	input:
    	path(trimmed_read)

    	output:
    	path("*")

    	script:
	sample_id = "${trimmed_read.baseName}".replaceAll(/trimmed_/, "")	
    	"""
    	NanoPlot -t "${task.cpus}" --title "${sample_id}" -c darkblue --fastq "${trimmed_read}"
    	"""
}
