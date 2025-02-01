process CheckQC {
    tag "Checking the quality of raw reads"
    container ""
    publishDir "$params.output_dir", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_result"

    script:
    """
    mkdir fastqc_${sample_id}_result
    fastqc -o fastqc_${sample_id}_result -t ${reads[0]} ${reads[1]}    
    """


}
process Trimming {
    tag "Trimming the poor quality reads"
    container ""

    input:
	tuple val(sample_id), path(reads)

	output:
	tuple val(sample_id), path("${sample_id}_R1_001.trimmed.fastq.gz"), path("${sample_id}_R2_001.trimmed.fastq.gz"), \
		emit: trim_reads
	path "${sample_id}_fastp.html", emit: html
	path "${sample_id}_fastp.json", emit: json

	script:
    """
	fastp -w ${task.cpus}\\
	 	-f ${params.trim_front_read_01}\\
		-t ${params.trim_tail_read_01}\\
		-F ${params.trim_front_read_02}\\
		-T ${params.trim_tail_read_02}\\
	  	-i ${reads[0]}\\
	   	-I ${reads[1]}\\
	    -o ${sample_id}_R1_001.trimmed.fastq.gz\\
		-O ${sample_id}_R2_001.trimmed.fastq.gz\\
		-h ${sample_id}_fastp.html\\
		-j ${sample_id}_fastp.json \\
	"""
}