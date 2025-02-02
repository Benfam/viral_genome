nextflow.enable.dsl = 2

params.input_files = "$projectDir/data/*_{R1,R2}_*.fastq.gz"
params.output_dir = "results/data"
params.database = "kaijuDb"

process Classfier {
    tag "virus-specific classification"
    publishDir "${params.output_dir}", pattern: "*.kaiju", mode: "copy"
    container "harbby1/taxatools"
    cpus 6
        
    input:
    tuple val(sample_id), path(reads)
    each path(database)

    output:
    tuple val(sample_id), path ("${sample_id}_viruses.kaiju"), emit: taxa

    script:
    """
    kaiju -z $task.cpus -t "${database}/nodes.dmp" \
      -f "${database}/kaiju_db_viruses.fmi" \
      -i ${reads[0]} \
      -j ${reads[1]} \
      -o "${sample_id}_viruses.kaiju"
    """
}
process Kronaformat {
    tag "Convert Kaiju output to Krona-compatible format"
    container "harbby1/taxatools"


    input:
    tuple val(sample_id), path (sample_id_kaiju_file)
    each path(database)
    
    
    output:
     tuple val(sample_id), path ("${sample_id}_viruses.krona"), emit: format 


    script:
    """
    kaiju2krona -t "$database/nodes.dmp" \
            -n "$database/names.dmp" \
            -i "$sample_id_kaiju_file" \
            -o "${sample_id}_viruses.krona"
    """
}
process Visualize {
    tag "Generating Krona HTML visualization..."
    container "harbby1/taxatools"
    publishDir "${params.output_dir}/html", pattern: "*.html", mode: "copy"

    input:
    tuple val(sample_id), path(sample_id_krona_file)
    output:
    path ("*.html", arity: '1..*')

    script:
    """
    perl /opt/KronaTools-2.8.1/scripts/ImportText.pl \
     "$sample_id_krona_file" \
     -o "${sample_id}_viruses.html"
    """
}

workflow {
     log.info """\
         Viral Metagenomics - N F   P I P E L I N E
         ==============================================
         Input files    : ${params.input_files}
        
         """
         .stripIndent()

    input_files_ch = Channel.fromFilePairs(params.input_files, checkIfExists: true)
    input_files_ch.view()
    database_ch = Channel.fromPath(params.database, checkIfExists: true, type: 'dir')
    database_ch.view()
    Classfier(input_files_ch, database_ch)
    Kronaformat(Classfier.out.taxa, database_ch)
    Visualize(Kronaformat.out.format)


}
