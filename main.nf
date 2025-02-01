nextflow.enable.dsl = 2

params.input_files = "$projectDir/data/*_{R1,R2}_001.fastq.gz"
params.output_dir = "results/data"
params.database = "kaijuDb"

process Classfier {
    tag "virus-specific classification"
    publishDir "${params.output_dir}", pattern: "*.kaiju", mode: "copy"
    container "harby/taxa_tool:latest"
        
    input:
    tuple val(sample_id), path(reads)
    each path(database)

    output:
    tuple val(sample_id) path ("${sample_id}_viruses.kaiju") emit: class 

    script:
    """
    mkdir kaiju_class
    kaiju -t "${database}/nodes.dmp" \
      -f "${database}/kaiju_db_viruses.fmi" \
      -i ${reads[0]} \
      -j ${reads[1]} \
      -o "${sample_id}_viruses.kaiju"
    """
}
process Kronaformat {
    tag "Convert Kaiju output to Krona-compatible format"
    container "harby/taxa_tool:latest"


    input:
    tuple val(sample_id) path (sample_id_kaiju_file)
    each path(database)
    
    
    output:
     tuple val(sample_id) path ("${sample_id}_viruses.krona") emit: format 


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
    container "harby/taxa_tool:latest"
    publishDir ${params.output_dir},  mode: "copy"

    input:
    tuple val(sample_id) path(sample_id_krona_file)
    output:
    path html_files

    script:
    """
    mkdir html_files
    perl /opt/KronaTools-2.8.1/scripts/ImportText.pl \
     "$sample_id_krona_file" \
     -o "html_files/${sample_id}_viruses.html"
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
    Kronaformat(Classfier.out.class, database_ch)
    Visualize(Kronaformat.out.format)


}