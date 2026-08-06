process VCF_ANNOTATOR {
    tag "annotate"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity','apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/vcf-annotator:1.0.0--1d57e1a9ec649578'
        : 'community.wave.seqera.io/library/vcf-annotator:1.0.0--742c39e931a8e8ce'}"

    input:
    tuple val(meta), path(vcf)
    path(genbank)

    output:
    path "annotated.vcf"

    script:
    """
    vcf-annotator \\
        ${vcf} \\
        ${genbank} \\
        --output annotated.vcf
    """
}