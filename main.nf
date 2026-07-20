include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_MERGE } from './modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_INDEX } from './modules/nf-core/bcftools/index/main'

process BCFTOOLS_FILLTAGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.gz"), emit: vcf
    tuple val(meta), path("*.tbi"), emit: index

    script:
    """
    bcftools +fill-tags \
      -o ${vcf.baseName}.filltags.vcf.gz ${vcf} -- \
      -t FORMAT/VAF

    bcftools index -t ${vcf.baseName}.filltags.vcf.gz
    """
}


workflow {
    main:
    // Load bam alignments
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], file(row.bam))
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }

    // Load reference
    reference_ch = channel.fromPath(params.reference, checkIfExists:true)
    | map { n -> tuple(n.baseName, n) }

    SAMTOOLS_FAIDX(
      reference_ch.map { meta, fasta -> tuple(meta, fasta, [])},
      []
    )

    indexed_reference_ch = reference_ch
    | join(SAMTOOLS_FAIDX.out.fai)

    input_ch = bam_ch
    | combine(indexed_reference_ch)
    | map { meta, bam, meta2, fasta, fai -> tuple(meta, bam, meta2, fasta, fai, [])}

    BCFTOOLS_MPILEUP(
      input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta, bam, [], []) },
      input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta2, fasta, fai) },
      input_ch.map { meta, bam, meta2, fasta, fai, other -> other }
    )

    BCFTOOLS_MPILEUP.out.vcf | BCFTOOLS_FILLTAGS

    norm_input_ch = BCFTOOLS_FILLTAGS.out.vcf | join (BCFTOOLS_FILLTAGS.out.index) | combine(reference_ch)

    BCFTOOLS_NORM(
      norm_input_ch.map{ meta, vcf, tbi, meta2, fasta -> tuple(meta, vcf, tbi)},
      norm_input_ch.map{ meta, vcf, tbi, meta2, fasta -> tuple(meta2, fasta)}
    )

    BCFTOOLS_INDEX (
      BCFTOOLS_NORM.out.vcf
    )

    merged_vcf_ch = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_INDEX.out.index)
    // .map { records ->
    //     def meta = [id: 'merged']
    //     def vcfs = records.collect { it[1] }
    //     def tbis = records.collect { it[2] }

    //     tuple(meta, vcfs, tbis)
    // }
    // merged_vcf_ch | view
}
