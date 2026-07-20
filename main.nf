include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'

workflow {
    main:
    // Load bam alignments
    if (params.samplesheet) {
      bam_ch = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple([id: row.sample], [file(row.bam)])
          }
    } else if (params.bam) {
        bam_ch = channel.fromPath(params.bam, checkIfExists: true)
          | map { bamfile -> tuple([id: bamfile.baseName], bamfile) }
    } else {
        error "Please specify either --samplesheet samplesheet.csv or --bam 'data/*.bam'"
    }
    bam_ch | view

    // Load reference
    reference_ch = channel.fromPath(params.reference, checkIfExists:true)
    | map { n -> tuple(n.baseName, n) }

    SAMTOOLS_FAIDX(
      reference_ch.map { meta, fasta -> tuple(meta, fasta, [])}
    )

    indexed_reference_ch = reference_ch
    | join(SAMTOOLS_FAIDX.out.fai)

    BCFTOOLS_MPILEUP(
      bam_ch.map { meta, bam -> tuple(meta, bam, [], []) },
      indexed_reference_ch,
      false
    )

}
