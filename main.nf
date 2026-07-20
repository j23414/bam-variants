include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'

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
    | view

    BCFTOOLS_MPILEUP(
      input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta, bam, [], []) },
      input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta2, fasta, fai) },
      input_ch.map { meta, bam, meta2, fasta, fai, other -> other }
    )

}
