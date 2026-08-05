include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM } from './modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_MERGE } from './modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_INDEX } from './modules/nf-core/bcftools/index/main'
include { BCFTOOLS_FILTER } from './modules/nf-core/bcftools/filter/main'

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

process BCFTOOLS_QUERY {

    input:
    path(vcf)

    output:
    path "sample_variants_results.tsv"

    script:
    """
    bcftools query \
      -HH \
      -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%Gene\\t%Product\\t%VariantType\\t%FeatureType\\t%IsPseudo\\t%IsGenic\\t%IsTransition\\t%IsSynonymous\\t%AminoAcidChange\\t%SNPCodonPosition\\t%CodonPosition\\t%AltAminoAcid\\t%RefAminoAcid\\t%AltCodon\\t%RefCodon[\\t%VAF]' \
      ${vcf} \
      | sed 's/:VAF//g' \
      | sed 's/^#CHROM/CHROM/g' \
      > sample_variants_results.tsv
    """
}

process VCF_ANNOTATOR {
    tag "annotate"

    input:
    path(vcf)
    path(genbank)

    output:
    path "annotated.vcf"

    script:
    """
    vcf-annotator \
        ${vcf} \
        ${genbank} \
        --output annotated.vcf
    """
}


workflow {
    main:
    // ============ Load Files ==================== //
    // BAM: Should be run through bam-filtered already
    // REF: Reference to call variants against
    // ============================================ //
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

    reference_ch = channel.fromPath(params.reference, checkIfExists:true)
    | map { n -> tuple(n.baseName, n) }

    SAMTOOLS_FAIDX(
      reference_ch.map { meta, fasta -> tuple(meta, fasta, [])},
      []
    )

    indexed_reference_ch = reference_ch
    | join(SAMTOOLS_FAIDX.out.fai)

    // === Optional debug
    // bam_ch.view { "bam_ch: $it" }
    // reference_ch { "reference_ch: $it" }
    // indexed_reference_ch { "indexed_reference_ch: $it" }

    // ============ MPILEUP and CALL ============== //
    // BAM = sorted bam against a reference
    // Compute genotype likelihoods and base counts at each genomic position
    // Internally the nf-core mmpileup module also runs a bcftools call, so the output is already vcf files
    // ============================================ //
    mpileup_input_ch = bam_ch
    | combine(indexed_reference_ch)
    | map { meta, bam, meta2, fasta, fai -> tuple(meta, bam, meta2, fasta, fai, [])}

    BCFTOOLS_MPILEUP(
      mpileup_input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta, bam, [], []) }, // BAM
      mpileup_input_ch.map { meta, bam, meta2, fasta, fai, other -> tuple(meta2, fasta, fai) }, // REF
      mpileup_input_ch.map { meta, bam, meta2, fasta, fai, other -> other } // SAVE_MPILEUP or FALSE here since bcftools call is in this module
    )

    mpileup_ch = BCFTOOLS_MPILEUP.out.vcf
    // === Optional debug
    // mpileup_ch.view { "mpileup_ch: $it" }

    // ============ FILLTAGS ====================== //
    // VCF = Variants
    // FIll in a FORMAT/VAF tag
    // ============================================ //
    BCFTOOLS_FILLTAGS(
      mpileup_ch
    )

    // ============ FILTER ====================== //
    // VCF = Variants
    // Filter variants by depth (So if you have a 1,000 minimum depth, 20 DP is cutoff for a 2% variant)
    // ============================================ //
    filter_input_ch = BCFTOOLS_FILLTAGS.out.vcf.join(BCFTOOLS_FILLTAGS.out.index)

    BCFTOOLS_FILTER(
      filter_input_ch
    )

    // ============ NORM ====================== //
    // VCF = Filtered variants
    // Standardizes and left-align indels, ensure match against ref genome, split multiallelic sites to separate lines, removes duplicate records
    // ============================================ //

    norm_input_ch = BCFTOOLS_FILTER.out.vcf | join(BCFTOOLS_FILTER.out.index) | combine(reference_ch)
    BCFTOOLS_NORM(
      norm_input_ch.map{ meta, vcf, tbi, meta2, fasta -> tuple(meta, vcf, tbi)},
      norm_input_ch.map{ meta, vcf, tbi, meta2, fasta -> tuple(meta2, fasta)}
    )

    // ============ MERGE ====================== //
    // VCF = Filtered and normalized variants
    // Standardizes and left-align indels, ensure match against ref genome, split multiallelic sites to separate lines, removes duplicate records
    // ============================================ //
    merge_input_ch = BCFTOOLS_NORM.out.vcf
    | join(BCFTOOLS_NORM.out.index)
    | map { meta, vcf, tbi -> tuple(meta, vcf, tbi)}
    //| view {"a: $it"}
    | collect(flat:false)
    //| view {"b: $it"}
    | map { 
        records ->
        def meta = [id: "merged"]

        def vcfs = records.collect{ it[1] }
        def tbis = records.collect{ it[2] }
        tuple(meta, vcfs, tbis)
      }
    | combine(indexed_reference_ch)

    //merge_input_ch.view { "merge_input_ch: $it" }

    // tuple val(meta), path(vcfs), path(tbis), path(bed)
    // tuple val(meta2), path(fasta), path(fai)
    BCFTOOLS_MERGE(
      merge_input_ch.map{meta, vcfs, tbis, meta2, fasta, fai -> tuple(meta, vcfs, tbis, [])},
      merge_input_ch.map{meta, vcfs, tbis, meta2, fasta, fai -> tuple(meta2, fasta, fai, [])}
    )

}
