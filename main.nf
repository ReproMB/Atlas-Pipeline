#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPROMB ATLAS 16S PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Uniform reprocessing of public 16S rRNA amplicon studies, one processing
    unit per study x primer pair x body-site group.

    Steps:
      1. Parse the study table, clean primer strings, deduplicate primer pairs
      2. Select samples from the per-study metadata and write QIIME 2 manifests
      3. Import FASTQ into QIIME 2 (paired-end or single-end)
      4. Optional primer removal by sequence (cutadapt)
      5. Quality profile and per-study truncation parameters
      6. Denoising: DADA2 (default) or Deblur, with QC and one adaptive retry (DADA2)
      7. One naive Bayes classifier per primer pair (SILVA 138.2)
      8. Keep d__Bacteria features and collapse to genus
      9. Export tables; per-study diversity and Lactobacillus summaries

    Usage:
      nextflow run main.nf \
          --input_csv studies.csv \
          --data_dir /path/to/fetchngs/data \
          --silva_seqs silva-138.2-ssu-nr99-seqs.qza \
          --silva_tax silva-138.2-ssu-nr99-tax.qza \
          --outdir results \
          -profile docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

// ── Helpers ─────────────────────────────────────────────────────────────────
// IUPAC-aware reverse complement. Needed so cutadapt can strip 3' read-through
// (a read running into the reverse complement of the opposite primer).
def IUPAC_COMPLEMENT = [
    'A':'T', 'C':'G', 'G':'C', 'T':'A', 'U':'A',
    'R':'Y', 'Y':'R', 'S':'S', 'W':'W', 'K':'M', 'M':'K',
    'B':'V', 'V':'B', 'D':'H', 'H':'D', 'N':'N',
]
def revcomp = { String primer ->
    if (!primer) return ''
    primer.toUpperCase().reverse().collect { c ->
        IUPAC_COMPLEMENT.getOrDefault(c as String, 'N')
    }.join('')
}

// Parse the per-study manual denoise params sheet (from prompt_denoise_params.py)
// into a Map<id, Map<field, String>>. Blank cells are dropped so they fall back
// to the automatic value downstream.
def parseManualParams = { path ->
    def lines = file(path).readLines().findAll { it?.trim() }
    if (!lines) return [:]
    def header = lines[0].split(',', -1)*.trim()
    def idIdx = header.indexOf('id')
    if (idIdx < 0) {
        log.error "manual_params file '${path}' has no 'id' column"
        System.exit(1)
    }
    def out = [:]
    lines.drop(1).each { line ->
        def cells = line.split(',', -1)*.trim()
        def id = cells[idIdx]
        if (!id) return
        def row = [:]
        header.eachWithIndex { col, i ->
            if (col != 'id' && i < cells.size() && cells[i]) row[col] = cells[i]
        }
        out[id] = row
    }
    return out
}

// ── Module includes ─────────────────────────────────────────────────────────
include { PARSE_AND_CLEAN_INPUT } from './modules/input_parsing/main'
include { PARSE_METADATA        } from './modules/input_parsing/main'
include { GENERATE_MANIFEST     } from './modules/input_parsing/main'
include { AGGREGATE_QC_REPORTS  } from './modules/denoise_qc/main'
include { EXPORT_BIOM           } from './modules/downstream/main'
include { EXPORT_TAXONOMY       } from './modules/downstream/main'
include { FILTER_BY_PREVALENCE  } from './modules/downstream/main'
include { ALPHA_DIVERSITY       } from './modules/downstream/main'
include { ALPHA_GROUP_SIGNIFICANCE } from './modules/downstream/main'
include { BETA_DIVERSITY        } from './modules/downstream/main'
include { BETA_GROUP_SIGNIFICANCE  } from './modules/downstream/main'
include { LACTOBACILLUS_ABUNDANCE  } from './modules/downstream/main'

// ── Subworkflow includes ────────────────────────────────────────────────────
include { CLASSIFIER_PREPARATION } from './subworkflows/classifier_preparation'
include { PER_STUDY_PROCESSING          } from './subworkflows/per_study_processing'
include { PER_STUDY_PROCESSING_DEBLUR   } from './subworkflows/per_study_processing_deblur'


// ── Help message ────────────────────────────────────────────────────────────
def helpMessage() {
    log.info"""
    =========================================
     REPROMB ATLAS 16S PIPELINE  v${workflow.manifest.version}
    =========================================

    Usage:
      nextflow run main.nf [options]

    Required:
      --input_csv           Path to the main study table (CSV)
      --data_dir            Base directory with fetchngs results
      --silva_seqs          SILVA reference sequences QZA
      --silva_tax           SILVA reference taxonomy QZA
      --study_metadata_dir  Directory containing per-study Excel metadata

    Optional:
      --outdir              Output directory [default: ./results]
      --denoise_mode        'auto' (per-study from quality stats, default) |
                            'manual' (per-study values from --manual_params)
      --manual_params       Path to manual_params.csv (required when
                            --denoise_mode manual). Generate it first with:
                              bin/prompt_denoise_params.py <input_csv> -o manual_params.csv
      --denoise_auto        Auto-determine denoise params [default: true]
      --classifier_cache_dir  Pre-built classifiers directory
      --run_diversity       Run diversity analysis [default: true]
      --help                Show this message

    Denoiser:
      --denoiser            'dada2' (default) | 'deblur' (studies whose FASTQ
                            files carry synthetic, flat quality scores)
      --deblur_trim_length          Global trim length (bp). Set per region instead.
      --deblur_trim_length_by_region  Map region→bp, single value per region group,
                            e.g. --deblur_trim_length_by_region.V3-V4 250
      --deblur_forward_regions      Regions processed forward-only (no merge)
      --deblur_min_reads    Deblur abundance gate [default: 50, > QIIME2's 10]

    Manual per-study denoise (two steps):
      1) bin/prompt_denoise_params.py <input_csv> -o manual_params.csv
      2) nextflow run main.nf ... --denoise_mode manual --manual_params manual_params.csv

    Profiles:
      -profile docker       Run with Docker containers
      -profile singularity  Run with Singularity containers
      -profile slurm        Submit to SLURM scheduler
      -profile test         Use reduced resources for testing
    """.stripIndent()
}


// ── Parameter validation ────────────────────────────────────────────────────
if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input_csv) {
    log.error "ERROR: --input_csv is required"
    helpMessage()
    exit 1
}

if (!params.data_dir) {
    log.error "ERROR: --data_dir is required"
    helpMessage()
    exit 1
}

if (!params.silva_seqs || !params.silva_tax) {
    log.error "ERROR: --silva_seqs and --silva_tax are required"
    helpMessage()
    exit 1
}

if (!params.study_metadata_dir) {
    log.error "ERROR: --study_metadata_dir is required"
    helpMessage()
    exit 1
}

if (!(params.denoiser in ['dada2', 'deblur'])) {
    log.error "ERROR: --denoiser must be 'dada2' or 'deblur' (got '${params.denoiser}')"
    exit 1
}

if (params.denoise_mode == 'manual' && !params.manual_params) {
    log.error """ERROR: --denoise_mode manual requires --manual_params <manual_params.csv>.
    Generate it first by running:  bin/prompt_denoise_params.py <input_csv> -o manual_params.csv"""
    exit 1
}

// Load the per-study manual override sheet once (empty map in auto mode).
manual_overrides = (params.denoise_mode == 'manual' && params.manual_params) \
    ? parseManualParams(params.manual_params) : [:]


// ── Main workflow ───────────────────────────────────────────────────────────

workflow {

    log.info"""
    =========================================
     REPROMB ATLAS 16S PIPELINE  v${workflow.manifest.version}
    =========================================
    input_csv     : ${params.input_csv}
    data_dir      : ${params.data_dir}
    outdir        : ${params.outdir}
    denoiser      : ${params.denoiser}
    denoise_auto  : ${params.denoise_auto}
    silva_seqs    : ${params.silva_seqs}
    silva_tax     : ${params.silva_tax}
    =========================================
    """.stripIndent()

    // ── Stage input files ───────────────────────────────────────────────
    input_csv_ch  = Channel.fromPath(params.input_csv, checkIfExists: true)
    data_dir_ch   = Channel.fromPath(params.data_dir, checkIfExists: true, type: 'dir')
    silva_seqs_ch = Channel.fromPath(params.silva_seqs, checkIfExists: true)
    silva_tax_ch  = Channel.fromPath(params.silva_tax, checkIfExists: true)
    study_metadata_dir_ch = Channel.fromPath(params.study_metadata_dir, checkIfExists: true, type: 'dir')


    // ═══════════════════════════════════════════════════════════════════
    // STEP 1: Parse and clean input table
    // ═══════════════════════════════════════════════════════════════════
    PARSE_AND_CLEAN_INPUT(input_csv_ch)

    // Parse the cleaned CSV into a channel of study entries
    // Group reproductive tract sites from the same study+primer into one unit
    studies_ch = PARSE_AND_CLEAN_INPUT.out.cleaned_csv
        .splitCsv(header: true, quote: '"')
        .map { row ->
            def study = row.StudyProject?.trim()
            def site = row.Site?.trim()
            def region = row.Region?.trim()
            def paired = row['Paired/Single']?.trim()
            def fwd = row.ForwardPrimerCanonical?.trim() ?: row.ForwardPrimerClean?.trim()
            def rev = row.ReversePrimerCanonical?.trim() ?: row.ReversePrimerClean?.trim()
            def read_len = row.ReadLength?.trim() ?: ''
            def primer_region_id = row.PrimerRegionID?.trim()
            def notes = row.Notes?.trim() ?: ''

            // Determine site group
            def repro_sites = ['Endometrium', 'Vagina', 'Cervix', 'Cervicovaginal']
            def site_group = site in repro_sites ? 'reproductive' : 'gut'

            // Create meta map
            def meta = [
                study_project:    study,
                site:             site,
                site_group:       site_group,
                region:           region,
                paired_single:    paired,
                fwd_primer:       fwd,
                rev_primer:       rev,
                read_length:      read_len,
                primer_len_f:     fwd?.length() ?: 0,
                primer_len_r:     rev?.length() ?: 0,
                primer_region_id: primer_region_id,
                notes:            notes,
                sample_filter:    null,
            ]

            return [study, site, site_group, region, paired, fwd, rev, primer_region_id, meta]
        }
        // Group reproductive tract sites from same study+primer
        .groupTuple(by: [0, 2, 3, 4, 5, 6, 7])
        .map { study, sites, site_group, region, paired, fwd, rev, primer_region_id, metas ->
            // Create combined meta
            def combined_sites = sites instanceof List ? sites : [sites]
            def id = "${study}_${site_group == 'reproductive' ? combined_sites.sort().join('_') : combined_sites[0]}"

            // Manual per-study denoise overrides (null in auto mode, or when this
            // id is absent from the sheet → graceful fallback to auto).
            // id is a GString; coerce to String for the map lookup (GString keys
            // do not match String keys in Groovy maps).
            def denoise_override = manual_overrides[id.toString()] ?: null
            if (params.denoise_mode == 'manual' && !denoise_override) {
                log.warn "Study ${id}: no row in --manual_params; falling back to automatic denoise params."
            }

            def meta = [
                id:               id,
                study_project:    study,
                sites:            combined_sites,
                site_group:       site_group,
                region:           region,
                paired_single:    paired,
                fwd_primer:       fwd,
                rev_primer:       rev,
                fwd_primer_rc:    revcomp(fwd),
                rev_primer_rc:    revcomp(rev),
                primer_len_f:     fwd?.length() ?: 0,
                primer_len_r:     rev?.length() ?: 0,
                primer_region_id: primer_region_id,
                denoise_override: denoise_override,
                sample_filter:    null,
            ]

            return [meta, combined_sites]
        }


    // ═══════════════════════════════════════════════════════════════════
    // STEP 2: Extract metadata and Generate QIIME2 manifests
    // ═══════════════════════════════════════════════════════════════════
    PARSE_METADATA(studies_ch, study_metadata_dir_ch.collect())

    manifest_input = PARSE_METADATA.out.sample_filter
        .join(studies_ch, by: 0) // join by meta.id
        .map { meta, filter_txt, sites -> 
            return [meta, sites, filter_txt]
        }

    GENERATE_MANIFEST(manifest_input, data_dir_ch.collect())

    // Read the actual mode from the .mode file emitted by generate_manifest.py.
    // This handles auto-fallback from PAIRED→SINGLE (e.g. Ion Torrent merged reads).
    manifest_with_mode = GENERATE_MANIFEST.out.manifest
        .map { meta, manifest, mode_file ->
            def actual_mode = mode_file.text.trim()
            if (actual_mode != meta.paired_single) {
                log.warn "Study ${meta.id}: mode changed from ${meta.paired_single} to ${actual_mode} (auto-detected)"
            }
            def updated_meta = meta + [paired_single: actual_mode]
            return [updated_meta, manifest]
        }


    // ═══════════════════════════════════════════════════════════════════
    // STEP 3: Build classifiers (parallel with study processing)
    // ═══════════════════════════════════════════════════════════════════
    // Extract unique primer combinations (already deduplicated by clean_primers.py)
    primer_combos_ch = studies_ch
        .map { meta, sites ->
            def primer_meta = [
                region_id:   meta.primer_region_id,
                region:      meta.region,
                fwd_primer:  meta.fwd_primer,
                rev_primer:  meta.rev_primer,
            ]
            return primer_meta
        }
        .unique { it.region_id }

    CLASSIFIER_PREPARATION(
        primer_combos_ch,
        silva_seqs_ch,
        silva_tax_ch
    )


    // ═══════════════════════════════════════════════════════════════════
    // STEP 4: Per-study processing pipeline
    // ═══════════════════════════════════════════════════════════════════
    // Route on params.denoiser. Both branches share import, primer removal,
    // classifiers, taxonomy and export; only the denoising core differs.
    if (params.denoiser == 'deblur') {
        PER_STUDY_PROCESSING_DEBLUR(
            manifest_with_mode,
            CLASSIFIER_PREPARATION.out.classifiers
        )
        genus_tables_ch = PER_STUDY_PROCESSING_DEBLUR.out.genus_tables
        taxonomy_ch     = PER_STUDY_PROCESSING_DEBLUR.out.taxonomy

        // STEP 5: no DADA2-style QC aggregation — Deblur emits DeblurStats, not
        // the filtered/merged/non-chimeric schema check_denoise_stats.py reads.
        // The guardrail is min-reads + downstream prevalence filtering.
    } else {
        PER_STUDY_PROCESSING(
            manifest_with_mode,
            CLASSIFIER_PREPARATION.out.classifiers
        )
        genus_tables_ch = PER_STUDY_PROCESSING.out.genus_tables
        taxonomy_ch     = PER_STUDY_PROCESSING.out.taxonomy

        // STEP 5: Aggregate QC reports (DADA2 path only)
        AGGREGATE_QC_REPORTS(
            PER_STUDY_PROCESSING.out.qc_reports
                .map { meta, report -> report }
                .collect(),
            PER_STUDY_PROCESSING.out.qc_warnings
                .map { meta, warnings -> warnings }
                .collect()
        )
    }


    // ═══════════════════════════════════════════════════════════════════
    // STEP 6: Downstream analysis (denoiser-agnostic)
    // ═══════════════════════════════════════════════════════════════════
    // Export tables (genus_tables are already bacteria-only — filtered at ASV/sOTU level upstream)
    EXPORT_BIOM(genus_tables_ch)
    EXPORT_TAXONOMY(taxonomy_ch)

    // Prevalence filtering
    FILTER_BY_PREVALENCE(genus_tables_ch)

    // Diversity (optional)
    ALPHA_DIVERSITY(FILTER_BY_PREVALENCE.out.filtered_table)
    BETA_DIVERSITY(FILTER_BY_PREVALENCE.out.filtered_table)

    // Join with metadata for stratified analysis.
    // Keyed on meta.id, NOT a plain .join() on the full meta map. PARSE_METADATA
    // emitted metadata_tsv against the PRE-mode-correction meta (paired_single
    // can still say PAIRED there), while every channel below carries the
    // CORRECTED meta from manifest_with_mode (paired_single auto-fixed to SINGLE
    // when GENERATE_MANIFEST detects a PAIRED→SINGLE fallback, e.g. merged/Ion
    // Torrent submissions). A bare .join() keys on full first-tuple-element
    // equality, so any study whose mode got auto-corrected silently drops out of
    // every downstream diversity/Lactobacillus join — no error, just missing
    // output (seen in PRJNA1079033_Endometrium_Vagina). The exported tables were
    // never affected: EXPORT_BIOM/EXPORT_TAXONOMY do not use this join.
    metadata_tsv_by_id = PARSE_METADATA.out.metadata_tsv
        .map { meta, tsv -> [meta.id, tsv] }

    alpha_sig_input = ALPHA_DIVERSITY.out.shannon
        .join(ALPHA_DIVERSITY.out.observed)
        .join(ALPHA_DIVERSITY.out.simpson)
        .map { meta, shannon, observed, simpson -> [meta.id, meta, shannon, observed, simpson] }
        .join(metadata_tsv_by_id, by: 0)
        .map { id, meta, shannon, observed, simpson, metadata_tsv -> [meta, shannon, observed, simpson, metadata_tsv] }

    ALPHA_GROUP_SIGNIFICANCE(alpha_sig_input)

    beta_sig_input = BETA_DIVERSITY.out.braycurtis
        .join(BETA_DIVERSITY.out.jaccard)
        .map { meta, braycurtis, jaccard -> [meta.id, meta, braycurtis, jaccard] }
        .join(metadata_tsv_by_id, by: 0)
        .map { id, meta, braycurtis, jaccard, metadata_tsv -> [meta, braycurtis, jaccard, metadata_tsv] }

    BETA_GROUP_SIGNIFICANCE(beta_sig_input)

    // Lactobacillus abundance
    lacto_input = EXPORT_BIOM.out.tsv
        .join(EXPORT_TAXONOMY.out.taxonomy_tsv)
        .map { meta, tsv, taxonomy_tsv -> [meta.id, meta, tsv, taxonomy_tsv] }
        .join(metadata_tsv_by_id, by: 0)
        .map { id, meta, tsv, taxonomy_tsv, metadata_tsv -> [meta, tsv, taxonomy_tsv, metadata_tsv] }

    LACTOBACILLUS_ABUNDANCE(lacto_input)
}


// ── Completion handler ──────────────────────────────────────────────────────
workflow.onComplete {
    log.info"""
    =========================================
     REPROMB ATLAS 16S PIPELINE — COMPLETE
    =========================================
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    =========================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
