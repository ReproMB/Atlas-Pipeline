/*
 * Subworkflow: classifier_preparation
 * Builds classifiers for each unique primer combination.
 * Supports caching: if classifier exists in cache dir, skip build.
 */

include { EXTRACT_READS   } from '../modules/classifier_build/main'
include { FIT_CLASSIFIER  } from '../modules/classifier_build/main'


workflow CLASSIFIER_PREPARATION {

    take:
    primer_combos_ch    // channel: [primer_meta]
    silva_seqs          // path: SILVA sequences QZA
    silva_tax           // path: SILVA taxonomy QZA

    main:

    // ── Check for cached classifiers ────────────────────────────────────
    if (params.classifier_cache_dir) {
        primer_combos_ch
            .branch { primer_meta ->
                def cache_file = file("${params.classifier_cache_dir}/classifier_${primer_meta.region_id}.qza")
                cached:    cache_file.exists()
                    return [primer_meta, cache_file]
                to_build:  true
                    return primer_meta
            }
            .set { branched }

        cached_classifiers = branched.cached
        to_build = branched.to_build
    } else {
        cached_classifiers = Channel.empty()
        to_build = primer_combos_ch
    }

    // ── Build new classifiers ───────────────────────────────────────────
    // Combine each primer_meta with the SILVA seqs
    to_build_with_seqs = to_build.combine(silva_seqs)

    // Extract reference reads for each primer pair
    EXTRACT_READS(to_build_with_seqs)

    // Combine extracted reads with taxonomy for fitting
    fit_input = EXTRACT_READS.out.ref_reads
        .combine(silva_tax)

    FIT_CLASSIFIER(fit_input)

    // ── Merge cached + newly built ──────────────────────────────────────
    all_classifiers = cached_classifiers
        .mix(FIT_CLASSIFIER.out.classifier)

    emit:
    classifiers = all_classifiers  // channel: [primer_meta, classifier.qza]
}
