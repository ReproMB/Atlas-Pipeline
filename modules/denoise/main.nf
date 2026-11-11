/*
 * Module: denoise
 * DADA2 denoising for paired-end and single-end data.
 * Supports auto-parametrization from quality stats JSON or manual overrides.
 */

process DADA2_DENOISE_PAIRED {
    tag "${meta.id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(demux_qza), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_table.qza"),      emit: table
    tuple val(meta), path("${meta.id}_rep_seqs.qza"),   emit: rep_seqs
    tuple val(meta), path("${meta.id}_stats.qza"),      emit: stats

    script:
    // Determine if we use auto params or manual overrides
    def use_auto = params.denoise_auto && !params.trunc_len_f && !params.trunc_len_r
    // Per-study manual overrides (denoise_mode == 'manual'); null otherwise.
    def ov = meta.denoise_override
    def use_manual = (ov != null)
    // When cutadapt removed primers upstream, do not trim again (trim_left = 0).
    def trim_f_default = params.run_primer_trim ? 0 : (meta.primer_len_f ?: 0)
    def trim_r_default = params.run_primer_trim ? 0 : (meta.primer_len_r ?: 0)
    // max_ee / trunc_q: per-study override if provided, else global.
    def ee = (ov?.max_ee) ?: params.max_ee
    def tq = (ov?.trunc_q) ?: params.trunc_q
    """
    #!/bin/bash
    set -euo pipefail

    # ── Determine denoise parameters ─────────────────────────────────────
    if [ "${use_manual}" = "true" ]; then
        # Manual per-study mode. Per field: sheet value → auto JSON → global → 0.
        PARAMS=\$(python3 -c "
import json
ov = {'trunc_len_f': '${ov?.trunc_len_f ?: ''}', 'trunc_len_r': '${ov?.trunc_len_r ?: ''}', 'trim_left_f': '${ov?.trim_left_f ?: ''}', 'trim_left_r': '${ov?.trim_left_r ?: ''}', 'min_overlap': '${ov?.min_overlap ?: ''}'}
q = {}
try:
    q = json.load(open('${quality_json}'))
except Exception:
    pass
def pick(key, jkey, glob):
    if ov.get(key, '') != '':
        return ov[key]
    if q.get(jkey) is not None:
        return q.get(jkey)
    return glob
trunc_f = pick('trunc_len_f', 'trunc_len_f', ${params.trunc_len_f ?: 0})
trunc_r = pick('trunc_len_r', 'trunc_len_r', ${params.trunc_len_r ?: 0})
trim_f = pick('trim_left_f', 'trim_left_f', ${trim_f_default})
trim_r = pick('trim_left_r', 'trim_left_r', ${trim_r_default})
min_overlap = pick('min_overlap', 'recommended_min_overlap', ${params.min_overlap})
print(f'{trunc_f} {trunc_r} {trim_f} {trim_r} {min_overlap}')
")
    elif [ "${use_auto}" = "true" ] && [ -f "${quality_json}" ]; then
        PARAMS=\$(python3 -c "
import json
q = json.load(open('${quality_json}'))
trunc_f = q.get('trunc_len_f', 0)
trunc_r = q.get('trunc_len_r', 0)
trim_f = q.get('trim_left_f', 0)
trim_r = q.get('trim_left_r', 0)
min_overlap = q.get('recommended_min_overlap', ${params.min_overlap})
print(f'{trunc_f} {trunc_r} {trim_f} {trim_r} {min_overlap}')
")
    else
        PARAMS="${params.trunc_len_f ?: 0} ${params.trunc_len_r ?: 0} ${params.trim_left_f ?: trim_f_default} ${params.trim_left_r ?: trim_r_default} ${params.min_overlap}"
    fi

    read TRUNC_F TRUNC_R TRIM_F TRIM_R MIN_OVERLAP <<< \$PARAMS

    echo "Denoise params [${use_manual ? 'manual' : (use_auto ? 'auto' : 'global')}]: trunc_f=\$TRUNC_F trunc_r=\$TRUNC_R trim_f=\$TRIM_F trim_r=\$TRIM_R min_overlap=\$MIN_OVERLAP max_ee=${ee} trunc_q=${tq}"

    # ── Run DADA2 ────────────────────────────────────────────────────────
    qiime dada2 denoise-paired \\
        --i-demultiplexed-seqs ${demux_qza} \\
        --p-trunc-len-f \$TRUNC_F \\
        --p-trunc-len-r \$TRUNC_R \\
        --p-trim-left-f \$TRIM_F \\
        --p-trim-left-r \$TRIM_R \\
        --p-max-ee-f ${ee} \\
        --p-max-ee-r ${ee} \\
        --p-trunc-q ${tq} \\
        --p-min-overlap \$MIN_OVERLAP \\
        --p-max-merge-mismatch ${params.max_merge_mismatch} \\
        --p-pooling-method ${params.pooling_method} \\
        --p-chimera-method ${params.chimera_method} \\
        --p-n-threads ${task.cpus} \\
        --p-n-reads-learn ${params.n_reads_learn} \\
        --p-hashed-feature-ids true \\
        --o-table ${meta.id}_table.qza \\
        --o-representative-sequences ${meta.id}_rep_seqs.qza \\
        --o-denoising-stats ${meta.id}_stats.qza
    """
}


// NB: qiime dada2 denoise-single accepts both SampleData[SequencesWithQuality]
// and SampleData[PairedEndSequencesWithQuality]; given a paired artifact it uses
// only the forward reads. The adaptive retry relies on this for the forward-only
// fallback (a paired study re-run single-end without re-importing).
process DADA2_DENOISE_SINGLE {
    tag "${meta.id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(demux_qza), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_table.qza"),      emit: table
    tuple val(meta), path("${meta.id}_rep_seqs.qza"),   emit: rep_seqs
    tuple val(meta), path("${meta.id}_stats.qza"),      emit: stats

    script:
    def use_auto = params.denoise_auto && !params.trunc_len_f
    // Per-study manual overrides (denoise_mode == 'manual'); null otherwise.
    def ov = meta.denoise_override
    def use_manual = (ov != null)
    // When cutadapt removed primers upstream, do not trim again (trim_left = 0).
    def trim_f_default = params.run_primer_trim ? 0 : (meta.primer_len_f ?: 0)
    // max_ee / trunc_q: per-study override if provided, else global.
    def ee = (ov?.max_ee) ?: params.max_ee
    def tq = (ov?.trunc_q) ?: params.trunc_q
    """
    #!/bin/bash
    set -euo pipefail

    # ── Determine denoise parameters ─────────────────────────────────────
    if [ "${use_manual}" = "true" ]; then
        # Manual per-study mode. Per field: sheet value → auto JSON → global → 0.
        PARAMS=\$(python3 -c "
import json
ov = {'trunc_len_f': '${ov?.trunc_len_f ?: ''}', 'trim_left_f': '${ov?.trim_left_f ?: ''}'}
q = {}
try:
    q = json.load(open('${quality_json}'))
except Exception:
    pass
def pick(key, jkey, glob):
    if ov.get(key, '') != '':
        return ov[key]
    if q.get(jkey) is not None:
        return q.get(jkey)
    return glob
trunc = pick('trunc_len_f', 'trunc_len_f', ${params.trunc_len_f ?: 0})
trim = pick('trim_left_f', 'trim_left_f', ${trim_f_default})
print(f'{trunc} {trim}')
")
    elif [ "${use_auto}" = "true" ] && [ -f "${quality_json}" ]; then
        PARAMS=\$(python3 -c "
import json
q = json.load(open('${quality_json}'))
trunc = q.get('trunc_len_f', 0)
trim = q.get('trim_left_f', 0)
print(f'{trunc} {trim}')
")
    else
        PARAMS="${params.trunc_len_f ?: 0} ${params.trim_left_f ?: trim_f_default}"
    fi

    read TRUNC TRIM <<< \$PARAMS

    echo "Denoise params [${use_manual ? 'manual' : (use_auto ? 'auto' : 'global')}]: trunc=\$TRUNC trim=\$TRIM max_ee=${ee} trunc_q=${tq}"

    # ── Run DADA2 ────────────────────────────────────────────────────────
    qiime dada2 denoise-single \\
        --i-demultiplexed-seqs ${demux_qza} \\
        --p-trunc-len \$TRUNC \\
        --p-trim-left \$TRIM \\
        --p-max-ee ${ee} \\
        --p-trunc-q ${tq} \\
        --p-pooling-method ${params.pooling_method} \\
        --p-chimera-method ${params.chimera_method} \\
        --p-n-threads ${task.cpus} \\
        --p-n-reads-learn ${params.n_reads_learn} \\
        --p-hashed-feature-ids true \\
        --o-table ${meta.id}_table.qza \\
        --o-representative-sequences ${meta.id}_rep_seqs.qza \\
        --o-denoising-stats ${meta.id}_stats.qza
    """
}
