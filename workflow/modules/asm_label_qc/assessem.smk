
rule prepare_region_cache:
    input:
        sizes = DIR_PROC.joinpath("asm_label_qc", "assembly_size", "{sample}.sizes.txt"),
        flagger = rules.normalize_flagger_annotation.output.bed_like,
        inspect_hifi = DIR_RES.joinpath("asm_label_qc", "norm_tables",
            "inspector", "{sample}.hifi.inspector-errors.tsv.gz"
        ),
        inspect_ont = DIR_RES.joinpath("asm_label_qc", "norm_tables",
            "inspector", "{sample}.ont.inspector-errors.tsv.gz"
        ),
        merqury = rules.normalize_merqury_kmer_tracks.output.bed_like,
        hifi_mq00 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.hifi.mq00.mosdepth-windowed.tsv.gz"
        ),
        hifi_mq60 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.hifi.mq60.mosdepth-windowed.tsv.gz"
        ),
        ont_mq00 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.ont.mq00.mosdepth-windowed.tsv.gz"
        ),
        ont_mq60 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.ont.mq60.mosdepth-windowed.tsv.gz"
        ),
        nucfreq = WORKDIR_EVAL.joinpath(
            "results/regions/{sample}",
            "{sample}.nucfreq-flagged.bed.gz"
        )
    output:
        hdf = DIR_RES.joinpath(
            "asm_label_qc", "region_cache",
            "{sample}.assessem-cache.h5"
        )
    log:
        DIR_LOG.joinpath(
            "asm_label_qc", "region_cache", "{sample}.assessem-cache.log"
        )
    conda:
        DIR_ENVS.joinpath("assessem.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_label_qc", "assessem-cache.py"),
        track_files = lambda wildcards, input: get_assessem_cli_parameters(input, "files"),
        track_labels = lambda wildcards, input: get_assessem_cli_parameters(input, "labels"),
        score_columns = lambda wildcards, input: get_assessem_cli_parameters(input, "columns"),
    shell:
        "{params.script} --verbose --data-cache {output.hdf} --genome-size {input.sizes} "
        " --track-files {params.track_files} --track-labels {params.track_labels} "
        " --score-columns {params.score_columns} > {log}"


rule run_all_assessem_jobs:
    input:
        cache = expand(
            rules.prepare_region_cache.output.hdf,
            sample=SAMPLES
        )

