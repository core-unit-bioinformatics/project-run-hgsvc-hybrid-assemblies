
rule align_hifiasm_to_verkko:
    input:
        verkko = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        vrk_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
        hifiasm = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.hsm-ps-sseq.asm-mrg-tag.fasta"
        ),
        hsm_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.hsm-ps-sseq.asm-mrg-tag.fasta.fai"
        )
    output:
        paf = DIR_PROC.joinpath(
            "asm_compare", "ava_align", "raw",
            "{sample}.hsm-to-vrk.ps-sseq.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: 65536 * attempt,
        time_hrs=lambda wildcards, attempt: 4 * attempt
    shell:
        "minimap2 -t {threads} -x asm5 --eqx --MD -c --secondary=no "
        "{input.verkko} {input.hifiasm} | pigz > {output.paf}"


rule normalize_ava_align_vrk_paf:
    input:
        paf = rules.align_hifiasm_to_verkko.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "asm_compare", "ava_align", "norm",
            "{sample}.hsm-to-vrk.ps-sseq.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=DIR_SCRIPTS.joinpath("normalize_paf.py").resolve(strict=True)
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule align_verkko_to_hifiasm:
    input:
        verkko = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        vrk_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
        hifiasm = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.hsm-ps-sseq.asm-mrg-tag.fasta"
        ),
        hsm_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{sample}.hsm-ps-sseq.asm-mrg-tag.fasta.fai"
        )
    output:
        paf = DIR_PROC.joinpath(
            "asm_compare", "ava_align", "raw",
            "{sample}.vrk-to-hsm.ps-sseq.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: 65536 * attempt,
        time_hrs=lambda wildcards, attempt: 4 * attempt
    shell:
        "minimap2 -t {threads} -x asm5 --eqx --MD -c --secondary=no "
        "{input.hifiasm} {input.verkko} | pigz > {output.paf}"


rule normalize_ava_align_hsm_paf:
    input:
        paf = rules.align_verkko_to_hifiasm.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "asm_compare", "ava_align", "norm",
            "{sample}.vrk-to-hsm.ps-sseq.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=DIR_SCRIPTS.joinpath("normalize_paf.py").resolve(strict=True)
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule evaluate_normalized_pafs:
    input:
        paf = DIR_PROC.joinpath(
            "asm_compare", "ava_align", "norm",
            "{sample}.{qry_to_trg}.ps-sseq.norm-paf.tsv.gz"
        ),
        gsize = select_genome_size_file
    output:
        stats = DIR_RES.joinpath(
            "asm_compare", "{qry_to_trg}",
            "statistics",
            "{sample}.{qry_to_trg}.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.stats.tsv"
        ),
        aln_regions = DIR_RES.joinpath(
            "asm_compare", "{qry_to_trg}",
            "regions",
            "{sample}.{qry_to_trg}.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aligned-regions.tsv.gz"
        ),
        sup_regions = DIR_RES.joinpath(
            "asm_compare", "{qry_to_trg}",
            "regions",
            "{sample}.{qry_to_trg}.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.support-regions.tsv.gz"
        ),
        mrg_regions = DIR_RES.joinpath(
            "asm_compare", "{qry_to_trg}",
            "regions",
            "{sample}.{qry_to_trg}.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.merged-regions.tsv.gz"
        )
    log:
        DIR_RES.joinpath(
            "asm_compare", "{qry_to_trg}",
            "statistics",
            "{sample}.{qry_to_trg}.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aln-filter.log"
        )
    wildcard_constraints:
        qry_to_trg="(vrk-to-hsm|hsm-to-vrk)",
        precision="(exact|struct)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=DIR_SCRIPTS.joinpath("eval_paf.py").resolve(strict=True),
        min_aln_len=lambda wildcards: to_int(wildcards.min_aln_len),
        min_seq_len=lambda wildcards: to_int(wildcards.min_seq_len),
        precision=lambda wildcards: parse_precision_param(wildcards.precision, "setting"),
        struct_t=lambda wildcards: parse_precision_param(wildcards.precision, "threshold"),
    shell:
        "{params.script} --verbose --paf-file {input.paf} "
        "--genome-size {input.gsize} "
        "--precision {params.precision} "
        "--retain-mapq-geq {wildcards.min_mapq} "
        "--retain-seqlen-geq {params.min_seq_len} "
        "--retain-alnlen-geq {params.min_aln_len} "
        "--struct-err-geq {params.struct_t} "
        "--split-seq-tag '.' "
        "--out-seq-stats {output.stats} "
        "--out-aligned-regions {output.aln_regions} "
        "--out-support-regions {output.sup_regions} "
        "--out-merged-regions {output.mrg_regions} "
        "&> {log}"


rule run_all_asm_ava_align:
    input:
        stats = expand(
            rules.evaluate_normalized_pafs.output.stats,
            sample=PLAIN_SAMPLES,
            qry_to_trg=["hsm-to-vrk", "vrk-to-hsm"],
            precision=["exact", "struct-50"],
            min_mapq=[1],
            min_seq_len=["0", "100k"],
            min_aln_len=["0", "10k"]
        )
