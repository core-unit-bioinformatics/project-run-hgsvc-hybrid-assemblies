
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


rule run_all_asm_ava_align:
    input:
        tsv_vrk = expand(
            rules.normalize_ava_align_vrk_paf.output.tsv,
            sample=PLAIN_SAMPLES
        ),
        tsv_hsm = expand(
            rules.normalize_ava_align_hsm_paf.output.tsv,
            sample=PLAIN_SAMPLES
        ),
