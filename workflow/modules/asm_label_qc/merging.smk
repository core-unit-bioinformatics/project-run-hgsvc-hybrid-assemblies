


DEFINED_MERGE_SPANS = {
    "wg": " | ",
    "ps": " | grep -v unassigned | ",
    "wg-no-ont": " | grep -v ISPCON | ",
    "ps-no-ont": " | grep -v unassigned | grep -v ISPCON | ",
}


rule merge_issue_labels:
    input:
        beds = expand(
            DIR_RES.joinpath(
                "asm_label_qc", "merge_tables", "by-sample",
                "{sample}", "{sample}.{annotation}.mrg-labels.bed"
            ),
            annotation=["sseqbreak", "flagger", "nucfreq", "merqury", "busco", "inspector"],
            allow_missing=True
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        grep=lambda wildcards: DEFINED_MERGE_SPANS[wildcards.span]
    shell:
        "cat {input.beds}"
        " {params.grep} "
        "sort -V -k1,1 -k2n,3n"
            " | "
        "bedtools merge -c 4 -o collapse -i /dev/stdin"
            " | "
        "gzip > {output.bed}"


rule run_all_merge_issues:
    input:
        bed = expand(
            rules.merge_issue_labels.output.bed,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        )
