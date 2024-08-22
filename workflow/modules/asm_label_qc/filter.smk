
localrules: subset_merged_issues_oneplus
rule subset_merged_issues_oneplus:
    input:
        table = rules.add_ngap_sizes.output.table
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.tsv.gz"
        )
        error_regions = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.errors.bed.gz"
        )
    run:
        import pandas as pd

        def count_distinct_labels(labels):
            distinct = set(l.split("::")[0] for l in labels.split(",")) - set(["no-labels"])
            return len(distinct)

        df = pd.read_csv(input.table, sep="\t", header=0)
        df["distinct_labels"] = df["labels"]
        sub = df.loc[df["distinct_labels"] > 1, :].copy()
        sub.to_csv(output.table, sep="\t", header=True, index=False)

        sub = sub[["seq", "start", "end"]].copy()
        sub.rename({"seq": "#seq"}, axis=1, inplace=True)
        sub.to_csv(error_regions, sep="\t", header=True, index=False)

        df.drop_duplicates("seq", inplace=True)
        df.insert(1, "start", 0)
        df = df[["seq", "start", "seq_length"]]
        df.rename({"seq": "#seq", "seq_length": "end"}, axis=1, inplace=True)
        df.to_csv(output.whole_genome, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


FILTER_GENOME = {
    "ps-no-ont": " | grep -v unassigned | ",
    "wg-no-ont": " | "
}


rule subtract_error_regions_from_genome:
    input:
        sizes = DIR_PROC.joinpath(
            "asm_label_qc", "assembly_size",
            "{sample}.sizes.txt"
        ),
        err = rules.subset_merged_issues_oneplus.output.error_regions
    output:
        clean = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.complement.bed.gz"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        filter_seq = lambda wildcards: FILTER_GENOME[wildcards.span]
    shell:
        "bedtools complement -i {input.err} -g {input.sizes}"
        "{params.filter_seq} gzip > {output.clean}"


localrules: summarize_region_stats
rule summarize_region_stats:
    input:
        regions = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.{regionset}.bed.gz"
        )
    output:
        table = clean = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.{regionset}.stats.tsv"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.regions, sep="\t", header=None, names=["seq", "start", "end"])
        num_regions = df.shape[0]
        df["length"] = df["end"] - df["start"]
        region_dist = df["length"].describe([0, 0.25, 0.50, 0.75, 1])
        print(region_dist)
        raise


rule run_all_summarize_error_regions:
    input:
        tables = expand(
            rules.summarize_region_stats.output.table,
            sample=SAMPLES,
            span=["ps-no-ont", "wg-no-ont"],
            regionset=["error", "complement"]
        )
