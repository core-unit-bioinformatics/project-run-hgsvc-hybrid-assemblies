
localrules: subset_merged_issues_oneplus
rule subset_merged_issues_oneplus:
    input:
        table = rules.add_ngap_sizes.output.table
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.tsv.gz"
        ),
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
        df["distinct_labels"] = df["labels"].apply(count_distinct_labels)
        sub = df.loc[df["distinct_labels"] > 1, :].copy()
        sub.to_csv(output.table, sep="\t", header=True, index=False)

        sub = sub[["seq", "start", "end"]].copy()
        sub.rename({"seq": "#seq"}, axis=1, inplace=True)
        sub.to_csv(output.error_regions, sep="\t", header=True, index=False)
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
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.subset-1p.{regionset}.stats.tsv"
        )
    run:
        import pandas as pd
        if wildcards.regionset == "complement":
            header = None
            names = ["seq", "start", "end"]
        else:
            header=0
            names = None
        df = pd.read_csv(input.regions, sep="\t", header=header, names=names)

        prefix = "error" if wildcards.regionset == "errors" else "clean"

        df["length"] = df["end"] - df["start"]
        region_dist = df["length"].describe([0.5, 0.25, 0.50, 0.75, 0.95])
        region_dist.rename(
            {
                "count": f"{prefix}_regions_num",
                "mean": f"{prefix}_region_size_mean",
                "std": f"{prefix}_region_size_stddev",
                "min": f"{prefix}_region_size_min",
                "5%": f"{prefix}_region_size_pct05",
                "25%": f"{prefix}_region_size_pct25",
                "50%": f"{prefix}_region_size_pct50",
                "75%": f"{prefix}_region_size_pct75",
                "95%": f"{prefix}_region_size_pct95",
                "max": f"{prefix}_region_size_max",
            }, inplace=True
        )
        region_dist.insert(0, "sample", wildcards.sample)
        region_dist.to_csv(output.table, sep="\t", header=0, index=False)
    # END OF RUN BLOCK


rule run_all_summarize_error_regions:
    input:
        tables = expand(
            rules.summarize_region_stats.output.table,
            sample=SAMPLES,
            span=["ps-no-ont", "wg-no-ont"],
            regionset=["errors", "complement"]
        )
