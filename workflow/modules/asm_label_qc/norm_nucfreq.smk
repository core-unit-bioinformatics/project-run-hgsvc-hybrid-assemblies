
rule bin_nucfreq_regions_by_coverage:
    input:
        bed = WORKDIR_EVAL.joinpath(
            "results/regions/{sample}",
            "{sample}.nucfreq.covann.tsv.gz"
        )
    output:
        bed = DIR_PROC.joinpath(
            "asm_label_qc", "norm_tables", "nucfreq",
            "{sample}.nucfreq-cov-bin.bed.gz"
        )
    run:
        import pandas as pd
        import sys

        df = pd.read_csv(input.bed, sep="\t", header=0)
        bins = [0, 50, 90, 110, 150, 200, sys.maxsize]
        df["score"] = pd.cut(
            df["hifi_pct_median_cov"].values,
            bins, right=False, include_lowest=True,
            labels=False
        )
        df = df[["contig", "start", "end", "score", "asm_unit", "num_hets"]]
        df.rename({"contig": "#seq"}, axis=1, inplace=True)
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_bin_nucfreq_regions:
    input:
        beds = expand(
            rules.bin_nucfreq_regions_by_coverage.output.bed,
            sample=SAMPLES
        )
