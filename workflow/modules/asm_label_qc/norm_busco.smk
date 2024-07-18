
localrules: normalize_busco_issue_annotation
rule normalize_busco_issue_annotation:
    input:
        bed = WORKDIR_EVAL.joinpath(
            "results/regions", "{sample}",
            "{sample}.busco.primates_odb10.issues.bed"
        )
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "busco",
            "{sample}.busco-issues.bed"
        )
    shell:
        import pandas as pd
        df = pd.read_csv(input.bed, sep="\t", header=0)
        df.drop("gene", axis=1, inplace=True)
        relabel = {
            "Fragmented": 1,
            "Duplicated": 2
        }
        df.insert(4, "score", 1)
        df["score"] = df["label"].replace(relabel, inplace=False)
        df["score"] = df["score"].astype(int)  # capture other potential labels
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_busco:
    input:
        bed = expand(
            rules.normalize_busco_issue_annotation.output.bed,
            sample=SAMPLES
        )

