
rule normalize_flagger_annotation:
    input:
        folder = FLAGGER_ROOT_FOLDER
    output:
        bed_like = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "flagger", "{sample}.flagger-labels.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pathlib as pl
        import pandas as pd
        label_map = {"Unk": 0, "Hap": 1, "Dup": 2, "Col": 3, "Err": 4}

        sample_id, assembler = wildcards.sample.split(".")
        assembler = {"vrk-ps-sseq": "verkko", "hsm-ps-sseq": "hifiasm"}[assembler]
        flagger_file = pl.Path(input.folder).joinpath(
            f"{assembler}", f"{sample_id}.{assembler}.alt_removed.flagger_final.bed"
        )
        if not flagger_file.is_file():
            if sample_id.startswith("NA"):
                sample_id = sample_id.replace("NA", "GM")
                flagger_file = pl.Path(input.folder).joinpath(
                    f"{assembler}", f"{sample_id}.{assembler}.alt_removed.flagger_final.bed"
                )
                if not flagger_file.is_file():
                    raise FileNotFoundError(flagger_file)
        table_header = ["chrom", "start", "end", "name", "score", "strand", "start2", "end2", "rgb"]
        use_header = table_header[:4]
        df = pd.read_csv(flagger_file, sep="\t", skiprows=1, header=None, names=table_header, usecols=use_header)
        df["score"] = df["name"].replace(label_map)
        df["score"] = df["score"].astype(int)
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)


rule run_normalize_flagger_results:
    input:
        tables = expand(
            rules.normalize_flagger_annotation.output.bed_like,
            sample=SAMPLES
        )
