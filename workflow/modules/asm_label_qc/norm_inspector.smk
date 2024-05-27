
rule normalize_inspector_output:
    input:
        folder = INSPECTOR_ROOT_FOLDER
    output:
        bed_like = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "inspector", "{sample}.inspector-errors.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import pathlib as pl
        id_only, assm = wildcards.sample.split(".")

        assembler = {"vrk-ps-sseq": "Verkko", "hsm-ps-sseq": "hifiasm"}[assm]
        error_files = pl.Path(input.folder).joinpath(
            f"{assembler}", f"{id_only}"
        ).glob(f"{id_only}_{read_type}*.bed")

        def norm_coordinate(coordinate, start):
            try:
                c = int(coordinate)
            except ValueError:
                if not ";" in coordinate:
                    raise
                coordinates = [int(c) for c in coordinate.split(";")]
                if start:
                    c = min(coordinates)
                else:
                    c = max(coordinates)
            return c

        concat = []
        for error_file in error_files:
            error_type = "struct" if "structural" in error_file.name else "small"
            read_type = "hifi" if "HiFi" in error_file.name else "ont"
            df = pd.read_csv(error_file, sep="\t", header=0)
            df["name"] = f"inspector_{read_type}_{error_type}_" + df["Type"].str.lower()
            df["start"] = df["Start_Position"].apply(norm_coordinate, args=(True,))
            df["end"] = df["End_Position"].apply(norm_coordinate, args=(False,))
            df["chrom"] = df["#Contig_Name"]
            df = df[["chrom", "start", "end", "name"]].copy()
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["chrom", "start"], ascending=True, inplace=True)
        concat.to_csv(output.bed_like, sep="\t", header=True, index=False)


rule run_normalize_inspector_results:
    input:
        tables = expand(
            rules.normalize_inspector_output.output.bed_like,
            sample=SAMPLES
        )
