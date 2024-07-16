
localrules: split_breakpoints_by_sample
rule split_breakpoints_by_sample:
    input:
        breakpoints = DIR_ANNOTATIONS.joinpath(
            "external", "hgsvc3_verkko_breakpoints.tsv"
        ).resolve(strict=True)
    output:
        tsv = DIR_PROC.joinpath(
            "asm_label_qc", "sseq_breaks",
            "{sample}.sseq-switch-breaks.tsv"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.breakpoints, sep="\t", header=0)
        df.rename({"window_min": "start_hpc", "window_max": "end_hpc"}, axis=1, inplace=True)
        df["sample"] = df["sample"] + ".vrk-ps-sseq"
        df = df.loc[df["sample"] == wildcards.sample, :].copy()
        if df.empty:
            with open(output.tsv, "w") as dump:
                pass
            with open(str(output.tsv) + ".EMPTY", "w") as dump:
                pass
        else:
            df.sort_values(["unitig", "start_hpc", "end_hpc"], inplace=True)
            df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: translate_hpc_breakpoint_coordinates
rule translate_hpc_breakpoint_coordinates:
    input:
        tsv = rules.split_breakpoints_by_sample.output.tsv,
        cmap = lambda wildcards: WORKDIR_ASSEMBLY.joinpath(
            "results/assemblies/verkko",
            f"{wildcards.sample.split('.')[0]}.ps-sseq",
            "aux",
            f"{wildcards.sample.split('.')[0]}.ps-sseq.graph-linear-hpc-map.tsv.gz"
        ),
        gsize = rules.prep_verkko_fasta_index.output.sizes
    output:
        bed = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "sseq_breaks", "{sample}.sseq-switch-breaks.bed"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import pathlib as pl
        import io

        out_columns = [
            "#seq", "start", "end", "name", "score",
            "unitig", "hpc_start", "hpc_end"
        ]

        df = pd.read_csv(input.tsv, sep="\t", header=0)
        if df.empty:
            buffer = io.StringIO()
            buffer.write("\t".join(out_columns) + "\n")
            with open(input.gsize) as listing:
                for line in listing:
                    seq, size = line.strip().split()
                    buffer.write(
                        f"{seq}\t0\t{size}\tMOCK\t0\tMOCK\t0\t0\n"
                    )
            with open(output.bed, "w") as dump:
                dump.write(buffer.getvalue())
            with open(str(output.bed)+".MOCK", "w"):
                pass
        else:
            cmap = pd.read_csv(input.cmap, sep="\t", header=0)
            cmap = cmap.loc[cmap["target_name"].isin(df["unitig"].unique()), :].copy()
            cmap = cmap.loc[cmap["tp_align_type"] != 2, :].copy()

            output = []
            for row in df.itertuples():
                select_tig = cmap["query_name"] == row.unitig
                select_start = row.start_hpc < cmap["query_end"]
                select_end = row.end_hpc >= cmap["query_start"]
                selector = select_tig & select_start & select_end
                if not selector.any():
                    raise RuntimeError(f"No coord mapping for {wildcards.sample} / {row}")
                cmap_sub = cmap.loc[selector, :]
                assert cmap_sub.shape[0] == 1
                # assume breakpoint is always somewhere "in the middle"
                # of the contig and the graph-to-linear alignment is surrounding
                # this location. Hence, the coordinates in cmap should
                # be smaller/larger
                # --- larger minus smaller
                offset_start = row.start_hpc - cmap_sub["query_start"]
                assert offset_start >= 0
                # --- larger minus smaller
                offset_end = cmap_sub["query_end"] - row.end_hpc
                assert offset_end >= 0

                expanded_start = cmap_sub["target_start_plain"] + offset_start
                expanded_end = cmap_sub["target_end_plain"] - offset_end
                assert expanded_start > expanded_end
                output.append(
                    (
                        cmap_sub["target_name"].iloc[0], expanded_start, expanded_end,
                        "SSEQBRKP", 1, row.unitig, row.start_hpc, row.end_hpc
                    )
                )
            output = pd.DataFrame.from_records(
                output, columns=out_columns
            )

            output.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_sseq_breakpoints:
    input:
        bed = expand(
            rules.translate_hpc_breakpoint_coordinates.output.bed,
            sample=SAMPLES
        )
