
rule split_centromere_annotation:
    input:
        listing = CENTROMERE_ANNOTATION
    output:
        beds = expand(
            DIR_RES.joinpath(
                "asm_label_qc", "norm_tables", "centromeres",
                "{sample}.active_asat_HOR_arrays_v2.bed"
            ),
            sample=SAMPLES
        )
    run:
        import pandas as pd
        header = ["sample_plain", "seq_loc", "chrom", "strand"]
        df = pd.read_csv(input.listing, sep="\t", header=None, names=header)
        def get_seq(seq_loc):
            seq, loc = seq_loc.split(":")
            start, end = loc.split("-")
            start = int(start)
            end = int(end)
            assert start < end
            return seq, start, end
        seqs = df["seq_loc"].apply()
        seqs = pd.DataFrame.from_records(seqs, index=df.index, columns=["seq", "start", "end"])
        df = df.join(seqs)
        df.drop("seq_loc", axis=1, inplace=True)
        def assign_assembler(seq):
            if any(x in seq for x in ["h1tg", "h2tg"]):
                return ".hsm-ps-sseq"
            elif any(x in seq for x in ["haplotype1", "haplotype2"]):
                return ".vrk-ps-sseq"
            else:
                raise ValueError(seq)
        df["assembler"] = df["seq"].apply(assign_assembler)
        df["sample"] = df["sample_plain"] + df["assembler"]
        df = df[["seq", "start", "end", "chrom", "strand", "sample"]]
        df.insert(4, "score", 1)
        df.sort_values(["sample", "seq", "start", "end"], inplace=True)
        df.rename({"seq": "#seq"}, axis=1, inplace=True)

        for sample, cens in df.groupby("sample"):
            out_file = DIR_RES.joinpath(
                "asm_label_qc", "norm_tables", "centromeres",
                f"{sample}.active_asat_HOR_arrays_v2.bed"
            )
            cens.to_csv(out_file, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_normalize_centromeres:
    input:
        beds = rules.split_centromere_annotation.output.beds
