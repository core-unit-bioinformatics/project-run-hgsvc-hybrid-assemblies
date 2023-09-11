
rule build_hla_cut_table:
    input:
        norm_paf = expand(
                WORKDIR_EVAL.joinpath(
                    "results/alignments/contig_to_ref/t2tv2/table",
                    "{sample}.asm-{asm_unit}.t2tv2.norm-paf.tsv.gz"
                ),
                sample=SAMPLES,
                asm_unit=["hap1", "hap2", "unassigned"]
            )
    output:
        cut_table = DIR_RES.joinpath(
            "regions", "hla", "assembly_cut_table.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        hla_start = int(28e6),
        hla_end = int(34e6)
    run:
        import pandas as pd
        import pathlib as pl
        cut_table = []

        for paf in input.norm_paf:
            source_file = pl.Path(paf).name
            source_asm_unit = source_file.rsplit(".", 4)[0]
            sample = source_file.split(".asm-")[0]

            df = pd.read_csv(paf, sep="\t", comment="#")
            df = df.loc[(df["target_name"] == "chr6") & (df["tp_align_type"] != 2), :]
            # contigs/alignments reaching less than 100 kbp into the window
            # are just ignored
            select_reach_left = (df["target_end"] - int(1e5)) > params.hla_start
            select_reach_right = (df["target_start"] - int(1e5) < params.hla_end)

            df = df.loc[select_reach_left & select_reach_right, :]
            for query, alignments in df.groupby("query_name"):
                min_q = alignments["query_start"].min()
                max_q = alignments["query_end"].max()
                min_t = alignments["target_start"].min()
                max_t = alignments["target_end"].max()

                query_length = alignments["query_length"].values[0]
                if query_length < int(1e6):
                    # use small query sequences as is, likely garbage anyway
                    cut_table.append(
                        (sample, query, min_q, max_q, min_t, max_t,
                         0, query_length, query_length, source_asm_unit, source_file))
                    continue

                # following: slack of 1 Mbp because rare
                # cases exist where the alignment starts a few
                # 100 kbp into the MHC window; assert that is does
                # not start/end in the middle of the MHC locus;
                # if that happens, need to investigate
                assert params.hla_start + int(1e6) > min_t, f"{sample} / {query} / {min_t}"
                offset_start = params.hla_start - min_t
                assert params.hla_end - int(1e6) < max_t, f"{sample} / {query} / {max_t}"
                offset_end = params.hla_end - max_t
                cut_begin = min_q + offset_start
                cut_end = max_q + offset_end
                cut_length = cut_end - cut_begin
                assert cut_length > int(5.5e6)
                assert cut_length < int(6.5e6)
                cut_table.append(
                    (sample, query, min_q, max_q, min_t, max_t,
                     min_q + offset_start, max_q + offset_end,
                     cut_length, source_asm_unit, source_file))

        cut_table = pd.DataFrame.from_records(
            cut_table, columns=[
                "sample", "query", "query_start", "query_end",
                "target_start", "target_end",
                "cut_query_begin", "cut_query_end",
                "cut_length", "source_assembly", "source_file"
            ]
        )
        cut_table.to_csv(output.cut_table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule extract_hla_sequences:
    input:
        fasta = expand(
                WORKDIR_EVAL.joinpath(
                    "results/assemblies", "{sample}",
                    "{sample}.asm-{asm_unit}.fasta.gz"
                ),
                sample=SAMPLES,
                asm_unit=["hap1", "hap2", "unassigned"]
            ),
        cut_table = rules.build_hla_cut_table.output.cut_table
    output:
        table_comp = DIR_RES.joinpath(
            "regions", "hla", "assembly_seq_composition.tsv"
        ),
        all_seqs = DIR_RES.joinpath(
            "regions", "hla", "assembly_all_hla.fasta.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("extract_region.py"),
        out_dir=lambda wildcards, output: pathlib.Path(output.table_comp).parent.joinpath("single_fasta")
    shell:
        "{params.script} --fasta {input.fasta} --cut-table {input.cut_table} "
            "--add-suffix HLA --dump-separately {params.out_dir} --clean-out-dir "
            "--dump-merged stdout --dump-stats {output.table_comp} "
            " | "
        "bgzip > {output.all_seqs}"
            " && "
        "samtools faidx {output.all_seqs}"


rule extract_hla_reference_sequence:
    input:
        fasta = WORKDIR_EVAL.joinpath(
            "global_ref", "chm13v2.0.fasta"
        ),
        faidx = WORKDIR_EVAL.joinpath(
            "global_ref", "chm13v2.0.fasta.fai"
        )
    output:
        fasta = DIR_LOCAL_REF.joinpath(
            "chm13v2.0.HLA.fasta.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        hla_start = int(28e6),
        hla_end = int(34e6)
    shell:
        "samtools faidx {input.fasta} chr6:{params.hla_start}-{params.hla_end} "
        " | gzip > {output.fasta}"


rule produce_hla_msa:
    input:
        ref_fasta = rules.extract_hla_reference_sequence.output.fasta,
        asm_fasta = rules.extract_hla_sequences.output.all_seqs,
    output:
        msa = DIR_RES.joinpath(
            "regions", "hla", "assembly_all_hla.msa.fasta"
        ),
        distmat = DIR_RES.joinpath(
            "regions", "hla", "assembly_all_hla.msa-distmat.txt"
        )
    log:
        DIR_LOG.joinpath("regions", "hla", "extract", "clustalo.log")
    conda:
        DIR_ENVS.joinpath("msa.yaml")
    threads: 12
    resources:
        mem_mb = lambda wildcards, attempt: 24576 * attempt,
        time_hrs = lambda wildcards, attempt: 47 * attempt
    shell:
        "zcat {input} | clustalo --in - --seqtype=DNA --infmt=fasta "
        "--percent-id --outfmt=fasta --threads 6 --verbose "
        "--full --distmat-out={output.distmat} "
        "--out {output.msa} &> {log}"
