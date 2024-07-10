

EVAL_ALIGN_PRECISION = ["exact", "struct-50"]


rule align_verkko_parent_to_child:
    input:
        child = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{child}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        chl_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{child}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
        parent = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{parent}.vrk-ps-sseq.asm-mrg-tag.fasta"
        ),
        par_idx = WORKDIR_EVAL.joinpath(
            "proc/05-preprocess/merge_tag_asm",
            "{parent}.vrk-ps-sseq.asm-mrg-tag.fasta.fai"
        ),
    output:
        paf = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "raw",
            "{parent}-to-{child}.vrk-ps-sseq.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("minimap.yaml")
    threads: 12
    resources:
        mem_mb=lambda wildcards, attempt: 65536 * attempt,
        time_hrs=lambda wildcards, attempt: 4 * attempt
    shell:
        "minimap2 -t {threads} -x asm5 --eqx --MD -c --secondary=no "
        "{input.child} {input.parent} | pigz > {output.paf}"


rule normalize_trio_align_parent_to_child:
    input:
        paf = rules.align_verkko_parent_to_child.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "norm",
            "{parent}-to-{child}.vrk-ps-sseq.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=DIR_SCRIPTS.joinpath("normalize_paf.py").resolve(strict=True)
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule evaluate_normalized_trio_pafs:
    input:
        paf = DIR_PROC.joinpath(
            "asm_compare", "trio_align", "norm",
            "{parent}-to-{child}.vrk-ps-sseq.norm-paf.tsv.gz"
        ),
        gsize = select_genome_size_file
    output:
        stats = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.stats.tsv"
        ),
        aln_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aligned-regions.tsv.gz"
        ),
        sup_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.support-regions.tsv.gz"
        ),
        mrg_regions = DIR_RES.joinpath(
            "asm_compare", "trio",
            "regions",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.merged-regions.tsv.gz"
        )
    log:
        DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.aln-filter.log"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_compare", "eval_paf.py").resolve(strict=True),
        min_aln_len=lambda wildcards: to_int(wildcards.min_aln_len),
        min_seq_len=lambda wildcards: to_int(wildcards.min_seq_len),
        precision=lambda wildcards: parse_precision_param(wildcards.precision, "setting"),
        struct_t=lambda wildcards: parse_precision_param(wildcards.precision, "threshold"),
    shell:
        "{params.script} --verbose --paf-file {input.paf} "
        "--genome-size {input.gsize} "
        "--precision {params.precision} "
        "--retain-mapq-geq {wildcards.min_mapq} "
        "--retain-seqlen-geq {params.min_seq_len} "
        "--retain-alnlen-geq {params.min_aln_len} "
        "--struct-err-geq {params.struct_t} "
        "--split-seq-tag '.' "
        "--out-seq-stats {output.stats} "
        "--out-aligned-regions {output.aln_regions} "
        "--out-support-regions {output.sup_regions} "
        "--out-merged-regions {output.mrg_regions} "
        "&> {log}"


rule postprocess_merged_trio_regions:
    input:
        mrg_regions = rules.evaluate_normalized_trio_pafs.output.mrg_regions,
        gsize = select_genome_size_file
    output:
        mrg_stats = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics",
            "{parent}-to-{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.{precision}.mrg-regions-stats.tsv"
        )
    run:
        import pandas as pd

        seq_lens = dict()
        with open(input.gsize, "r") as listing:
            for line in listing:
                seq_name, length = line.split()[:2]
                seq_lens[seq_name] = int(length)
                # genome size file is fasta index by default;
                # that contains tagged sequence names
                seq_name, _ = seq_name.rsplit(".", 1)
                seq_lens[seq_name] = int(length)

        df = pd.read_csv(input.mrg_regions, sep="\t", header=0)
        df["length"] = df["end"] - df["start"]

        pct_stats = df["length"].describe(percentiles=[0.25,0.5,0.75,0.95])
        pct_stats.rename(
            {
                "count": "num_regions", "mean": "length_mean",
                "std": "stddev", "min": "length_pct_0",
                "25%": "length_pct_25", "50%": "length_pct_50", "75%": "length_pct_75",
                "95%": "length_pct_95", "max": "length_pct_100"
            },
            inplace=True
        )
        num_seqs = df["seq_name"].nunique()
        known_seqs = set(df["seq_name"].unique())
        total_length = sum(seq_lens[sn] for sn in known_seqs)
        support_length = df["length"].sum()
        support_pct = round(support_length/total_length * 100, 2)

        stats_row = [
            wildcards.parent, wildcards.child, wildcards.precision,
            num_seqs, total_length, support_length, support_pct,
            int(pct_stats["num_regions"]), int(round(pct_stats["length_mean"], 0)),
            int(round(pct_stats["length_pct_25"],0)), int(round(pct_stats["length_pct_50"],0)),
            int(round(pct_stats["length_pct_75"],0)), int(round(pct_stats["length_pct_95"],0)),
            int(round(pct_stats["length_pct_100"],0))
        ]

        header = [
            "parent", "child", "precision",
            "num_contigs", "total_length", "supported_length", "supported_pct",
            "num_support_regions", "region_length_mean", "region_length_pct25",
            "region_length_pct50", "region_length_pct75", "region_length_pct95",
            "region_length_pct100"
        ]

        with open(output.mrg_stats, "w") as dump:
            _ = dump.write("\t".join(header) + "\n")
            _ = dump.write("\t".join(map(str, stats_row)) + "\n")
    # END OF RUN BLOCK


TRIO_MAP = {
    "HG00514": ["HG00512", "HG00513"],
    "HG00733": ["HG00731", "HG00732"],
    "NA19240": ["NA19238", "NA19239"]
}


rule concat_merged_trio_region_stats:
    input:
        tables = lambda wildcards: expand(
            rules.postprocess_merged_trio_regions.output.mrg_stats,
            parent=TRIO_MAP[wildcards.child],
            precision=EVAL_ALIGN_PRECISION,
            allow_missing=True
        )
    output:
        summary = DIR_RES.joinpath(
            "asm_compare", "trio",
            "statistics", "summary",
            "{child}.vrk-ps-sseq.mapq-{min_mapq}.seq-{min_seq_len}.aln-{min_aln_len}.merged-regions.tsv"
        )
    run:
        import pandas as pd
        import pathlib as pl

        concat = []
        for table_file in input.tables:
            df = pd.read_csv(table_file, header=0, sep="\t")
            concat.append(df)
        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.sort_values(["parent", "precision"], inplace=True)

        concat.to_csv(output.summary, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_asm_trio_align:
    input:
        summary = expand(
            rules.concat_merged_trio_region_stats.output.summary,
            child=["HG00733", "HG00514", "NA19240"],
            min_mapq=[1],
            min_seq_len=["0", "100k"],
            min_aln_len=["0", "10k"]
        )
