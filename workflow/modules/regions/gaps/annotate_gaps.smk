
rule segdups_to_bed:
    input:
        sd95 = rules.split_segdup_annotation.output.sd95,
        sd98 = rules.split_segdup_annotation.output.sd98
    output:
        bed_95 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-095.bed"
        )),
        bed_98 = temp(DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-098.bed"
        ))
    run:
        import pandas as pd
        import hashlib as hl

        def compute_row_id(row):

            row_hash = hl.md5(f"{row.seq}{row.start}{row.end}".encode("utf-8")).hexdigest()
            row_name = f"{row_hash}|{row.label}"

        for input_file, output_file, region_label in zip(input, output, ["SD95", "SD98"]):
            df = pd.read_csv(input_file, sep="\t", header=0)
            df["label"] = region_label
            df["name"] = df.apply(compute_row_id, axis=1)
            df.rename({"seq": "#contig"}, axis=1, inplace=True)
            df.sort_values(["#contig", "start", "end"], inplace=True)
            df.to_csv(output_file, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_overlapping_segdups:
    input:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.bed"
        ),
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "label_segdup",
            "{sample}.sd-{pct_id}.mrg.bed"
        ),
    wildcard_constraints:
        pct_id = "(095|098)"
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    shell:
        "bedtools merge -c 5 -o first -i {input.bed} > {output.bed}"


rule annotate_gaps_with_segdups:
    """Annotate assemblies (assembly coordinate space)
    with segdup annotation. The input file in ref coordinates
    is just here for enforcing.
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        segdups = rules.merge_overlapping_segdups.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-sd{pct_id}.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.qry_view} -b {input.segdups} > {output.isect}"


localrules: hprc_gaps_to_bed
rule hprc_gaps_to_bed:
    input:
        tsv = HPRC_COMMON_GAPS
    output:
        bed = DIR_PROC.joinpath(
            "regions", "gaps", "norm_tables",
            "hprc_common_gaps.bed"
        )
    run:
        import pandas as pd
        df = pd.read_csv(input.tsv, sep="\t", header=0)
        df.sort_values(["chrom", "start", "end"], inplace=True)
        df.rename({"chrom": "#chrom"}, axis=1, inplace=True)
        df.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule annotate_gaps_with_hprc_gaps:
    """NB: HPRC gaps only exist in T2T space
    """
    input:
        qry_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = WORKDIR_EVAL.joinpath(
            "results", "regions", "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
        gaps = rules.hprc_gaps_to_bed.output.bed
    output:
        isect = DIR_PROC.joinpath(
            "regions", "gaps", "intersections",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.isect-hprcgaps.tsv"
        )
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wao -a {input.trg_view} -b {input.gaps} > {output.isect}"


rule simplify_hprc_gap_intersection:
    input:
        table = rules.annotate_gaps_with_hprc_gaps.output.isect
    output:
        tsv = DIR_RES.joinpath(
            "regions", "gaps", "norm_tables",
            "{sample}",
            "{sample}.asm-{asm_unit}.{refgenome}.ref-coord.hprc-gaps.tsv"
        )
    run:
        import pandas as pd
        header = ["aln_seq", "aln_start", "aln_end", "aln_label", "aln_length", "aln_info"]
        header += ["aln_base_block", "aln_coarse_block"]
        header += ["gap_chrom", "gap_start", "gap_end", "gap_id", "hprc_haps"]
        header += ["gap_sd_assoc", "gap_cov_drop", "gap_win_start", "gap_win_end", "overlap_bp"]

        df = pd.read_csv(input.table, sep="\t", header=None, names=header)
        df = df.loc[df["overlap_bp"] > 0, :].copy()
        df.drop(["hprc_haps", "gap_win_start", "gap_win_end", "aln_length", "aln_base_block", "aln_coarse_block"], axis=1, inplace=True)
        df["sample"] = wildcards.sample
        df["asm_unit"] = wildcards.asm_unit
        df["gap_length"] = df["gap_end"] - df["gap_start"]
        df["overlap_pct"] = (df["gap_length"] / df["overlap_bp"] * 100).round(2)
        df.sort_values(["aln_seq", "aln_start", "aln_end"], inplace=True)
        df.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OR RUN BLOCK


rule run_all_annotate_gaps:
    input:
        tables = expand(
            rules.annotate_gaps_with_segdups.output.isect,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["hg38", "t2tv2"],
            pct_id=["095", "098"]
        ),
        gaps = expand(
            rules.simplify_hprc_gap_intersection.output.tsv,
            sample=SAMPLES,
            asm_unit=MAIN_ASSEMBLY_UNITS,
            refgenome=["t2tv2"]
        )
