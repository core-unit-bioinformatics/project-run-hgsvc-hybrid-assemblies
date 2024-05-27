
rule merge_merqury_kmer_tracks:
    input:
        bed_files = collect_merqury_kmer_tracks
    output:
        bed_like = temp(DIR_PROC.joinpath(
            "asm_label_qc", "norm_merqury",
            "{sample}.asmonly-kmer.bed"
        ))
    conda:
        DIR_ENVS.joinpath("bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "cat {input.bed_files} | sort -V -k1,1 -k2,2n"
            " | "
        "bedtools merge -i /dev/stdin > {output.bed_like}"


rule normalize_merqury_kmer_tracks:
    input:
        bed_like = rules.merge_merqury_kmer_tracks.output.bed_like
    output:
        bed_like = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables",
            "merqury", "{sample}.merqury-asmonly-kmer.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        df = pd.read_csv(input.bed_like, sep="\t", header=None, names=["chrom", "start", "end"])
        df["name"] = "merqury_asmonly_kmer"
        df.to_csv(output.bed_like, sep="\t", header=True, index=False)


rule run_normalize_merqury_results:
    input:
        tables = expand(
            rules.normalize_merqury_kmer_tracks.output.bed_like,
            sample=SAMPLES
        )
