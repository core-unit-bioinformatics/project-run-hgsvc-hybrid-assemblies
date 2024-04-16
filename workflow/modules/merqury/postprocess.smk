

localrules: normalize_merqury_qv_estimates
rule normalize_merqury_qv_estimates:
    input:
        summary_file = lambda wildcards: find_merqury_output_file(wildcards.sample, "qv_summary"),
        detail_file = lambda wildcards: find_merqury_output_file(wildcards.sample, "qv_detail")
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-qv-est.tsv"
        )
    run:
        import pandas as pd

        processed_files = []

        summary = pd.read_csv(
            input.summary_file[0], sep="\t", header=None,
            names=["entity", "error_bp", "total_bp", "qv_est", "err_rate"]
        )
        processed_files.append(shorten_merqury_file_path(input.summary_file[0]))
        assert summary.shape[0] == 3
        summary["sample"] = wildcards.sample
        summary["asm_unit"] = "na"
        summary["sequence"] = "total"
        au = []
        for row in summary.itertuples():
            if "hap1" in row.entity or "h1" in row.entity:
                au.append("hap1")
            elif "hap2" in row.entity or "h2" in row.entity:
                au.append("hap2")
            elif "both" in row.entity.lower():
                au.append("wg")
            else:
                raise ValueError(row)
        summary["asm_unit"] = au

        for hap_file in input.detail_file:
            asm_unit = None
            if "hap1" in hap_file.name:
                asm_unit = "hap1"
            elif "hap2" in hap_file.name:
                asm_unit = "hap2"
            else:
                raise ValueError(hap_file)
            detail = pd.read_csv(
                hap_file, sep="\t", header=None,
                names=["sequence", "error_bp", "total_bp", "qv_est", "err_rate"]
            )
            processed_files.append(shorten_merqury_file_path(hap_file))
            detail["sample"] = wildcards.sample
            detail["asm_unit"] = asm_unit
            summary = pd.concat([summary, detail], axis=0, ignore_index=False)

        summary = summary[["sample", "asm_unit", "sequence", "error_bp", "total_bp", "qv_est", "err_rate"]]
        with open(output.tsv, "w") as table:
            for f in processed_files:
                table.write(f"# {f}\n")
            summary.to_csv(table, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_merqury_normalized_qv_estimates:
    input:
        tables = expand(
            rules.normalize_merqury_qv_estimates.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}.qv-est.tsv"
        )
    run:
        import pandas as pd
        import io

        merge = []
        for table_file in input.tables:
            table_buffer = io.StringIO()
            with open(table, "r") as table:
                for line in table:
                    if line.startswith("#"):
                        continue
                    table_buffer.write(line)
            table_buffer.seek(0)
            df = pd.read_csv(table_buffer, sep="\t", header=0)
            merge.append(df)
        merge = pd.concat(merge, axis=0, ignore_index=False)
        merge.sort_values(["sample", "asm_unit", "sequence"], inplace=True)
        merge.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


localrules: normalize_merqury_kmer_completeness
rule normalize_merqury_kmer_completeness:
    input:
        table = lambda wildcards: find_merqury_output_file(wildcards.sample, "completeness")
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}", "{sample}",
            "{sample}.merqury-kmer-completeness.tsv"
        )
    run:
        import pandas as pd

        rows = []
        with open(input.table[0], "r") as listing:
            for line in listing:
                columns = line.strip().split()
                assert columns[1] == "all"  # no clue ...
                pct_complete = float(columns[-1])
                kmer_total = int(columns[-2])
                kmer_found = int(columns[-3])
                asm_unit = None
                if "hap1" in columns[0]:
                    asm_unit = "hap1"
                elif "hap2" in columns[0]:
                    asm_unit = "hap2"
                elif "both" in columns[0].lower():
                    asm_unit = "wg"
                else:
                    raise ValueError(line.strip())
                assert asm_unit is not None
                rows.append(
                    (wildcards.sample, asm_unit, kmer_found, kmer_total, pct_complete)
                )
        df = pd.DataFrame.from_records(
            rows, columns=["sample", "asm_unit", "kmer_present", "kmer_total", "kmer_present_pct"]
        )
        with open(output.tsv, "w") as dump:
            _ = dump.write(f"# {shorten_merqury_file_path(input.table[0])}\n")
            df.to_csv(dump, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule merge_merqury_normalized_kmer_completeness:
    input:
        tables = expand(
            rules.normalize_merqury_kmer_completeness.output.tsv,
            sample=SAMPLES,
            allow_missing=True
        )
    output:
        tsv = DIR_RES.joinpath(
            "merqury", "{assembler}.kmer-completeness.tsv"
        )
    run:
        import pandas as pd

        merge = []
        for table in input.tables:
            df = pd.read_csv(table, sep="\t", header=0, comment="#")
            merge.append(df)
        merge = pd.concat(merge, axis=0, ignore_index=False)
        merge.sort_values(["sample", "asm_unit"], inplace=True)
        merge.to_csv(output.tsv, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_merqury_postprocess:
    input:
        qv_est = expand(
            rules.merge_merqury_normalized_qv_estimates.output.tsv,
            assembler=[ASSEMBLER]
        ),
        kmer_completeness = expand(
            rules.merge_merqury_normalized_kmer_completeness.output.tsv,
            assembler=[ASSEMBLER]
        )
