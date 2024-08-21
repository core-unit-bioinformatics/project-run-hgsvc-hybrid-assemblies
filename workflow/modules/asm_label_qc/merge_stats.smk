
rule compute_merged_label_count_statistics:
    input:
        regions = rules.add_ngap_sizes.output.table
    output:
        dump = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.augmented.count-stats.pck"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import collections as col
        import pickle as pck

        df = pd.read_csv(input.regions, sep="\t", header=0)

        combo_counter = col.Counter()
        combo_counter["total_regions_num"] = df.shape[0]
        combo_counter["total_regions_bp"] = (df["end"] - df["start"]).sum()

        assembly_length = df.drop_duplicates("seq", inplace=False)["seq_length"].sum()
        ngap_length = df.drop_duplicates("seq", inplace=False)["ngap_length"].sum()
        combo_counter["total_assembly_length"] = assembly_length
        combo_counter["total_ngap_length"] = ngap_length
        combo_counter["total_gapless_length"] = assembly_length - ngap_length

        def count_labels(row, counter):

            if row.labels == "no-labels":
                # skip - no annotation
                return
            region_length = row.end - row.start
            total_labels = row.labels.split(",")
            distinct_labels = sorted(set(l.split("::")[0] for l in total_labels))
            num_distinct = len(distinct_labels)

            recorded_labels = set()
            # in the following, note that total_labels can contain
            # the same label several times if it was merged several
            # times into a larger region; hence, processing state
            # for each label is tracked via recorded_labels
            for label_bp in total_labels:
                label, bp = label_bp.split("::")
                bases = int(bp)

                # total stats
                if label not in recorded_labels:
                    counter[(label, "unpaired", "total", "region_num")] += 1
                    counter[(label, "unpaired", "total", "region_bp")] += region_length
                counter[(label, "unpaired", "total", "label_bp")] += bases

                if num_distinct == 1:
                    # singleton stats
                    counter[(label, "unpaired", "single", "region_num")] += 1
                    counter[(label, "unpaired", "single", "region_bp")] += region_length
                    counter[(label, "unpaired", "single", "label_bp")] += bases
                else:
                    # merged stats
                    if label not in recorded_labels:
                        counter[(label, "unpaired", "merged", "region_num")] += 1
                        counter[(label, "unpaired", "merged", "region_bp")] += region_length
                    counter[(label, "unpaired", "merged", "label_bp")] += bases

                recorded_labels.add(label)

            if num_distinct > 1:
                # combination stats - merged equals total
                for (a,b) in itt.combinations(distinct_labels, 2):
                    counter[("pair", (a, b), "region_num")] += 1
                    counter[("pair", (a, b), "region_bp")] += region_length
                counter[("combination", tuple(distinct_labels), "region_num")] += 1
                counter[("combination", tuple(distinct_labels), "region_bp")] += region_length

            return

        _ = df.apply(count_labels, axis=1, args=(combo_counter,))
        with open(output.dump, "wb") as cache:
            _ = pck.dump(combo_counter, cache)
    # END OF RUN BLOCK


rule compute_association_label_annotation:
    input:
        regions = rules.add_ngap_sizes.output.table
    output:
        table = DIR_RES.joinpath(
            "asm_label_qc", "merge_tables", "by-sample",
            "{sample}", "{sample}.merged-issues.{span}.augmented.assoc-tests.tsv"
        )
    conda:
        DIR_ENVS.joinpath("assessem.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_label_qc", "label_assoc.py")
    shell:
        "{params.script} --augmented-regions {input.regions} --output {output.table}"


rule run_all_merged_region_stats:
    input:
        dumps = expand(
            rules.compute_merged_label_count_statistics.output.dump,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        ),
        tables = expand(
            rules.compute_association_label_annotation.output.table,
            sample=SAMPLES,
            span=["wg", "ps", "wg-no-ont", "ps-no-ont"]
        )
