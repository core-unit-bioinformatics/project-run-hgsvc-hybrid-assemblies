import pathlib

def collect_merqury_kmer_tracks(wildcards):

    sample_id, assm = wildcards.sample.split(".")
    assembler = {"vrk-ps-sseq": "verkko", "hsm-ps-sseq": "hifiasm"}[assm]

    kmer_tracks = sorted(MERQURY_RESULT_ROOT_ALL.joinpath(
        f"{sample_id}-{assembler}"
        ).glob("**/*only.bed")
    )
    if assembler == "verkko":
        kmer_tracks += sorted(MERQURY_RESULT_ROOT_ALL.joinpath(
            f"{sample_id}-{assembler}-unassigned"
            ).glob("**/*only.bed")
        )

    files = []
    for kmer_track in kmer_tracks:
        assert assembler in kmer_track.name or assm in kmer_track.name, [p.name for p in kmer_tracks]
        files.append(kmer_track)
    if assembler == "verkko":
        if not len(files) == 3:
            raise RuntimeError(f"Merqury verkko error: {wildcards}")
    if assembler == "hifiasm":
        if not len(files) == 2:
            raise RuntimeError(f"Merqury hifiasm error: {wildcards}")
    return sorted(files)


def get_assessem_cli_parameters(input_files, get_list):

    track_files = []
    track_labels = []
    score_columns = []
    out_list = None
    for input_file in input_files:
        file_name = pathlib.Path(input_file).name
        if file_name.endswith(".sizes.txt"):
            continue
        track_files.append(input_file)
        if file_name.endswith(".flagger-labels.tsv.gz"):
            track_labels.append("flagger")
            score_columns.append("score")
        elif file_name.endswith(".hifi.inspector-errors.tsv.gz"):
            track_labels.append("inspect_hifi")
            score_columns.append("binary")
        elif file_name.endswith(".ont.inspector-errors.tsv.gz"):
            track_labels.append("inspect_ont")
            score_columns.append("binary")
        elif file_name.endswith(".merqury-asmonly-kmer.tsv.gz"):
            track_labels.append("merqury")
            score_columns.append("binary")
        elif file_name.endswith(".mosdepth-windowed.tsv.gz"):
            # NA18989.vrk-ps-sseq.hifi.mq00.mosdepth-windowed.tsv.gz
            parts = file_name.split(".")
            read_type = parts[-5]
            assert read_type in ["hifi", "ont"]
            mapq = parts[-4]
            track_labels.append(f"rd_{read_type}_{mapq}")
            score_columns.append("cov")
        elif file_name.endswith(".nucfreq-flagged.bed.gz"):
            track_labels.append("nucfreq")
            score_columns.append("num_hets")
        else:
            raise ValueError(f"Cannot process filename: {file_name}")
    if get_list == "files":
        out_list = track_files
    if get_list == "labels":
        out_list = track_labels
    if get_list == "columns":
        out_list = score_columns
    assert out_list is not None
    return out_list



