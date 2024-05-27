
def collect_merqury_kmer_tracks(wildcards):

    sample_id, assm = wildcards.sample.split(".")
    assembler = {"vrk-ps-sseq": "verkko", "hsm-ps-sseq": "hifiasm"}[assm]
    kmer_tracks = MERQURY_RESULT_ROOT_ALL.joinpath(
        f"{sample_id}-{assembler}"
    ).glob("**/*only.bed")

    files = []
    for kmer_track in kmer_tracks:
        assert assembler in kmer_track.name
        files.append(kmer_track)
    if assembler == "verkko":
        if not len(files) == 3:
            raise RuntimeError(f"Merqury verkko error: {wildcards}")
    if assembler == "hifiasm":
        if not len(files) == 2:
            raise RuntimeError(f"Merqury hifiasm error: {wildcards}")
    return sorted(files)
