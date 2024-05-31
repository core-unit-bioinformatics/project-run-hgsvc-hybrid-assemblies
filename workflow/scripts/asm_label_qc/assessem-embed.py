#!/usr/bin/env python3

import argparse as argp
import collections as col
import pathlib as pl
import sys
import timeit

import pandas as pd
import pacmap
import numpy as np
import scipy.stats as stats


def parse_command_line():

    parser = argp.ArgumentParser()

    parser.add_argument(
        "--data-cache", "-dc",
        type=lambda x: pl.Path(x).resolve(strict=True),
        required=True,
        dest="data_cache",
        help="Path to existing data cache file (HDF format)."
    )

    parser.add_argument(
        "--bin-size", "-bs",
        type=int,
        default=int(1e4),
        dest="bin_size",
        help="Bin size in bp to aggregate cached data. Default: 10 000"
    )

    parser.add_argument(
        "--binned-data", "-ob", "--out-binned",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="binned_data",
        help="Path to output file holding the binned data."

    )

    parser.add_argument(
        "--transformed-data", "-ot", "--out-trans",
        type=lambda x: pl.Path(x).resolve(strict=False),
        dest="transformed_data",
        help="Path to output file holding the transformed data."
    )

    parser.add_argument(
        "--verbose", "-vb",
        action="store_true",
        default=False,
        dest="verbose",
        help="Print logging and runtime information. Default: False"
    )

    args = parser.parse_args()


    return args


def print_runtime_summary(timer):

    for timed_op, timings in timer.items():
        delta = round(timings[3], 5)
        msg = (
            "\n==================\n"
            f"{timed_op}\n"
            f"{timings[0]} walltime: ~{delta} sec.\n"
            "==================\n"
        )
        sys.stdout.write(msg)
    return


def load_cached_genome_infos(data_cache):

    with pd.HDFStore(data_cache, "r") as hdf:
        try:
            genome_sizes = hdf["genome_size"]
        except KeyError:
            raise RuntimeError("Invalid data cache file - no entry 'genome_size'")
    total_length = genome_sizes.loc[genome_sizes["sequence"] == "genome", "size"].iloc[0]
    offsets = dict(
        (row.sequence, row.start) for row in genome_sizes.itertuples()
    )
    return total_length, offsets


def load_cached_track_infos(data_cache):

    with pd.HDFStore(data_cache, "r") as hdf:
        try:
            track_infos = hdf["data_tracks"]
            track_labels = set(track_infos.index)
        except KeyError:
            track_infos = None
            track_labels = set()
    return track_infos, track_labels


def load_region_track(data_cache, track_label, num_regions):

    with pd.HDFStore(data_cache, "r") as hdf:
        regions = hdf[f"tracks/{track_label}"]
    if not num_regions == regions.shape[0]:
        err_msg = (
            f"Size mismatch for number of regions: {track_label}\n"
            f"In metadata: {num_regions}\n"
            f"Loaded from cache: {regions.shape[0]}\n"
        )
        raise ValueError(err_msg)
    return regions


def load_cached_track_data(data_cache, track_info, genome_size, bin_size, blunt_end):

    data_type = {"int8": np.int8, "uint8": np.uint8}[track_info.data_range]

    regions = load_region_track(data_cache, track_info.Index, track_info.num_regions)

    data = np.zeros(genome_size, dtype=data_type)
    for row in regions.itertuples():
        data[row.start:row.end] = row.score
    data = data[:blunt_end]

    data = data.reshape((-1, bin_size))
    score_type = track_info.score_type
    if score_type == "binary":
        data = (np.count_nonzero(data, axis=1)/bin_size).round(3).astype(np.float16)
    elif score_type == "categorical":
        data = np.apply_along_axis(lambda row: stats.mode(row).mode, 1, data)
    elif score_type == "quantitative":
        data = np.median(data, axis=1)
    else:
        raise ValueError(f"Unhandled score type: {track_info.Index}")

    # the aggregate operations will result in floats
    # but high precision is pointless here
    data = data.astype(np.float16, casting="same_kind", copy=False)
    return data


def build_dataset(data_cache, bin_size, timer):

    genome_size, _ = load_cached_genome_infos(data_cache)
    num_bins = genome_size // bin_size
    blunt_end = num_bins * bin_size

    track_infos, _ = load_cached_track_infos(data_cache)
    num_tracks = track_infos.shape[0]

    num_tracks = len(track_infos)
    timer["alloc"] = [
        f"Allocating memory: {num_tracks}X{num_bins}",
        timeit.default_timer(),
        None, None
    ]
    dataset = np.zeros((num_tracks, num_bins), dtype=np.float16)
    timer["alloc"][2] = timeit.default_timer()
    timer["alloc"][3] = timer["alloc"][2] - timer["alloc"][1]

    label_order = []
    for row_num, row in enumerate(track_infos.itertuples(), start=0):
        track_label = row.Index
        timer[f"load:{track_label}"] = [
            f"Loading track data: {row.source_file}",
            timeit.default_timer(),
            None, None
        ]
        dataset[row_num, :] = load_cached_track_data(data_cache, row, genome_size, bin_size, blunt_end)
        label_order.append(track_label)
        timer[f"load:{track_label}"][2] = timeit.default_timer()
        timer[f"load:{track_label}"][3] = timer[f"load:{track_label}"][2] - timer[f"load:{track_label}"][1]

    dataset = dataset.transpose()
    dataset = pd.DataFrame(dataset, columns=label_order)
    return dataset


def compute_embedding_transformation(dataset):

    assert dataset.shape[0] > dataset.shape[1]

    # FOLLOWS default parameterization from quickstart
    # initializing the pacmap instance
    # Setting n_neighbors to "None" leads to a default choice shown below in "parameter" section
    embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0)
    # fit the data (The index of transformed data corresponds to the index of the original data)
    transformed = embedding.fit_transform(dataset, init="pca")
    return transformed


def main():

    timer = col.OrderedDict()
    timer["assessem-embed"] = ["Total runtime", timeit.default_timer(), None, None]
    args = parse_command_line()

    dataset = build_dataset(args.data_cache, args.bin_size, timer)
    args.binned_data.parent.mkdir(exist_ok=True, parents=True)
    #dataset.to_feather(args.binned_data, compression='lz4')
    dataset.to_csv(args.binned_data, sep="\t", header=True, index=False)

    timer["embed"] = ["PaCMAP embedding", timeit.default_timer(), None, None]
    transformed = compute_embedding_transformation(dataset)
    timer["embed"][2] = timeit.default_timer()
    timer["embed"][3] = timer["embed"][2] - timer["embed"][1]

    args.transformed_data.parent.mkdir(exist_ok=True, parents=True)
    np.save(args.transformed_data, transformed, allow_pickle=False)

    timer["assessem-embed"][2] = timeit.default_timer()
    timer["assessem-embed"][3] = timer["assessem-embed"][2] - timer["assessem-embed"][1]

    if args.verbose:
        print_runtime_summary(timer)

    return 0


if __name__ == "__main__":
    main()
