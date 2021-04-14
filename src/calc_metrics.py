import argparse
import numpy as np
import os
import pandas as pd
import re
import subprocess
import sys

from io import StringIO
from pybedtools import BedTool
from scipy import spatial  # cosine similarity function
from ubiquerg import is_command_callable

SEGMENTATION_LENGTH = -1

if not is_command_callable("ailist"):
    print("Command 'ailist' not callable. Install ailist to continue.")
    sys.exit(1)


def run_ailist(f1, f2, database_output=True):
    """
    f:          the query file
    num_files:  the number of files used in the segmentation.
                this will be used to find the segmentation file in regionMap
    database_output: if the ailist results are the database regions (True) or the query regions (False)

    returns pandas dataframe of the database regions that overlapped the query regions

    you need ailist in your path
    """

    query = f1
    database = f2

    if not database_output:
        ailist_run = subprocess.run(
            ["ailist", database, query], encoding="utf-8", stdout=subprocess.PIPE
        ).stdout
        ailist_result = ailist_run.rsplit("\n", 3)[0].split("\n", 2)[2]
    else:
        ailist_run = subprocess.run(
            ["ailist", query, database], encoding="utf-8", stdout=subprocess.PIPE
        ).stdout
        ailist_result = ailist_run.rsplit("\n", 3)[0].split("\n", 2)[2]
    df = pd.read_csv(StringIO(ailist_result), sep="\t", header=None)
    df[0] = df[0].apply(lambda x: x.rstrip(":"))
    df[1].astype(int)
    df[2].astype(int)
    df[3].astype(int)
    return df


def ailist_vectorize(f, segmentation_file):
    """
    creates a vector from the ailist results
    """

    ailist_df = run_ailist(f, segmentation_file)
    ailist_df.columns = ["chrom", "start", "end", "overlaps"]
    ailist_df.sort_values(["chrom", "start", "end"], inplace=True)

    segmentation_df = pd.read_csv(
        segmentation_file,
        sep="\t",
        header=0,
        usecols=[0, 1, 2],
        dtype={0: str, 1: int, 2: int},
    )
    if not str(segmentation_df.iloc[0, 1]).isdigit():
        segmentation_df.columns = segmentation_df.iloc[0]
        segmentation_df = segmentation_df[1:]
    segmentation_df.columns = ["chrom", "start", "end"]
    segmentation_df["start"].astype(int)
    segmentation_df["end"].astype(int)
    vector_df = segmentation_df.merge(
        ailist_df, how="left", on=["chrom", "start", "end"]
    )
    # print("number of overlaps: ", vector_df[vector_df['overlaps'].notnull()].shape)
    # print("percentage of overlaps: ", ailist_df.shape[0]/segmentation_df.shape[0])
    vector_df["overlaps"] = vector_df["overlaps"].apply(lambda x: 1 if x >= 1 else 0)
    return list(vector_df["overlaps"])


def cosine_similarity(v1, v2):
    return 1 - spatial.distance.cosine(v1, v2)


def euclidean_dist(v1, v2):
    # with a vector of length L that can only have values 0 to 1, sqrt(L) is the max distance
    return 1 - np.linalg.norm(np.array(v1) - np.array(v2)) / SEGMENTATION_LENGTH


def jaccard(f1, f2):
    num_regions_f1 = 0
    with open(f1) as f:
        for i, l in enumerate(f):
            pass
        num_regions_f1 = i + 1
    num_regions_f2 = 0
    with open(f2) as f:
        for i, l in enumerate(f):
            pass
        num_regions_f2 = i + 1
    if num_regions_f2 > num_regions_f1:
        ailist_df = run_ailist(f2, f1)
    else:
        ailist_df = run_ailist(f1, f2)
    overlaps = ailist_df[ailist_df[3] != 0].shape[0]
    norm = overlaps / (num_regions_f1 + num_regions_f2 - overlaps)
    return norm


def coverage(f1, f2):
    a = BedTool(f1)
    b = BedTool(f2)
    res1 = a.coverage(b)
    res2 = b.coverage(a)
    df1 = pd.read_csv(res1.fn, sep="\t", header=None)
    df2 = pd.read_csv(res2.fn, sep="\t", header=None)
    return (
        df1[len(df1.columns) - 1].mean()
        + df2[len(df2.columns) - 1].mean()
    ) / 2.0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-f", help="original file", required=True)
    # parser.add_argument("-s", help="segmentation/universe file", required=True)
    # parser.add_argument("-b", help="bedshifted files directory", required=True)
    parser.add_argument("-i", help="Input directory with perturbed BED files", required=True)
    parser.add_argument("-o", help="Output directory for scores", required=True)
    parser.add_argument("-p", help="PEP config file", required=True)
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="verbose", default=False
    )

    args = parser.parse_args()

    verbose = args.verbose
    pep_config = args.p
    output_dir = args.o
    # bedshifted_dir = args.b
    # original = args.f
    # segmentation = args.s

    # For interactive testing
    # verbose = True
    # pep_config = "project_config.yaml"

    import peppy
    import glob

    prj = peppy.Project(pep_config)

    # Set up arrays to store scores
    jaccard_distances = {s.sample_name: [] for s in prj.samples}
    coverage_distances = {s.sample_name: [] for s in prj.samples}
    euc_distances = {s.sample_name: [] for s in prj.samples}
    cos_distances = {s.sample_name: [] for s in prj.samples}

    for s in prj.samples:
        print(f"Processing sample {s.sample_name}")
        # With this new approach we compute these segmentation and original file paths
        # for each sample
        original_file_path = os.path.join(s.base_path, s.file)
        segmentation_file_path = s.universe
        with open(segmentation_file_path) as f:
            for i, l in enumerate(f):
                pass
            SEGMENTATION_LENGTH = np.sqrt(i + 1)
        original_vector = ailist_vectorize(original_file_path, segmentation_file_path)
        fp = os.path.join(args.i, s.sample_name, s.sample_name + "_rep*.bed")

        rep_files = sorted(glob.glob(fp))
        if len(rep_files) < 1:
            print("Error: no replicate files found")

        for f in rep_files:
            if verbose:
                print("Processing ", f)

            ## JACCARD
            jaccard_score = jaccard(original_file_path, f)
            jaccard_distances[s.sample_name].append(jaccard_score)
            if verbose:
                print("Jaccard:\t", jaccard_score)

            ## COVERAGE
            coverage_score = coverage(original_file_path, f)
            if verbose:
                print("Coverage:\t", coverage_score)
            coverage_distances[s.sample_name].append(coverage_score)

            ## EUCLIDEAN
            perturbed_vector = ailist_vectorize(f, segmentation_file_path)
            euclidean_score = euclidean_dist(original_vector, perturbed_vector)
            if verbose:
                print("Euclidean:\t", euclidean_score)
            euc_distances[s.sample_name].append(euclidean_score)

            ## COSINE
            cosine_score = cosine_similarity(original_vector, perturbed_vector)
            if verbose:
                print("Cosine:\t\t", cosine_score)
            cos_distances[s.sample_name].append(cosine_score)

    if verbose:
        print("Number of elements calculated:")
        print("Jaccard: ", [(k, len(v)) for k, v in jaccard_distances.items()])
        print("Coverage: ", [(k, len(v)) for k, v in coverage_distances.items()])
        print("Euclidean: ", [(k, len(v)) for k, v in euc_distances.items()])
        print("Cosine: ", [(k, len(v)) for k, v in cos_distances.items()])

    df1 = pd.DataFrame.from_dict(jaccard_distances)
    df2 = pd.DataFrame.from_dict(coverage_distances)
    df3 = pd.DataFrame.from_dict(euc_distances)
    df4 = pd.DataFrame.from_dict(cos_distances)
    df1.to_csv(os.path.join(output_dir, "jaccard_results.txt"), index=False)
    df2.to_csv(os.path.join(output_dir, "coverage_results.txt"), index=False)
    df3.to_csv(os.path.join(output_dir, "euclidean_results.txt"), index=False)
    df4.to_csv(os.path.join(output_dir, "cosine_results.txt"), index=False)
