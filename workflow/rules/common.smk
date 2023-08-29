import pandas as pd

import os
import sys

# --------  Load sample sheet -------- #
df = pd.read_table(
    config["manifest"], dtype=str, index_col="sample"
).fillna("N/A")  # sample\tquery_fofn\treference_fofn\tcomparison_type


# --------  Constraints -------- #
wildcard_constraints:
    input_type="query|reference",


# --------  Input functions -------- #
def get_final_output(wildcards):
    final_outputs = [
        "results/read_qv/{sample}/query-reference_kqv.txt.gz",
        "results/plots/{sample}/kde-before_filter.png",
        "results/plots/{sample}/kde-after_filter.png",
        "results/reads_filtered/{sample}/filtered_out/kraken2/summary.tsv.gz",
    ]

    if config["new_fastq"]:
        final_outputs.append("results/reads_filtered/{sample}/fastq.fofn")

    return expand(final_outputs, sample=df.index)

def get_kraken2_summaries(wildcards):
    fofn_df = get_query_fastq(sample_name=wildcards.sample)
    return expand(
        "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-summary.tsv.gz",
        sample=wildcards.sample,
        cell_name=fofn_df.index.tolist()
    )

def get_reads(which_one="fastq_files"):
    def inner(wildcards):
        if which_one == "fastq_files":
            fofn_path = df.at[wildcards.sample, f"{wildcards.input_type}_fofn"]
            return pd.read_table(fofn_path, names=["filepath"]).filepath.tolist()
        elif which_one == "cell":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            return fofn_df.at[wildcards.cell_name, "filepath"]
        elif which_one == "fai":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            fofn_df["fai"] = fofn_df["filepath"] + ".fai"
            return fofn_df.at[wildcards.cell_name, "fai"]

    return inner

def get_comparison_type(wildcards):
    ct = df.at[wildcards.sample, "comparison_type"]
    if ct == "N/A":
        return "self"
    else:
        return ct

def get_query_fastq(sample_name):
    fofn_path = df.at[sample_name, "query_fofn"]
    fofn_df = pd.read_table(fofn_path, names=["filepath"])

    # Get the cell name minus the extension
    def get_cell_name(fn):
        bn = os.path.basename(fn).split(".")[:-1]
        return ".".join(bn)

    fofn_df["cell_name"] = fofn_df.filepath.map(get_cell_name)

    fofn_df.set_index("cell_name", inplace=True)
    return fofn_df


def get_query_outs(which_one):
    def inner(wildcards):
        if which_one == "qv":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            cell_specific_qv = fofn_df.apply(
                lambda row: f"results/read_qv/{wildcards.sample}/{row.name}-reference_qv.txt.gz",
                axis=1,
            ).tolist()

            return cell_specific_qv
        elif which_one == "new_fastqz":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            new_fastqz = fofn_df.apply(
                lambda row: f"results/reads_filtered/{wildcards.sample}/fastq/{row.name}-subset.fastq.gz",
                axis=1,
            ).tolist()

            return new_fastqz
        elif which_one == "filtered_yak":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            filtered_yak = fofn_df.apply(
                lambda row: f"results/reads_filtered/{wildcards.sample}/{row.name}-reference_qv-filtered.txt.gz",
                axis=1,
            ).tolist()

            return filtered_yak
        else:
            raise ValueError(
                f"Unsupported param in get_query_outs (current: {which_one}) at the moment. Choose from (qv, new_fastqz, filtered_yak)"
            )

    return inner
