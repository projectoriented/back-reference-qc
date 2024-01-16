import os
import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

# --------  Load sample sheet -------- #
df = pd.read_table(
    config["manifest"], dtype=str
).set_index(["sample"], drop=False)

df.fillna("N/A", inplace=True)

validate(df, schema="../schemas/manifest.schema.yaml")

# --------  Constraints -------- #
wildcard_constraints:
    input_type="query|reference",
    suffix="with_reference_help|from_raw"


# --------  Input functions -------- #
def get_final_output(wildcards):
    final_outputs = []

    for row in df.itertuples():
        final_outputs.append(f"results/reads_filtered/{row.sample}/filtered_out/kraken2/summary.tsv.gz"),
        if row.reference_fofn != "N/A":
            final_outputs.append(f"results/read_qv/{row.sample}/query-reference_kqv.txt.gz"),
            final_outputs.append(f"results/plots/{row.sample}/kde-before_filter.png"),
            final_outputs.append(f"results/plots/{row.sample}/kde-after_filter.png"),

        if config["new_fastq"]:
            final_outputs.append(f"results/reads_filtered/{row.sample}/fastq.fofn")

    return expand(final_outputs)

def get_kraken2_summaries(wildcards):
    fofn_df = get_query_fastq(sample_name=wildcards.sample)
    return expand(
        "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-summary.tsv.gz",
        sample=wildcards.sample,
        cell_name=fofn_df.index.tolist()
    )

def calc_mem_gb(wildcards, input, attempt, threads):
    mb = max(1.5 * input.size_mb, 1000)
    gb = int(mb / 1000)

    if threads != 1:
        gb = int(max(gb / threads, 2))

    return gb * attempt

def get_kraken2_inputs(which_one):
    def inner(wildcards):
        kraken_param = ""
        if df.at[wildcards.sample, "reference_fofn"] == "N/A":
            fofn_df = get_query_fastq(sample_name=wildcards.sample)
            seq = fofn_df.at[wildcards.cell_name, "filepath"]
            if seq.endswith(".gz"):
                kraken_param = "--gzip-compressed"
        else:
            seq = "results/reads_filtered/{sample}/filtered_out/fasta/{cell_name}_undesirable-quality.fa"

        if which_one == "kraken_param":
            return kraken_param
        elif which_one == "seq":
            return seq
        else:
            raise ValueError(f"Unsupported arg: {which_one} for get_kraken2_inputs")
    return inner

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

    # Check for duplicates
    if fofn_df.cell_name.duplicated().any():
        # Add this column to offset duplicated namings
        fofn_df["idx"] = range(fofn_df.shape[0])

        fofn_df["cell_name"] = fofn_df["cell_name"] + "-" + fofn_df.idx.astype(str)

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

            suffix = "with_reference_help"
            if df.at[wildcards.sample, "reference_fofn"] == "N/A":
                suffix = "from_raw"

            new_fastqz = fofn_df.apply(
                lambda row: f"results/reads_filtered/{wildcards.sample}/fastq/{row.name}_target-reads_{suffix}-subset.fastq.gz",
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
