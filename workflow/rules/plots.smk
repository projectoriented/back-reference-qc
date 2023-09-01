rule generate_kde_plot:
    input:
        read_qv_files=get_query_outs(which_one="qv"),
    output:
        sample_qv_kde="results/plots/{sample}/kde-before_filter.png",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import gzip

        df = pd.DataFrame()
        for f in input.read_qv_files:
            with gzip.open(f, "rt") as infile:
                fn = os.path.basename(f)

                data = infile.readlines()
                read_qvs = []

                for line in data:
                    if line.startswith("SQ"):
                        line_split = line.split()
                        read_name = line_split[1]
                        qv = line_split[5]
                        read_qvs.append({"read_name": read_name, "qv": qv, "cell": fn})

                one_df = pd.DataFrame.from_records(read_qvs)
            df = pd.concat([df, one_df])

        df["cell"] = df["cell"].str.split(".", expand=True)[0]
        # Make sure qv is float
        df = df.astype({"qv": "float64"})

        # Exclude the right-side extremes, yak outputs 99.00 if read vs. short read kmer libraries are 100% matching
        df = df.loc[df["qv"] < 80]

        ax = sns.kdeplot(data=df, x="qv", hue="cell", palette="colorblind")

        # Move legend outside of plot
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

        ax.set(xlabel="QV", ylabel="Gaussian Kernel Density Estimate")
        plt.title(f"Read QVs using its Illumina K-mers\nsample={wildcards.sample}")
        plt.tight_layout()
        plt.savefig(output.sample_qv_kde, dpi=300)


rule kde_filtered:
    input:
        read_qv_files=get_query_outs(which_one="filtered_yak"),
    output:
        sample_qv_kde="results/plots/{sample}/kde-after_filter.png",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        import matplotlib.pyplot as plt
        import seaborn as sns

        df = pd.DataFrame()
        for f in input.read_qv_files:
            one_df = pd.read_table(f, header=0)
            one_df["filepath"] = os.path.basename(f)
            one_df["cell"] = one_df["filepath"].str.split(".", expand=True)[0]
            one_df = one_df[["read_name", "qv", "cell"]]
            df = pd.concat([df, one_df])

        # Make sure qv is float
        df = df.astype({"qv": "float64"})
        # Exclude the right-side extremes, yak outputs 99.00 if read vs. short read kmer libraries are 100% matching
        df = df.loc[df["qv"] < 80]

        ax = sns.kdeplot(data=df, x="qv", hue="cell", palette="colorblind")

        # Move legend outside of plot
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

        ax.set(xlabel="QV", ylabel="Gaussian Kernel Density Estimate")
        plt.title(
            f"Read QVs (filtered) using its Illumina K-mers\nsample={wildcards.sample}"
        )
        plt.tight_layout()
        plt.savefig(output.sample_qv_kde, dpi=300)
