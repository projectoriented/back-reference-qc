rule extract_reads:
    input:
        yak_out_txt="results/read_qv/{sample}/{cell_name}-reference_qv.txt.gz",
        fastq_fai=get_reads(which_one="fai"),
    output:
        yak_filtered="results/reads_filtered/{sample}/{cell_name}-reference_qv-filtered.txt.gz",
        undesirable_reads="results/reads_filtered/{sample}/filtered_out/{cell_name}_undesirable-reads.txt",
        target_reads=temp(
            "results/reads_filtered/{sample}/{cell_name}_target-reads.txt"
        ),
        new_fai="results/reads_filtered/{sample}/fastq/{cell_name}-subset.fastq.gz.fai",
    params:
        comparison_type = get_comparison_type
    log:
        "results/reads_filtered/{sample}/log/{cell_name}-extract_undesirable_reads.log",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=72,
    run:
        from scipy import stats

        # logging
        sys.stdout = open(log[0], "w")

        print(f"comparison_type: {params.comparison_type} for sample: {wildcards.sample}")
        print(f"z_filter: {FILTER_Z} for sample: {wildcards.sample}")
        yak_out_df = pd.read_csv(
            input.yak_out_txt,
            sep="\t",
            comment="C",
            header=None,
            usecols=[1, 5],
            names=["read_name", "qv"],
        ).dropna()
        total_records = yak_out_df.shape[0]

        # Exclude the right-side extremes, yak outputs 99.00 if read vs. short read kmer libraries are 100% matching
        extreme_df = yak_out_df.copy()
        extreme_df = extreme_df.query("qv == 99")
        prcntge = "{:.2f}".format(extreme_df.shape[0] / total_records)

        print(f"QV-99\t{extreme_df['qv'].median()}\t{wildcards.cell_name}\t{prcntge}")

        # Exclude the left-side extremes, yak outputs 0.00
        zeros_df = yak_out_df.copy()
        zeros_df = zeros_df.query("qv == 0")
        prcntge = "{:.2f}".format(zeros_df.shape[0] / total_records)
        print(f"QV-0\t{zeros_df['qv'].median()}\t{wildcards.cell_name}\t{prcntge}")

        # Calculate Z score
        combine_extreme_indicies = zeros_df.index.tolist() + extreme_df.index.tolist()
        yak_out_df = yak_out_df.loc[~yak_out_df.index.isin(combine_extreme_indicies)]
        yak_out_df["z"] = stats.zscore(yak_out_df["qv"])

        if params.comparison_type == "self":
            # Write out a filtered yak output to plot KDE later on
            filter_df = yak_out_df[yak_out_df["z"] > FILTER_Z]

            # Get undesirable reads to filter out
            yak_out_df = yak_out_df[yak_out_df["z"] < FILTER_Z]

            # Table to write out in addition to the undesirable ones
            which_df = zeros_df
        elif params.comparison_type == "other":
            # Write out a filtered yak output to plot KDE later on
            filter_df = yak_out_df[yak_out_df["z"] < abs(FILTER_Z)]

            # Get undesirable reads to filter out
            yak_out_df = yak_out_df[yak_out_df["z"] > abs(FILTER_Z)]

            # Table to write out in addition to the undesirable ones
            which_df = extreme_df
        else:
            raise ValueError(f"Invalid argument for comparison (current: {params.comparison_type}). Choose from either (other, self).")


        filter_df.to_csv(
            output.yak_filtered,
            header=True,
            index=False,
            columns=["read_name", "qv"],
            sep="\t",
        )
        del filter_df

        prcntge = "{:.2f}".format(yak_out_df.shape[0] / total_records)

        print(f"QV-low\t{yak_out_df['qv'].median()}\t{wildcards.cell_name}\t{prcntge}")

        # Write out the low quality reads as well as the ones with zero QV.
        yak_out_df = pd.concat([yak_out_df, which_df]).fillna("N/A")
        yak_out_df.to_csv(output.undesirable_reads,index=False,header=True,sep="\t")

        # Filter out from fai and output the target reads
        fastq_fai_df = pd.read_csv(
            input.fastq_fai,
            sep="\t",
            header=None,
            names=[
                "read_name",
                "length",
                "offset",
                "linebases",
                "linewidth",
                "qualoffset",
            ],
        )

        # Remove the low quality reads from the resulting fastq
        fastq_fai_df = fastq_fai_df[
            ~fastq_fai_df["read_name"].isin(
                yak_out_df.read_name
            )
        ]

        fastq_fai_df.to_csv(
            output.target_reads,
            header=False,
            index=False,
            columns=["read_name"],
            sep="\t",
        )
        fastq_fai_df.to_csv(output.new_fai, header=False, index=False, sep="\t")


if config["new_fastq"]:
    parts = 15

    scattergather:
        sg_parts=parts,

    rule split_fastq:
        input:
            target_reads="results/reads_filtered/{sample}/{cell_name}_target-reads.txt",
        output:
            splitted=temp(
                scatter.sg_parts(
                    "results/reads_filtered/{{sample}}/temp/{{cell_name}}_target-reads_{scatteritem}.txt"
                )
            ),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 16,
            hrs=72,
        run:
            import numpy as np

            df = pd.read_table(input.target_reads, header=None)

            dfs = np.array_split(df, parts)
            for idx, entry in enumerate(dfs):
                entry.to_csv(output.splitted[idx], index=False, header=False)

    rule subset_fastq:
        input:
            target_reads="results/reads_filtered/{sample}/temp/{cell_name}_target-reads_{scatteritem}.txt",
            fastq=get_reads(which_one="cell"),
        output:
            subsetted_fastq=temp(
                "results/reads_filtered/{sample}/temp/{cell_name}-subset_{scatteritem}.fastq"
            ),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 16,
            hrs=72,
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            f"seqtk/{SEQTK_VERSION}",
        shell:
            """
            seqtk subseq {input.fastq} {input.target_reads} > {output.subsetted_fastq}
            """

    rule gather_fastq:
        input:
            subset_fastqs=gather.sg_parts(
                "results/reads_filtered/{{sample}}/temp/{{cell_name}}-subset_{scatteritem}.fastq"
            ),
        output:
            new_fastqz="results/reads_filtered/{sample}/fastq/{cell_name}-subset.fastq.gz",
            merged_fastqs=temp(
                "results/reads_filtered/{sample}/fastq/{cell_name}-subset.fastq"
            ),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 16,
            hrs=72,
        envmodules:
            "modules",
            "modules-init",
            "modules-gs/prod",
            "modules-eichler/prod",
            "samtools/1.14",
        shell:
            """
            cat {input.subset_fastqs} > {output.merged_fastqs} \
                && bgzip --keep {output.merged_fastqs} \
                && samtools fqidx {output.new_fastqz}
            """

    rule make_new_fofn:
        input:
            fastq=get_query_outs(which_one="new_fastqz"),
        output:
            new_fofn="results/reads_filtered/{sample}/fastq.fofn",
        threads: 1
        resources:
            mem=lambda wildcards, attempt: attempt * 16,
            hrs=72,
        shell:
            """
            readlink -f {input.fastq} > {output.new_fofn}
            """
