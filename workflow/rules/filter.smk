rule extract_reads:
    input:
        yak_out_txt="results/read_qv/{sample}/{cell_name}-reference_qv.txt.gz",
        fastq_fai=get_reads(which_one="fai"),
    output:
        yak_filtered="results/reads_filtered/{sample}/{cell_name}-reference_qv-filtered.txt.gz",
        undesirable_reads="results/reads_filtered/{sample}/filtered_out/{cell_name}_undesirable-reads.txt",
        target_reads=temp(
            "results/reads_filtered/{sample}/{cell_name}_target-reads_with_reference_help.txt"
        ),
    params:
        comparison_type = get_comparison_type
    log:
        "results/reads_filtered/{sample}/log/{cell_name}-extract_undesirable_reads.log",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 45,
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

rule extract_reads_from_raw_kraken_output:
    input:
        kraken2_out = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-out.txt.gz"
    output:
        target_reads = "results/reads_filtered/{sample}/{cell_name}_target-reads_from_raw.txt"
    params:
        target_taxa=config.get("ignore_taxid", [9606])
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    run:
        df = pd.read_table(input.kraken2_out,header=None,names=["is_classified", "qname", "taxonomy", "length(bp)"])
        target_taxa_list = params.target_taxa

        # Write out the target read names, taxid 9606 is Homo sapiens
        target_taxas = '|'.join([f"\(taxid {x}\)" for x in target_taxa_list])
        df.loc[
            df["taxonomy"].str.contains(fr"{target_taxas}")
        ].to_csv(
            output.target_reads,columns=["qname"],header=False,index=False
        )
