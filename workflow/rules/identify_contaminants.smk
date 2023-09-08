rule undesirable_seq:
    input:
        cell_fq = get_reads(which_one="cell"),
        undesirable_reads = "results/reads_filtered/{sample}/filtered_out/{cell_name}_undesirable-reads.txt"
    output:
        undesirable_fa = temp("results/reads_filtered/{sample}/filtered_out/fasta/{cell_name}_undesirable-quality.fa"),
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
        seqtk subseq \
            {input.cell_fq} <( tail -n +2 {input.undesirable_reads} | cut -f1 ) \
            | seqtk seq -a /dev/stdin > {output.undesirable_fa}
        """

rule kraken2:
    input:
        undesirable_quality_fa = "results/reads_filtered/{sample}/filtered_out/fasta/{cell_name}_undesirable-quality.fa",
    output:
        sequence_taxonomy = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-out.txt.gz",
    params:
        kraken2_db = KRAKEN2_DB
    threads: 16
    resources:
        mem=lambda wildcards, attempt: attempt * 8,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"kraken2/{KRAKEN2_VERSION}",
    shell:
        """
        kraken2 \
            --db {params.kraken2_db} {input.undesirable_quality_fa} \
            --use-names \
            --memory-mapping \
            --threads {threads} \
            | cut -f1,2,3,4 \
            | gzip -c > {output.sequence_taxonomy}
        """

rule summarize_kraken2:
    input:
        kraken2_out = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-out.txt.gz",
        cell_fai = get_reads(which_one="fai")
    output:
        kraken2_summary = temp("results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-summary.tsv.gz")
    threads: 1
    resources:
        mem= lambda wildcards,attempt: attempt * 16,
        hrs=72,
    run:
        df = pd.read_table(input.kraken2_out, header=None, names=["is_classified", "qname", "taxonomy", "length(bp)"])
        df["cell_name"] = wildcards.cell_name
        df["sample"] = wildcards.sample
        df["cell_total_reads"] = pd.read_table(input.cell_fai, header=None)[0].shape[0]
        if not df.empty:
            df.groupby(["taxonomy", "cell_name", "sample", "cell_total_reads"])["qname"]\
                .count()\
                .reset_index()\
                .rename(columns={"qname": "reads_mapped"})\
                .to_csv(output.kraken2_summary, index=False, header=True, columns=["sample", "reads_mapped", "taxonomy", "cell_name", "cell_total_reads"], sep='\t')
        else:
            pd.DataFrame(columns=["sample", "reads_mapped", "taxonomy", "cell_name", "cell_total_reads"]).to_csv(output.kraken2_summary, index=False, header=True, sep='\t')

rule merge_summaries:
    input:
        kraken2_summaries = get_kraken2_summaries
    output:
        merged = "results/reads_filtered/{sample}/filtered_out/kraken2/summary.tsv.gz"
    threads: 1
    resources:
        mem=lambda wildcards,attempt: attempt * 16,
        hrs=72,
    run:
        df = pd.concat([pd.read_table(x, header=0) for x in input.kraken2_summaries])
        if not df.empty:
            # Get sum of total reads in all cells
            total_cell_reads = sum(df.cell_total_reads.unique())

            df2 = df.groupby(["taxonomy"])["reads_mapped"].sum().reset_index()
            df3 = df.groupby("taxonomy")["cell_name"].apply(",".join).reset_index()
            df2 = df2.merge(df3,on="taxonomy")
            df2["sample"] = wildcards.sample

            # Get sum of undesirable reads mapped to taxonomies
            total_reads_mapped = sum(df2.reads_mapped)

            # Add a total row
            df2 = pd.concat([df2, pd.DataFrame({"sample": ["total"], "reads_mapped": total_reads_mapped, "taxonomy": ["N/A"], "cell_name": total_cell_reads})])

            df2.to_csv(output.merged, index=False, header=True, columns=["sample", "reads_mapped", "taxonomy", "cell_name"], sep='\t')
        else:
            df.to_csv(output.merged, index=False, header=True, columns=["sample", "reads_mapped", "taxonomy", "cell_name"], sep='\t')