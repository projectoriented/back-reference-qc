rule lowQ_seq:
    input:
        cell_fq = get_reads(which_one="cell"),
        lowQ_reads = "results/reads_filtered/{sample}/filtered_out/{cell_name}_lowQV-reads.txt"
    output:
        lowQ_fa = temp("results/reads_filtered/{sample}/filtered_out/fasta/{cell_name}_low-quality.fa"),
    threads: config["default"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"seqtk/{SEQTK_VERSION}",
    shell:
        """
        seqtk subseq {input.cell_fq} <( tail -n +2 {input.lowQ_reads} | cut -f1 ) | seqtk seq -a /dev/stdin > {output.lowQ_fa}
        """

rule kraken2:
    input:
        low_quality_fa = "results/reads_filtered/{sample}/filtered_out/fasta/{cell_name}_low-quality.fa",
    output:
        sequence_taxonomy = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-out.txt.gz",
    params:
        kraken2_db = KRAKEN2_DB
    threads: config["kraken2"]["heavy"]["threads"]
    resources:
        mem=lambda wildcards, attempt: attempt * config["kraken2"]["heavy"]["mem"],
        hrs=config["kraken2"]["hrs"],
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"kraken2/{KRAKEN2_VERSION}",
    shell:
        """
        kraken2 --db {params.kraken2_db} {input.low_quality_fa} --use-names --memory-mapping --threads {threads} | cut -f1,2,3,4 | gzip -c > {output.sequence_taxonomy}
        """

rule summarize_kraken2:
    input:
        kraken2_out = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-out.txt.gz"
    output:
        kraken2_summary = "results/reads_filtered/{sample}/filtered_out/kraken2/{cell_name}_kraken2-summary.tsv.gz"
    threads: config["default"]["threads"]
    resources:
        mem= lambda wildcards,attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"],
    run:
        df = pd.read_table(input.kraken2_out, header=None, names=["is_classified", "qname", "taxonomy", "length(bp)"])
        df["cell_name"] = wildcards.cell_name
        if not df.empty:
            df.groupby(["taxonomy", "cell_name"])["qname"]\
                .count()\
                .reset_index()\
                .rename(columns={"qname": "reads_mapped"})\
                .to_csv(output.kraken2_summary, index=False, header=True, columns=["reads_mapped", "taxonomy", "cell_name"], sep='\t')
        else:
            pd.DataFrame(columns=["reads_mapped", "taxonomy", "cell_name"]).to_csv(output.kraken2_summary, index=False, header=True, columns=["taxonomy", "reads_mapped", "cell_name"], sep='\t')

rule merge_summaries:
    input:
        kraken2_summaries = get_kraken2_summaries
    output:
        merged = "results/reads_filtered/{sample}/filtered_out/kraken2/summary.tsv.gz"
    threads: config["default"]["threads"]
    resources:
        mem= lambda wildcards,attempt: attempt * config["default"]["mem"],
        hrs=config["default"]["hrs"],
    run:
        df = pd.concat([pd.read_table(x) for x in input.kraken2_summaries])
        if not df.empty:
            df2 = df.groupby(["taxonomy"])["reads_mapped"].sum().reset_index()
            df3 = df.groupby("taxonomy")["cell_name"].apply(",".join).reset_index()
            df2 = df2.merge(df3,on="taxonomy")
            df2["sample"] = wildcards.sample
            df2.to_csv(output.merged, index=False, header=True, columns=["sample", "reads_mapped", "taxonomy", "cell_name"], sep='\t')
        else:
            df.to_csv(output.merged, index=False, header=True, columns=["sample", "reads_mapped", "taxonomy", "cell_name"], sep='\t')