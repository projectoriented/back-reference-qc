import pandas as pd

if config["new_fastq"]:
    parts = 15

    scattergather:
        sg_parts=parts,

    rule split_fastq:
        input:
            target_reads="results/reads_filtered/{sample}/{cell_name}_target-reads.txt",
            kraken2_target_reads = "results/reads_filtered/{sample}/retain/kraken2/{cell_name}_target-reads.txt"
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

            if os.path.getsize(input.kraken2_target_reads) > 0:
                kraken_df = pd.read_table(input.kraken2_target_reads, header=None)
            else:
                kraken_df = pd.DataFrame()

            df = pd.read_table(input.target_reads, header=None)
            df = pd.concat([df, kraken_df])

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
