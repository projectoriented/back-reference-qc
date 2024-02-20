rule buildmer:
    input:
        reads=get_reads(which_one="fastq_files"),
    output:
        hash_table=temp("results/hash_table/{sample}/{input_type}.yak"),
    params:
        sr_endedness=ILLUMINA_ENDEDNESS,
        data_type=lambda wildcards: config[f"{wildcards.input_type}_data_type"]
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"yak/{YAK_VERSION}",
    shell:
        """
        if [[ {params.data_type} == "long" ]]; then

            # build k-mer hash tables for high-coverage reads; discard singletons
            yak count -b 37 -t {threads} -o {output.hash_table} {input.reads}

        elif [[ {params.data_type} == "short" ]]; then

            if [[ {params.sr_endedness} == "single" ]]; then
                yak count -b 37 -t {threads} -o {output.hash_table} <(zcat {input.reads})
            else
                # for paired end: to provide two identical streams
                yak count -b 37 -t {threads} -o {output.hash_table} <(zcat {input.reads}) <(zcat {input.reads})
            fi
        else
            echo "Invalid arguments data_type: {params.data_type} and short read endedness: {params.sr_endedness}"
            exit 1

        fi
        """

rule compute_read_qv:
    input:
        query_read=get_reads(which_one="cell"),
        reference_hash_table="results/hash_table/{sample}/reference.yak",
    output:
        qv_txt="results/read_qv/{sample}/{cell_name}-reference_qv.txt.gz",
    threads: 16
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"yak/{YAK_VERSION}",
    shell:
        """
        # compute read QV
        rsync -av {input.query_read} {resources.tmpdir}
        rsync -av {input.reference_hash_table} {resources.tmpdir}
        yak qv -t {threads} -p {resources.tmpdir}/$( basename {input.reference_hash_table} ) {resources.tmpdir}/$( basename {input.query_read} ) | gzip -c > {output.qv_txt}
        """


rule compute_kmer_qv:
    input:
        query_hash_table="results/hash_table/{sample}/query.yak",
        reference_hash_table="results/hash_table/{sample}/reference.yak",
    output:
        kmer_qv="results/read_qv/{sample}/query-reference_kqv.txt.gz",
    threads: 1
    resources:
        mem=calc_mem_gb,
        hrs=72,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"yak/{YAK_VERSION}",
    shell:
        """
        # make sure directories exist
        mkdir -p $( dirname {output.kmer_qv} )

        # compute k-mer QV for reads
        yak inspect {input.query_hash_table} {input.reference_hash_table} | gzip -c > {output.kmer_qv}
        """
