import os
import hashlib
from datetime import datetime

def _get_checksum(file_name):
    """
    Get the MD5 checksum of a file.
    :param file_name: File to check.
    :return: MD5 checksum.
    """
    # Code from Stack Overflow:
    # http://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file

    hash_md5 = hashlib.md5()

    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)

    return hash_md5.hexdigest()

rule make_md5sum_tab:
    input:
        fofn = rules.make_new_fofn.output.new_fofn,
    output:
        tab = "results/reads_filtered/{sample}/fastq.md5.tsv"
    threads: 1,
    resources:
         mem=120,
         hrs=48,
    run:
        md5_df = pd.read_csv(input.fofn,header=None,names=["CLEANED_PATH"])
        md5_df["CLEANED_MD5"] = md5_df["CLEANED_PATH"].apply(_get_checksum)
        md5_df["CLEANED_SIZE"] = md5_df["CLEANED_PATH"].apply(os.path.getsize)
        md5_df.to_csv(output.tab, sep="\t", index=False)

rule overwrite_sample_fastq_files:
    input:
        cleaned_md5_tab = rules.make_md5sum_tab.output.tab,
    output:
        tab="results/overwrite_records/{sample}.overwrite_records.tab.gz"
    threads: 12,
    resources:
        mem=12,
        hrs=48,
    run:
        fofn_df = get_query_fastq(sample_name=wildcards.sample)
        original_fofn_df = fofn_df.reset_index()
        cleaned_df = pd.read_csv(input.cleaned_md5_tab,sep="\t",header=0)
        cleaned_df["CELL"] = cleaned_df["CLEANED_PATH"].apply(lambda x: x.split('/')[-1].split('.fastq_target')[0])
        original_fofn_df = original_fofn_df.rename(columns={"filepath":"ORIGINAL_PATH"})
        original_fofn_df["CELL"] = original_fofn_df["cell_name"].str.replace(".fastq", "", regex=False)
        original_fofn_df.drop(columns=["cell_name"],inplace=True)

        overwrite_df = pd.merge(original_fofn_df, cleaned_df, on="CELL", how="outer")

        print("Getting MD5s")
        overwrite_df["ORIGINAL_MD5"] = overwrite_df["ORIGINAL_PATH"].apply(_get_checksum)
        overwrite_df["ORIGINAL_SIZE"] = overwrite_df["ORIGINAL_PATH"].apply(os.path.getsize)
        overwrite_df["DATE"] = None
        overwrite_df["STATUS"] = None
        overwrite_df["SAMPLE"] = wildcards.sample
        overwrite_df = overwrite_df[["SAMPLE","CELL", "ORIGINAL_PATH", "ORIGINAL_SIZE", "ORIGINAL_MD5", "CLEANED_PATH", "CLEANED_SIZE", "CLEANED_MD5", "DATE", "STATUS"]]

        for index, row in overwrite_df.iterrows():
            original_path = row.ORIGINAL_PATH
            original_md5 = row.ORIGINAL_MD5
            original_tmp_path = f"{row.ORIGINAL_PATH}.tmp"
            cleaned_path = row.CLEANED_PATH
            cleaned_md5 = row.CLEANED_MD5
            if original_md5 == cleaned_md5:
                overwrite_df["STATUS"] = "Skipped:Identical"
                overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue

            if (overwrite_df['ORIGINAL_PATH'].isnull() | overwrite_df['ORIGINAL_PATH'].isnull()).any():
                overwrite_df["STATUS"] = "Skipped:Missing"
                overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                continue

            # Data overwritting rename (tmp) -> hardlink -> remove tmp

            print(f"mv {original_path} {original_tmp_path}")
            print(f"ln {cleaned_path} {original_path}")
            print(f"rm -f {original_tmp_path}")

#           os.system(f"mv {original_path} {original_tmp_path}")
#           os.system(f"ln {cleaned_path} {original_path}")
#           os.system(f"rm -f {original_tmp_path}")

            overwrite_df["STATUS"] = "Overwritten"
            overwrite_df["DATE"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        overwrite_df.to_csv(output.tab, sep="\t", index=False, compression="gzip")
        
