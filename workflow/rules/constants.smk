# Versions
YAK_VERSION = "f37704a"
KRAKEN2_VERSION="2.1.0"
SEQTK_VERSION="1.4"

FILTER_Z=config.get("z_filter", -2)
KRAKEN2_DB=config["kraken2_db"]

# Config parameters
ILLUMINA_ENDEDNESS=config.get("illumina","pair")
