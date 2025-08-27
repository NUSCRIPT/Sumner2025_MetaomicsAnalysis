# For configuring files 

BAS_root = "/path/to/my_data/BAS_manuscript/"
BAS_data = "/path/to/genomics_cluster/HartmannLab/jack/bas_pipeline/mlm2/results/"

# METAGENOMICS
FILE_mgx_ko = paste0(BAS_data, "MGX/final_tables/mgx_humann_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv")
FILE_mgx_pfam = paste0(BAS_data, "MGX/final_tables/mgx_humann_genefamilies_cpm_PFAM_unstratified_named.tsv")
FILE_mgx_path = paste0(BAS_data, "MGX/final_tables/mgx_humann_pathabundance_cpm_unstratified.tsv")
FILE_mgx_mpa_all = paste0(BAS_data, "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt")
FILE_mgx_reads = paste0(BAS_data, "MGX/kneaddata/mgx_kneaddata_count_table.tsv")
FILE_mgx_dna_meta = "~/SCRIPT_01_Metagenomics.xlsx"

# METAGENOMICS - vOTUS
FILE_mgx_virome_table = "../virome/vOTU_feature_table.tsv"
FILE_mgx_virome_taxa = "../virome/vOTU_taxa_table.tsv"

# METATRANSCRIPTOMICS
# FILE_mtx_ko = paste0(BAS_data, "MTX/final_tables/mtx_humann_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv")
# FILE_mtx_pfam = paste0(BAS_data, "MTX/final_tables/mtx_humann_genefamilies_cpm_PFAM_unstratified_named.tsv")
# FILE_mtx_pfam = paste0(BAS_data, "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM_unstratified_named.tsv")
# FILE_mtx_bracken = paste0(BAS_data, "MTX/kraken_uniq/merged_bracken_uniq_report_profile.tsv")

FILE_mtx_ko = paste0(BAS_data, "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv")
FILE_mtx_pfam = paste0(BAS_data, "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv")
# FILE_mtx_pfam = paste0(BAS_data, "MTX/final_tables/mtx_humann_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv")
FILE_mtx_path = paste0(BAS_data, "MTX/final_tables/humann_closed/mtx_humann_closed_pathabundance_cpm_unstratified.tsv")
FILE_mtx_mpa_all = paste0(BAS_data, "MTX/metaphlan/mtx_metaphlan_abundance_table_all.txt")
FILE_mtx_bracken = paste0(BAS_data, "MTX/kraken2_standard/merged_bracken_report_profile.tsv")
FILE_mtx_reads = paste0(BAS_data, "MTX/kneaddata/mtx_kneaddata_count_table.tsv")

DIRECTORY_mtx_kraken_uniq = paste0(BAS_data, "MTX/kraken_uniq/")
# AMPLICON


# METADATA

