mkdir data
cp -R MGX/final_tables/ data/
mv data/final_tables/ data/mgx
cp MGX/metaphlan/mgx_metaphlan_abundance_table_*.txt data/mgx/
cp MGX/kneaddata/mgx_* data/mgx/
cp -R MTX/final_tables/ data
mv data/final_tables/ data/mtx
cp MTX/kneaddata/mtx_* data/mtx/
cp MTX/metaphlan/mtx_metaphlan_abundance_table_* data/mtx/
cp MTX/humann/mtx_humann_* data/mtx/
cp MTX/humann_MGXmpa/mtx_humann_* data/mtx/

cp MTX/kraken2_standard/merged_bracken_report_profile.tsv data/mtx/
cp MGX/humann/mgx_humann_* data/mgx/