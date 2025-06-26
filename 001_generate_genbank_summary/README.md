
1. Setup the env: install `uv` and then use `uv run script.py`. Check [UV](https://docs.astral.sh/uv/) to see how to use this tool

2. `001_generate_genbank_summary/gene_name_sum_all.py` will generate the gene name from the Genbank file, you need to run the script as follows 
```bash
python gene_name_sum_all.py -i sequence.gb -o sequence.gene_name.txt -s sequence.strange.txt
```

1. you will get `sequence.gene_name.txt` , which contain the gene name in the genbank file. This is important due to the fact that people name the same gene in different way, such as COI, CO1, COX1 all represent the Cytochrome c oxidase I (COX1) gene. 
2. Then you need to edit the `sequence.gene_name.txt` to format like bellow. which will be used as input for `001_generate_genbank_summary/Super_MG_process_all.py`

```
>rrnS
12s ribosomal rna
s-rrna
>rrnL
16s ribosomal rna
l-rrna
```

4. you then need to create two file which describe the information that you want to keep in the final output

- This is used to define columns in the csv file (`001_generate_genbank_summary/col_names.txt`)

```
GI
Orgn
Taxa_info
Description
Journal
PMED
title
Author
rrnS_strand
rrnS_seq
rrnL_strand
rrnL_seq
```

- This is used to define how the fasta file will be structured. At this moment, you can keep two fields in the header of the fasta file, GI and Orgn is used in this example bellow (`001_generate_genbank_summary/for_fasta.txt`)

```
GI
Orgn
rrnS_seq
rrnL_seq
```

5. Final step: run the following code to get the csv table (`001_generate_genbank_summary/test_dir/test_table.csv`) that summary the genbank file and the multi-seq fasta file (`001_generate_genbank_summary/test_dir/*.fasta`) of the gene that you selected

```bash
python Super_MG_process_all.py -i sequence.gb -g sequence.gene_name.txt -t test_table -o test_dir -c col_names.txt -f for_fasta.txt
```

