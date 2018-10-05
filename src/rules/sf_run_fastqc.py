"""Run fastqc on rna seq fastq files and output to interim data directory
"""

include: 'const.py'

rule run_fastqc:
    input: expand(DATA + 'raw/fq_files/{sample}_{read}.fq.gz', sample=('S1','S5','S11','S12','S13','S15','S16','S22','S23','S31','S33','S34'), read=(1,2))
    output: directory(DATA + 'interim/fastqc_files/')
    shell: 'fastqc -o {output} {input}'

rule run_multiqc:
    input: directory(DATA + 'interim/fastqc_files/')
    output: directory(DATA + 'interim/fastqc_html/')
    shell: 'multiqc -o {output} {input}' 

rule prep_deseq2:
    input:
        fq_dir = directory(DATA + 'raw/fq_files/'),
        gene_ann = DATA + 'interim/annotations/{gene_in}.txt',
        samples = DATA + 'interim/annotations/directory_sample_filtered.txt',
    output: 
        csv = DATA + 'interim/deseq2_output/condition_treated_results_filtered_{gene_in}.csv',
        pca = DATA + 'interim/deseq2_output/rna_seq_PCA_filtered_{gene_in}.pdf',
        cluster = DATA + 'interim/deseq2_output/rna_seq_cluster_filtered_{gene_in}.pdf'
    shell: 'Rscript {SCRIPTS}HumanRNAseq_step1.R {input.fq_dir} {input.samples} {input.gene_ann} {output.csv} {output.pca} {output.cluster}'

rule run_all_deseq2:
    input: expand(DATA + 'interim/deseq2_output/rna_seq_cluster_filtered_{gene_in}.pdf', gene_in = ('tx2gene_full_mouse', 'tx2gene_mouse'))
