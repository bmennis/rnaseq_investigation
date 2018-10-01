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
    shell: 'Rscript {SCRIPTS}HumanRNAseq_step1.R'
