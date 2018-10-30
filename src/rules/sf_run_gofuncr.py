
rule run_gofuncr:
    input: DATA + 'interim/deseq2_output/condition_treated_results_filtered_tx2gene_full_mouse.csv'
    output:
        results = DATA + 'interim/gofuncr_output/gofuncr_results.csv'
        res_inpt = DATA + 'interim/gofuncr_output/gofuncr_results_input.csv'
        res_anno =  DATA + 'interim/gofuncr_output/gofuncr_results_anno.csv'
        deseq2_filt = DATA + 'interim/deseq2_results_filtered.csv'
    shell: 'rscript prep_gofuncr.R {input} {output.results} {output.res_inpt} {output.res_anno} {output.deseq2_filt}'
