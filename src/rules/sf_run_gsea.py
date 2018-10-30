
rule prep_gsea:
    input: 
        deseq2_res = DATA + 'interim/deseq2_results_filtered.csv'
        homologs = DATA + 'interim/annotations/human_mouse_homolog.out' 
    output: DATA + 'interim/mice_human_homolog_gsea_input.rnk'
    shell: 'rscript prep_gsea.R {input.deseq2_res} {input.homologs} {output}'

rule run_gsea:
    input:
        rnk = DATA + 'interim/mice_human_homolog_gsea_input.rnk'
        gmt = DATA + 'interim/annotations/critical_gsea.gmt'
    output: DATA + 'processed/'
    shell: 'java -Xmx512m xtools.gsea.GseaPreranked -gmx {input.gmt} -norm meandiv -nperm 1000 -rnk {input.rnk} -scoring_scheme weighted -rpt_label mice_gsea -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report true -out {output} -gui false'
