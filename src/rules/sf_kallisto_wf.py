##########################################
#Specify directory and kallisto index file
##########################################

FASTQDIR=config['FASTQDIR']
KINDEX=config['KINDEX']
OUTPUTFILE=config['OUTPUTFILE']
SCRIPTSDIR=config['SCRIPTSDIR']
NUMTHREADS=config['NUMTHREADS']

##########################################

#Define directories
IDS, = glob_wildcards(FASTQDIR+"{id}_1.fq.gz")
workdir: FASTQDIR

# a pseudo-rule that collects the target files
rule all:
    input: OUTPUTFILE

# Running kallisto
rule runKallisto:
    input:  pair1 = "{name}_1.fq.gz",  pair2="{name}_2.fq.gz"
    output: "{name}/abundance.tsv"
    benchmark: "benchmarks/{name}.json"
    threads:int(NUMTHREADS)
    shell:  "kallisto quant -t "+NUMTHREADS+" -i "+KINDEX+" -o "+FASTQDIR+"{wildcards.name} "+FASTQDIR+"{input.pair1} "+FASTQDIR+"{input.pair2}"

# Rename file and move it up a directory so we can collect it
rule rename:
    input: "{name}/abundance.tsv"
    output: "{name}.tsv"
    shell: "mv {input} {output}"

# Collect all the files and create FPKM Matrix
rule createTpmDF:
    input: expand("{name}.tsv", name=IDS)
    output: OUTPUTFILE
    shell : "Rscript "+SCRIPTSDIR+"merge_kallisto_res.R "+FASTQDIR+" "+OUTPUTFILE

