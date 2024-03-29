checkm_res = "lineage_wf"

FULLSTATSCOLS = ['Bin Id',
 'Genome',
 'Marker lineage',
 '# genomes',
 '# markers',
 '# marker sets',
 'Completeness',
 'Contamination',
 'Strain heterogeneity',
 'Genome size (bp)',
 '# ambiguous bases',
 '# scaffolds',
 '# contigs',
 'N50 (scaffolds)',
 'N50 (contigs)',
 'Mean scaffold length (bp)',
 'Mean contig length (bp)',
 'Longest scaffold (bp)Longest contig (bp)',
 'GC',
 'GC std (scaffolds > 1kbp)',
 'Coding density',
 'Translation table',
 '# predicted genes',
 '0',
 '1',
 '2',
 '3',
 '4',
 '5+']

STATSCOLS = ['Bin Id',
 'Genome',
 'Marker lineage',
 '# genomes',
 '# markers',
 '# marker sets',
 '0',
 '1',
 '2',
 '3',
 '4',
 '5+',
 'Completeness',
 'Contamination',
 'Strain heterogeneity']

TAXCOLS = ['Bin Id',
 'Genome',
 '# unique markers (of 43)',
  '# multi-copy',
   'Taxonomy']
    
        
rule cm_target_checkm:
    output:
        temp(os.path.join(RESDIR,"checkm-res","checkm.done")),
    input:
        expand(os.path.join(RESDIR , "checkm-res","tables","checkM_{table}.tsv"),table=["taxonomy","statistics_full","statistics"]),os.path.join(RESDIR , "checkm-res", "tmp", "checkm.done")
    params:
        tmp = os.path.join(RESDIR , "checkm-res", "tmp"),
    shell:
        'touch {output}'
        #'rm -r {params.tmp} && 

def aggregate_batches_checkm(wildcards):    
    return expand(
        os.path.join(RESDIR , "checkm-res", "tmp","{batch}",checkm_res,"checkM_{{table}}.tsv"),
        batch=GENOMESBATCH.keys()
    )


rule cm_concat_tables:
    output:
        os.path.join(RESDIR , "checkm-res","tables","checkM_{table}.tsv")
    input:
        aggregate_batches_checkm,         
    shell:
        "head -n 1 {input[0]} > {output} && "
        "tail -q -n +2 {input} >> {output} "

rule cm_CHECKM_taxonomy:
    """ Create an extended taxonomy. """
    output:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp", "{batch}", checkm_res ,"checkM_taxonomy.tsv"),
    input:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp", "{batch}", checkm_res ,"checkM_statistics.tsv"),
    conda: 
        os.path.join("..", "envs","checkm.yaml")
    params:
        checkm_dir = lambda wildcards, output: os.path.dirname(output.tsv),
        checkm_data = config["config-genomes"]["CheckM_data"],
    log:
        os.path.join(RESDIR , "checkm-res", 'logs', "{batch}", checkm_res , "taxonomy.log")
    shell:
        "echo \"{params.checkm_data}\" | checkm data setRoot {params.checkm_data} &> {log}; "
        "checkm tree_qa "
        "--tab_table "              # output as tab-separated file
        "-f {output.tsv} "          # output filename
        "{params.checkm_dir} "      # output folder
        "&>> {log} "


rule cm_CHECKM_qa_full:
    """ Create an extended output with more statistics based on previous
    computation.
    
    It only uses `checkm qa` so it will be fast and not recompute everything.

    """
    output:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp","{batch}", checkm_res , "checkM_statistics_full.tsv"),
    input:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp", "{batch}", checkm_res , "checkM_statistics.tsv"),
    conda: 
        os.path.join("..",  "envs","checkm.yaml")
    threads: 
        10
    params:
        checkm_dir = lambda wildcards, output: os.path.dirname(output.tsv),
        checkm_data = config["config-genomes"]["CheckM_data"],
        markerfile =  lambda wildcards, output: os.path.dirname(output.tsv) + "/lineage.ms",
    log:
        os.path.join(RESDIR , "checkm-res", 'logs', "{batch}", checkm_res , "stats_full.log")
    benchmark:
        os.path.join(RESDIR , "checkm-res", "benchmarks","qa_full","{batch}.txt")
    shell:
        "echo \"{params.checkm_data}\" | checkm data setRoot {params.checkm_data} &>> {log};"
        "checkm qa "
        "-o 2 "                             # extended summary of bin stats
        "--tab_table "                      # output as tab-separated file
        "-f {output.tsv} "                  # output filename
        "-t {threads} "
        "{params.markerfile} "              # <marker file>
        "{params.checkm_dir} "              # <output folder>
        "&>> {log} "


def get_extension_from_file(wildcards,input):
    with open(str(input)) as f:
        return "."+os.path.basename(f.read().strip()).split('.')[-1]

rule cm_CHECKM_lineage_wf:
    """ Main checkM rule.
        Compute CheckM on a given batch of bins.
    """
    output:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp", "{batch}","lineage_wf","checkM_statistics.tsv"),
    input:             
        os.path.join(RESDIR , "checkm-res", "tmp", 'bins', '{batch}', 'batch.tsv'),
    conda: 
        os.path.join(".." , "envs" , "checkm.yaml")
    resources:
        mem= 16000 if config["config-genomes"]["low_memory"] else 50000,
        mem_mb= 20000 if config["config-genomes"]["low_memory"] else 70000,
        # time= 24:00:00, # time limit for each job
        nodes= 1,
        cpus_per_task = 10 if config["config-genomes"]["low_memory"] else  15 ,           
    threads: 
        10
    params:            
        checkm_dir = lambda wildcards, output: os.path.dirname(output.tsv),
        checkm_data = config["config-genomes"]["CheckM_data"],
        low_memory = "--reduced_tree" if config["config-genomes"]["low_memory"] else "",
        tmp = os.path.join(RESDIR , "checkm-res", "tmp", '{batch}', 'tmp' ),
        ext = lambda wildcards,input: get_extension_from_file(wildcards,input),
    log:
        os.path.join(RESDIR , "checkm-res", 'logs', "{batch}", checkm_res , "lineage_wf.log")
    benchmark:
        os.path.join(RESDIR , "checkm-res" , "benchmarks","lineage","{batch}.txt") 
    shell:
        #"mkdir -p {params.tmp}; "
        "echo {params.checkm_data}; " #export LC_ALL=C ; 
        "checkm data setRoot {params.checkm_data} &>> {log}; "
        "checkm lineage_wf "
        "--tab_table "          # output as a tab-separated file
        "-f {output.tsv} "      # filename for the output
        "-t {threads} "
        "{params.low_memory} "        
        "-x {params.ext} "               # extension of the bin files
        "{input} "  
        "{params.checkm_dir} "  # (output) directory where to store the results
        "&>> {log}; "
        


# rule cm_bins_into_batches:
#     """ Gunzip a `fasta.gz` into the appropriate batch directory. """
#     output:
#         temp(os.path.join(RESDIR , "checkm-res", "tmp", 'bins', '{batch}', '{bin}.fna'))
#     input:
#         lambda wildcards: GENOMESBATCH[wildcards.batch][wildcards.bin],
#     params:
#         cmd = lambda wildcards,input: "gunzip -c " if str(input).endswith('.gz') else "cat ",
#     threads: 
#         1 
#     shell:
#         "{params.cmd} {input} > {output}"

rule cm_bins_into_batches:
    output:
        os.path.join(RESDIR , "checkm-res", "tmp", 'bins', '{batch}', 'batch.tsv')
    params:
        batch = lambda wildcards: GENOMESBATCH[wildcards.batch],
    run:
        with open(str(output),'w') as fh:
            for label,genome_path in params.batch.items():                
                fh.write('{}\t{}\n'.format(
                    label, genome_path
                ))

def crop_nth_iteration (str,n,symbol):
    pos=str.find(symbol)
    while pos>=0 and n>1:
        pos=str.find(symbol,pos+len(symbol))
        n-=1
    if pos>=0:
        return(str[:pos])
    else:
        return(str)

rule add_genome_col:
    output:
        touch(temp(os.path.join(RESDIR , "checkm-res", "tmp", "checkm.done"))),
    input: 
        expand(os.path.join(RESDIR , "checkm-res","tables","checkM_{table}.tsv"),table=["taxonomy","statistics_full","statistics"]),
    run:
        for f in input:
            pd.set_option('display.max_columns', None)
            df=pd.read_csv(f,sep='\t')
            df['Genome']=df['Bin Id'].apply(lambda x:crop_nth_iteration(x,2,'_'))
            df.insert(1,'Genome',df.pop('Genome'))
            df.to_csv(f,sep='\t',index=False)

