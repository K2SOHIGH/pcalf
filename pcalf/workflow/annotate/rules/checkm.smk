checkm_res = "lineage_wf"

FULLSTATSCOLS = ['Bin Id',
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
 '# unique markers (of 43)',
  '# multi-copy',
   'Taxonomy']



# rule cm_empty_checkm_tables:
#     output:
#         taxo   = temp(os.path.join(RESDIR , "checkm-res","tables","checkM_taxonomy.tsv")),
#         stats  = temp(os.path.join(RESDIR , "checkm-res","tables","checkM_statistics.tsv")),
#         fstats = temp(os.path.join(RESDIR , "checkm-res","tables","checkM_statistics_full.tsv")),
#     params:
#         genomes = INPUT
#     run:
#         import pandas as pd
#         tdf = pd.DataFrame(
#             columns = TAXCOLS[1:],
#             index = list(params.genomes.keys())
#             )
#         tdf.index.name = TAXCOLS[0]
#         tdf.to_csv(str(output.taxo),sep="\t",header=True)

#         sdf = pd.DataFrame(
#             columns = STATSCOLS[1:],
#             index = list(params.genomes.keys())
#             )
#         sdf.index.name = STATSCOLS[0]
#         sdf.to_csv(str(output.stats),sep="\t",header=True)

#         fsdf = pd.DataFrame(
#             columns = FULLSTATSCOLS[1:],
#             index = list(params.genomes.keys())
#             )
#         fsdf.index.name = FULLSTATSCOLS[0]
#         fsdf.to_csv(str(output.fstats),sep="\t",header=True)


# rule cm_target_quick_checkm:
#     output:
#         touch(temp(os.path.join(RESDIR,"checkm-res","checkm.quick.done")))
#     input:
#         rules.cm_empty_checkm_tables.output
    
        
rule cm_target_checkm:
    output:
        temp(os.path.join(RESDIR,"checkm-res","checkm.done")),
    input:
        expand(os.path.join(RESDIR , "checkm-res","tables","checkM_{table}.tsv"),table=["taxonomy","statistics_full","statistics"]),
    params:
        tmp = os.path.join(RESDIR , "checkm-res", "tmp"),
    shell:
        'rm -r {params.tmp} && touch {output}'

def aggregate_batches_checkm(wildcards):    
    return expand(
        os.path.join(RESDIR , "checkm-res", "tmp","{batch}",checkm_res,"checkM_{{table}}.tsv"),
        batch=GENOMESBATCH.keys()
    )


rule cm_concat_tables:
    output:
        protected(os.path.join(RESDIR , "checkm-res","tables","checkM_{table}.tsv"))
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


def get_bins(wildcards):
    return expand(
        os.path.join(RESDIR , "checkm-res", "tmp", 'bins', '{batch}', '{bin}.fna'),
        batch=wildcards.batch , bin = list(GENOMESBATCH[wildcards.batch].keys()) 
    )
    
rule cm_CHECKM_lineage_wf:
    """ Main checkM rule.
        Compute CheckM on a given batch of bins.
    """
    output:
        tsv = os.path.join(RESDIR , "checkm-res", "tmp", "{batch}","lineage_wf","checkM_statistics.tsv"),
    input:
        # Function that returns the locations of genomes associated with {batch}        
        get_bins,
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
        batch_dir = lambda wildcards, input: os.path.dirname(input[0]),
        low_memory = "--reduced_tree" if config["config-genomes"]["low_memory"] else "",
        #tmp = "--tmpdir %s" % RESDIR , "checkm-res", "tmp" if RESDIR , "checkm-res", "tmp" else "",
    log:
        os.path.join(RESDIR , "checkm-res", 'logs', "{batch}", checkm_res , "lineage_wf.log")
    benchmark:
        os.path.join(RESDIR , "checkm-res" , "benchmarks","lineage","{batch}.txt") 
    shell:
        "echo {params.checkm_data}; " #export LC_ALL=C ; 
        "checkm data setRoot {params.checkm_data} &>> {log}; "
        "checkm lineage_wf "
        "--tab_table "          # output as a tab-separated file
        "-f {output.tsv} "      # filename for the output
        "-t {threads} "
        "{params.low_memory} "
        "-x fna "               # extension of the bin files
        "{params.batch_dir}/ "  # (input) directory containing the bin files
        "{params.checkm_dir} "  # (output) directory where to store the results
        "&>> {log} "

def get_genome_file_for_checkm(wildcards):
    return glob.glob(
        os.path.join(RESDIR,"datas","assemblies", "ncbi_dataset", wildcards.bin, "GC*_genomic.fna.gz" )
        )[0]

rule cm_bins_into_batches:
    """ Gunzip a `fasta.gz` into the appropriate batch directory. """
    output:
        temp(os.path.join(RESDIR , "checkm-res", "tmp", 'bins', '{batch}', '{bin}.fna'))
    input:
        lambda wildcards: GENOMESBATCH[wildcards.batch][wildcards.bin],
    params:
        cmd = lambda wildcards,input: "gunzip -c " if str(input).endswith('.gz') else "cat ",
    threads: 
        1 
    shell:
        "{params.cmd} {input} > {output}"
