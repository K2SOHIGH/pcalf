

rule gm_target_gtdbtk:
    # GTDBTK target 
    output:
        temp(os.path.join(RESDIR , "gtdbtk-res","gtdb-tk.done")),
    input:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv"),
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.classify.tree"),
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.bac120.classify.tree"),
        os.path.join(RESDIR , "gtdbtk-res","identify","batch_file.tsv"),           
        # os.path.join(RESDIR , "gtdbtk-res","benchmark.txt"),
    shell:
        "touch {output}"


# def aggregate_benchmarks_files(wildcards):
#     batches = genomes_into_batches(wildcards)
#     return expand(os.path.join(RESDIR , "gtdbtk-res", "benchmarks","align","{batch}.txt"),batch=batches.keys()) + \
#         expand(os.path.join(RESDIR , "gtdbtk-res", "benchmarks","identify","{batch}.txt"),
#             batch=batches.keys())        

# rule gm_expand_gtdbtk_benchmark:
#     output:                
#         os.path.join(RESDIR , "gtdbtk-res", "benchmark.txt")
#     input:
#         aggregate_benchmarks_files,
#         os.path.join(RESDIR , "gtdbtk-res", "benchmarks","classify","classify.txt"),
#     params:
#         benchmarkdir = os.path.join(RESDIR , "gtdbtk-res", "benchmarks")
#     shell:
#         "tail -q -n +2 {input} > {output} "


rule gm_GTDBTK_clean_table:
    output:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv")
    input:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.tsv")
    run:            
        df = pd.read_csv(str(input),sep="\t",header=0,index_col=None)        
        df.user_genome = df.user_genome.str.replace("USER_","")
        df.to_csv(str(output),sep="\t",header=True,index=False)

rule gm_concatenate_archaea_and_bacteria_results:
    output:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.tsv")                                            
    input:
        bac = os.path.join(RESDIR , "gtdbtk-res","gtdbtk.bac120.summary.tsv"),
                        #"classify", 
        ar = os.path.join(RESDIR , "gtdbtk-res", "gtdbtk.ar53.summary.tsv"), 
                        #"classify", 
    shell:
        "cat {input.bac} "                # keep header from first input
        "<(tail -q -n +2 {input.ar}) "  # remove header from other input
        "> {output} "

rule gm_GTDBTK_Classify:
    """ GTDBTk classify step on a set of genomes.

    From documentation:
        Finally, the classify step uses pplacer to find the maximum-likelihood 
        placement of each genome in the GTDB-Tk reference tree. GTDB-Tk 
        classifies each genome based on its placement in the reference tree, 
        its relative evolutionary divergence, and/or average nucleotide 
        identity (ANI) to reference genomes.
    """
    output:
        temp(os.path.join(RESDIR , "gtdbtk-res",  "gtdbtk.ar53.summary.tsv")),
        os.path.join(RESDIR , "gtdbtk-res", "gtdbtk.ar53.classify.tree"),
        temp(os.path.join(RESDIR , "gtdbtk-res", "gtdbtk.bac120.summary.tsv")),
        os.path.join(RESDIR , "gtdbtk-res", "gtdbtk.bac120.classify.tree"),
    input: 
        ##align
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.ar53.msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.ar53.user_msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.bac120.msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.bac120.user_msa.fasta.gz"),           
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.bac120.filtered.tsv"),
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.ar53.filtered.tsv"), 
        ##identify - not used
        os.path.join(RESDIR , "gtdbtk-res","identify", "gtdbtk.bac120.markers_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res","identify", "gtdbtk.translation_table_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res","identify", "gtdbtk.ar53.markers_summary.tsv"),
        ##batch file
        os.path.join(RESDIR , "gtdbtk-res","identify","batch_file.tsv"),        
    conda: 
        os.path.join("..",  "envs","gtdbtk_2.1.yaml")
    resources:
        mem= 100000,
        mem_mb= 100000,
        time= "10-12", # time limit for each job
        nodes= 1,
        cpus_per_task= 15,            
    threads: 
        15
    log:
        os.path.join(RESDIR , "gtdbtk-res","logs","gtdbtk_classify_all.log")
    params:
        extension = config["config-genomes"]["genome_extension"],
        align_dir = os.path.join(RESDIR , "gtdbtk-res"),
        outdir = os.path.join(RESDIR , "gtdbtk-res"),
        batchfile = os.path.join(RESDIR , "gtdbtk-res","identify","batch_file.tsv"),
        gtdbtk_data = config["config-genomes"]["GTDB"],
    benchmark:
        os.path.join(RESDIR , "gtdbtk-res" , "benchmarks","classify","classify.txt")
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_data}; "
        "gtdbtk classify "
        "--mash_db {params.gtdbtk_data} "
        "--batchfile {params.batchfile} "
        "--extension {params.extension} "
        "--align_dir {params.align_dir} "
        "--out_dir {params.outdir} "
        "--cpus {threads} &>> {log} && touch {output} "

rule gm_keepbatch:
    output:
        os.path.join(RESDIR , "gtdbtk-res","identify","batch_file.tsv"),
    input:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","batch_file.tsv"),
    params: 
        p = os.path.join(config["config-genomes"]['merge_with'],"identify",'batch_file.tsv') if config["config-genomes"]['merge_with'] \
            else "" ,
    shell:
        "cat {input} > {output} ; "
        "if [ ! -z {params.p} ] ; then "
        "   cat {params.p} >> {output} ; "
        "fi"


def aggregate_identify_tsv(wildcards):
    batches = genomes_into_batches(wildcards)
    return expand(
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify", "gtdbtk.{{resultfile}}.tsv"),
        batch=batches.keys())



rule gm_merge_batches_tsv_identify:
    output:
        os.path.join(RESDIR , "gtdbtk-res","identify", "gtdbtk.{resultfile}.tsv"),
    input:
        aggregate_identify_tsv,
        # expand(os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify", "gtdbtk.{{resultfile}}.tsv"),
        # batch=batches.keys())
    threads: 
        1
    params:
        first = lambda wildcards, input: input[0] if len(input) > 1 else input,
        p = os.path.join(config["config-genomes"]['merge_with'],"identify","gtdbtk.{resultfile}.tsv") if config["config-genomes"]['merge_with'] \
            else "" ,    
    shell:
        "head -n 1 {params.first} > {output} && "
        "tail -q -n +2 {input} >> {output} && "
        "if [ ! -z {params.p} ] ; then "
        "   tail -q -n +2  {params.p} >> {output} ; "
        "fi && touch {output}"            

def aggregate_align_tsv(wildcards):
    batches = genomes_into_batches(wildcards)
    return expand(
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align", "gtdbtk.{{resultfile}}.tsv"),
        batch=batches.keys())

rule gm_merge_batches_tsv_align:
    output:
        os.path.join(RESDIR , "gtdbtk-res","align", "gtdbtk.{resultfile}.tsv"),
    input:
        aggregate_align_tsv,
        # expand(
        # os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align", "gtdbtk.{{resultfile}}.tsv"),
        # batch=batches.keys())
    threads: 
        1
    params:
        first = lambda wildcards, input: input[0] if len(input) > 1 else input,
        p = os.path.join(config["config-genomes"]['merge_with'],"align","gtdbtk.{resultfile}.tsv") if config["config-genomes"]['merge_with'] \
            else "" ,            
    shell:
        "head -n 1 {params.first} > {output} && "
        "tail -q -n +2 {input} >> {output}  && "     
        "if [ ! -z {params.p} ] ; then "
        "   tail -q -n +2  {params.p} >> {output} ; "
        "fi && touch {output}"

def aggregate_align_fastas(wildcards):
    batches = genomes_into_batches(wildcards)
    return expand(
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.{{resultfile}}.fasta.gz"),
        batch=batches.keys())


rule gm_merge_batches_fasta_align:
    """ Merge fasta files from previous steps, so we only have
    one tree for all the genomes (and not one tree per batch) after
    the Classify step.
    """
    output:
        os.path.join(RESDIR , "gtdbtk-res","align","gtdbtk.{resultfile}.fasta.gz"),
    input:
        aggregate_align_fastas,
        # expand(
        # os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.{{resultfile}}.fasta.gz"),
        # batch=batches.keys())
    threads: 
        1
    params: 
        p = os.path.join(config["config-genomes"]['merge_with'],"align",'gtdbtk.{resultfile}.fasta.gz') if config["config-genomes"]['merge_with'] \
            else "" ,
    shell:
        "cat {input} > {output} ; "
        "if [ ! -z {params.p} ] ; then "
        "   cat {params.p} >> {output} ; "
        "fi"

rule gm_GTDBTK_Align:
    """ GTDBTk align step on a given batch of genomes.

    From documentation:
        The align step concatenates the aligned marker genes and filters the 
        concatenated Multiple Sequence Alignments (MSA) to approximately 
        5,000 amino acids.
    """
    output:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.ar53.msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.ar53.user_msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.ar53.filtered.tsv"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.bac120.msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.bac120.user_msa.fasta.gz"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","align","gtdbtk.bac120.filtered.tsv"),
    input:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.bac120.markers_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.translation_table_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.ar53.markers_summary.tsv"),
    conda: 
        os.path.join("..",  "envs","gtdbtk_2.1.yaml")
    resources:
        mem= 100000,
        mem_mb= 100000,
        time= "7-12", # time limit for each job
        nodes= 1,
        cpus_per_task= 10,           
    threads: 
        10
    log:
        os.path.join(RESDIR , "gtdbtk-res","logs","align","gtdbtk_align_{batch}.log")
    params:
        batch_dir = os.path.join( RESDIR , "gtdbtk-res", "tmp", "{batch}" ),
        gtdbtk_data = config["config-genomes"]["GTDB"],
    benchmark:
        os.path.join(RESDIR , "gtdbtk-res" , "benchmarks","align","{batch}.txt")        
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_data}; "
        "gtdbtk align "
        "--identify_dir {params.batch_dir} "
        "--out_dir {params.batch_dir} "
        "--cpus {threads} "
        "&>> {log} "
        "&& touch {output} " #touch output to avoid missing file i.e ar53 missing files in most cases


rule gm_GTDBTK_Identify:
    """ GTDBTk identify step on a given batch of genomes.
    From documentation:
        The identify step calls genes using Prodigal, and uses HMM models and 
        the HMMER package to identify the 120 bacterial and 122 archaeal marker 
        genes used for phylogenetic inference. Multiple sequence alignments 
        (MSA) are obtained by aligning marker genes to their respective HMM 
        model. 
    """
    output:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.bac120.markers_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.translation_table_summary.tsv"),
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}","identify","gtdbtk.ar53.markers_summary.tsv"),
    input:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}", 'batch_file.tsv')
    conda: 
        os.path.join("..",  "envs","gtdbtk_2.1.yaml")
    resources:
        mem = 75000,
        mem_mb = 75000,
        # time= "7-12", # time limit for each job
        nodes= 1,
        cpus_per_task= 10,          
    threads: 
        20   
    params:
        extension = config["config-genomes"]["genome_extension"],
        batch_dir = os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}"),
        gtdbtk_data = config["config-genomes"]["GTDB"]#os.path.abspath(os.path.join("Database","release89")),
    benchmark:
        os.path.join(RESDIR , "gtdbtk-res" , "benchmarks","identify","{batch}.txt")
    log:
        os.path.join(RESDIR , "gtdbtk-res","logs","taxonomy", "gtdbtk_identify_{batch}.log")
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_data}; "
        "gtdbtk identify "
        "--batchfile {input} "
        "--extension {params.extension} "
        "--out_dir {params.batch_dir} "
        "--cpus {threads} "
        "&>> {log} "


def aggregate_batches_gtdbtk(wildcards):
    batches = genomes_into_batches(wildcards)
    return expand(os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}", 'batch_file.tsv'),batch=batches.keys())


rule gm_expand_batch_file:
    output:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","batch_file.tsv")
    input:
        aggregate_batches_gtdbtk, #expand(os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}", 'batch_file.tsv'),batch=BATCHES.keys())
        # expand(os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}", 'batch_file.tsv'),batch=batches.keys())
    shell:
        "cat {input} > {output}"

def get_genome_file_for_gtdbtk(wildcards):
    genomes = {}
    start_dir = os.path.join(RESDIR,"datas","assemblies", "ncbi_dataset")
    pattern   = "GC*_genomic.fna.gz"     
    for dir,_,_ in os.walk(start_dir,followlinks=True):
        g = glob.glob(os.path.join(dir,pattern))
        if g:            
            gid = os.path.basename(os.path.dirname(g[0]))
            genomes[gid] = g  
    
    batches = _2batches(genomes, config["config-genomes"]["batch_size"])    
    return batches[wildcards.batch]

rule gm_gtdb_bins_into_batches:
    output:
        os.path.join(RESDIR , "gtdbtk-res", "tmp","{batch}", 'batch_file.tsv')
    input:
        lambda wildcards : list(get_genome_file_for_gtdbtk(wildcards).values() ) ,
        # lambda wildcards : list(batches[wildcards.batch].values()),
        #list(get_genome_file_for_gtdbtk(wildcards).values()) #lambda wildcards: [BINS[i] for i in BATCHES[wildcards.batch]],
    threads: 
        1
    params:
        bins_tuple = get_genome_file_for_gtdbtk,
    run:
        with open(str(output),'w') as outfile:
            for bin_id,bin_path in params.bins_tuple.items():
                outfile.write("{}\tUSER_{}\n".format(bin_path,bin_id))
