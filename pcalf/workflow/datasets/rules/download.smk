
localrules: dm_ncbi_dehydrated_dataset , dm_ncbi_rehydrate

wildcard_constraints:
    taxid = "[0-9]*"



checkpoint dm_download_checkpoint:
    output:
        touch(os.path.join(RESDIR, "download.checkpoint")),
    input:
        os.path.join(RESDIR,  "download.done"),
        #os.path.join(RESDIR, "fallback.done"),

# rule dm_fallback:
#     """
#         For each assembly that couldn't be downloaded with ncbi-dataset-cli, 
#         get URL from NCBI report if exists and try to wget the genome file.
#         fallback.txt is a two columns width tabular file with :
#             <assembly>.<version> \t <download src>
#         where download src is one of :
#             - NCBI Datasets CLI // (download from CLI)
#             - NCBI FTP server // (download directly from ftp server)
#             - not avalaible for download
#     """
#     output:
#         fb = os.path.join(RESDIR, "fallback.done"),
#     input:
#         flag       = os.path.join(RESDIR, "download.done"),
#         ftp_report = os.path.join(RESDIR,"reports",
#              '{}.report.tsv'.format(date.today().strftime('%Y-%m'))),
#         cli_report = os.path.join( RESDIR, "ncbi_dataset", "data","assembly_data_report.tsv" ),
#         accession  = os.path.join(RESDIR,'accession.filtered.txt'), 
#     params:
#         assembly_dir = os.path.join(RESDIR, "ncbi_dataset","data")
#     conda:
#         os.path.join("..","envs","biopython_1.79.yaml")
#     script:
#         os.path.join("..","scripts","fallback.py")



rule dm_ncbi_rehydrate:
    """
    Rehydrate ncbi dataset in gzip format. 
    404 errors may occured.
    """
    output:
        os.path.join(RESDIR,  "download.done"),
    input:
        os.path.join(RESDIR,  "ncbi_dataset", "fetch.txt"),
    params:
        ddir =  RESDIR,
    log:
        os.path.join( RESDIR,"rehydration.log"),
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "datasets rehydrate --gzip --directory {params.ddir} &> {log} || "
        "echo 'missing assemblies detected' ;  "
        "touch {output} "#&& rm -rf {params.ncbi_dataset} "
        #"mv {params.ncbi_dataset}/data/GC* {params.ddir} && "
        


rule dm_ncbi_filter_fetch_summary:
    output:
        os.path.join( RESDIR, "ncbi_dataset", "fetch.txt"),
    input:
        os.path.join( RESDIR,  "ncbi_dataset", "tmp.fetch.txt"),
    params:
        exclude = "| grep -f {} -v ".format(config["exclude_accession"]) if config["exclude_accession"] else "",         
    shell:
        "cat {input} {params.exclude} > {output}"
        

rule dm_ncbi_dataset_table:
    """
        Convert assembly report to table
    """
    output:
        os.path.join( RESDIR, 
            "ncbi_dataset", "data", "assembly_data_report.tsv" ),
    input:
        os.path.join( RESDIR, 
            "ncbi_dataset", "data", "assembly_data_report.jsonl" )
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "dataformat tsv genome --force --inputfile {input} > {output} || touch {output}"

rule dm_tmp_fetch:
    output:
        os.path.join( RESDIR,  "ncbi_dataset", "tmp.fetch.txt"),
    input:
        os.path.join( RESDIR,  "ncbi_dataset","data","assembly_data_report.jsonl" ),            
    params:
        fetch = os.path.join( RESDIR,   "ncbi_dataset", "fetch.txt" ),
    shell:
        "mv {params.fetch} {output}"

rule dm_ncbi_dehydrated_dataset:
    """
        Download dehydrated dataset from NCBI based on assembly accession.
        (One for cultured assemblies, anotherone for MAGs.)
    """
    output:
        os.path.join( RESDIR,  "ncbi_dataset","data","assembly_data_report.jsonl" ),        
    input:        
        os.path.join(RESDIR,"tmp", 'accession.filtered.txt'),
    params:        
        outfile  = os.path.join( RESDIR, "ncbi_dataset.zip"),
        outdir   = RESDIR,
        outfetch = os.path.join( RESDIR, "fetch.txt" ),
    retries:
        5
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "cat {input} | cut -f 1 | sort | uniq | "
        "datasets download genome accession "
        "--inputfile - "
        "--dehydrated "
        "--filename {params.outfile} --assembly-version all --include genome,cds && "
        "unzip -o {params.outfile} -d {params.outdir};"
        
        #rm {params.outfile}; 

