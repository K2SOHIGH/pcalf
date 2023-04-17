
localrules: dm_ncbi_dehydrated_dataset , dm_ncbi_rehydrate, dm_fallback

wildcard_constraints:
    taxid = "[0-9]*"

# rule dm_ncbi_download_checkpoint: 
#     """
#         dm_ncbi_download_checkpoint it's also a checkpoint because we can loose some assemblies from dm_ncbi_checkpoint.
#         "NCBI dataset" command line tool might failed to retrieve one or more assemblies.
#     """
#     output:
#         touch(
#             os.path.join(RESDIR,"datas","download.checkpoint"),
#             ),
#     input:
#         os.path.join(RESDIR, "datas" , "assemblies",  "report.tsv"),


checkpoint dm_download_checkpoint:
    output:
        touch(os.path.join(RESDIR,"datas", "download.checkpoint")),
    input:
        os.path.join(RESDIR,"datas", "fallback.done"),
    # output:
    #     report = os.path.join(RESDIR, "datas", "report.tsv"),
    #     fail = os.path.join(RESDIR, "datas",  "unreachable.tsv"),
    # input:
    #     report = os.path.join( RESDIR,"datas", "ncbi_dataset", "data","assembly_data_report.tsv" ),
    #     fb = os.path.join(RESDIR,"datas", "fallback.done"),
    # params:
    #     assembly_dir = os.path.join(RESDIR,"datas",  "ncbi_dataset","data"),
    # script:
    #     os.path.join("..","scripts","datas.py")


rule dm_fallback:
    """
        For each assembly that couldn't be downloaded with ncbi-dataset-cli, 
        get URL from NCBI report if exists and try to wget the genome file.
        fallback.txt is a two columns width tabular file with :
            <assembly>.<version> \t <download src>
        where download src is one of :
            - NCBI Datasets CLI // (download from CLI)
            - NCBI FTP server // (download directly from ftp server)
            - not avalaible for download
    """
    output:
        fb = os.path.join(RESDIR,"datas", "fallback.done"),
    input:
        flag       = os.path.join(RESDIR,"datas", "download.done"),
        ftp_report = os.path.join(RESDIR,"reports",
             '{}.report.tsv'.format(date.today().strftime('%Y-%m'))),
        cli_report = os.path.join( RESDIR,"datas", "ncbi_dataset", "data","assembly_data_report.tsv" ),
        accession  = os.path.join(RESDIR,"datas",'accession.filtered.txt'), 
    params:
        assembly_dir = os.path.join(RESDIR,"datas", "ncbi_dataset","data")
    conda:
        os.path.join("..","envs","biopython_1.79.yaml")
    script:
        os.path.join("..","scripts","fallback.py")



rule dm_ncbi_rehydrate:
    """
    Rehydrate ncbi dataset in gzip format. 
    404 errors may occured.
    """
    output:
        os.path.join(RESDIR,"datas",  "download.done"),
    input:
        os.path.join(RESDIR,"datas",  "ncbi_dataset", "fetch.txt"),
    params:
        ddir =  os.path.join( RESDIR, "datas"),
    log:
        os.path.join( RESDIR,"datas","rehydration.log"),
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "datasets rehydrate --gzip --directory {params.ddir} &> {log} || "
        "echo 'missing assemblies detected' ;  "
        "touch {output} "#&& rm -rf {params.ncbi_dataset} "
        #"mv {params.ncbi_dataset}/data/GC* {params.ddir} && "
        


rule dm_ncbi_filter_fetch_summary:
    output:
        os.path.join( RESDIR,"datas" , "ncbi_dataset", "fetch.txt"),
    input:
        os.path.join( RESDIR,"datas",  "ncbi_dataset", "tmp.fetch.txt"),
    params:
        exclude = "| grep -f {} -v ".format(config["exclude_accession"]) if config["exclude_accession"] else "",         
    shell:
        "cat {input} {params.exclude} > {output}"
        


    # run:
    #     l_exclude = []
    #     if params.EXCLUDE_ACCS and os.path.exists(params.EXCLUDE_ACCS):
    #         with open(str(params.EXCLUDE_ACCS)) as stream:
    #             for line in stream.readlines():
    #                 l_exclude.append(line.strip())
    #     l_include = []                     
    #     with open(str(output),'w') as streamout:
    #         with open(str(input),'r') as streamin:
    #             for line in 


# rule dm_ncbi_dataset_fetch_summary: 
#     """
#         Concatenante 'fetch' files from cultured and mag datasets for download.
#     """
#     output:
#         os.path.join( RESDIR,"datas" , "assemblies", "ncbi_dataset", "raw.fetch.txt"),
#     input:
#         expand(
#             os.path.join(RESDIR,"datas",  "ncbi_dataset", "fetch.txt"),
#             mag=["mag","cultured"]
#         )
#     shell:
#         "cat {input} > {output}"
    

# rule dm_ncbi_expand_dataset_table:
#     """
#         Merging datas from cultured and mags datasets.
#     """
#     output:
#         os.path.join(RESDIR, "datas" , "assemblies",  "assembly_data_report.tsv"),
#     input:        
#         expand(
#             os.path.join( RESDIR,"datas", "ncbi_dataset","data","assembly_data_report.tsv" ),
#             mag=["mag","cultured"]
#         ), 
#     run:
#         import pandas as pd
#         ldf = []
#         for f in input:
#             if os.stat(f).st_size != 0:
#                 df = pd.read_csv(f,sep="\t",header=0)
#                 mag = False
#                 if "mag" in f:
#                     mag = True
#                 df["is_mag"] = mag
#                 ldf.append(df)
#         cdf = pd.concat(ldf,axis=0)
#         cdf.set_index("Assembly Accession",inplace=True)        
#         cdf.to_csv(str(output),sep="\t",header=True)



rule dm_ncbi_dataset_table:
    """
        Convert assembly report to table
    """
    output:
        os.path.join( RESDIR,"datas", 
            "ncbi_dataset", "data", "assembly_data_report.tsv" ),
    input:
        os.path.join( RESDIR,"datas", 
            "ncbi_dataset", "data", "assembly_data_report.jsonl" )
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "dataformat tsv genome --force --inputfile {input} > {output} || touch {output}"

rule dm_tmp_fetch:
    output:
        os.path.join( RESDIR,"datas",  "ncbi_dataset", "tmp.fetch.txt"),
    input:
        os.path.join( RESDIR,"datas",  "ncbi_dataset","data","assembly_data_report.jsonl" ),            
    params:
        fetch = os.path.join( RESDIR,"datas",  "ncbi_dataset", "fetch.txt" ),
    shell:
        "mv {params.fetch} {output}"

rule dm_ncbi_dehydrated_dataset:
    """
        Download dehydrated dataset from NCBI based on assembly accession.
        (One for cultured assemblies, anotherone for MAGs.)
    """
    output:
        os.path.join( RESDIR,"datas",  "ncbi_dataset","data","assembly_data_report.jsonl" ),        
    input:        
        os.path.join(RESDIR,"datas", 'accession.filtered.txt'),
    params:        
        outfile = os.path.join( RESDIR,"datas", "ncbi_dataset.zip"),
        outdir  = os.path.join( RESDIR,"datas"),
        outfetch = os.path.join( RESDIR,"datas", "fetch.txt" ),
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "cat {input} | cut -f 1 | sort | uniq | "
        "datasets download genome accession "
        "--inputfile - "
        "--dehydrated "
        "--filename {params.outfile} --assembly-version all --include genome,cds && "
        "unzip -o {params.outfile} -d {params.outdir} || "
        "touch {output} ; "
        
        #rm {params.outfile}; 

