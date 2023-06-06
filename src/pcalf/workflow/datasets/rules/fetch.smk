localrules: dm_ncbi_accession_from_ncbi_cli_taxid ,dm_ncbi_get_children_taxid, dm_ncbi_get_reports

# remove accessions based on input file : config["exclude_accession"]
rule dm_ncbi:
    output:
        os.path.join(RESDIR,"tmp", 'accession.filtered.txt'),
    input:
        os.path.join(RESDIR,"tmp",'accession.txt'),
    params:
        EXCLUDE_ACCS = config["exclude_accession"],        
    run:
        l_exclude = []
        if params.EXCLUDE_ACCS and os.path.exists(params.EXCLUDE_ACCS):
            with open(str(params.EXCLUDE_ACCS)) as stream:
                for line in stream.readlines():
                    l_exclude.append(line.strip())
        l_include = []                 
        with open(str(input),'r') as stream:
            for line in stream.readlines():
                acc = line.strip()
                cacc = acc.replace("GCA","GCF") if acc.startswith("GCA") else acc.replace("GCF","GCA")

                if acc not in l_exclude:
                    l_include.append(acc)
                if cacc not in l_exclude:
                    l_include.append(cacc)                            
        with open(str(output),'w') as streamout :             
            streamout.write("\n".join(list(set(l_include))))

def get_accession(wildcards):
    files = []
    if config["include_accession"]:
        if os.path.exists(config["include_accession"]):
            files.append( config["include_accession"])
    if config["taxid"]:
        #files.append( os.path.join( RESDIR,"tmp","accession.taxid.txt" ))
        files.append( os.path.join( RESDIR,"tmp","accession.report.txt"))
        files.append( os.path.join( RESDIR,"tmp","accession.cli.txt") )
    return files 

rule dm_ncbi_accession:
    """
        Make accession file depending on workflow input (acc file/taxids). 
        GC(A|F)_XXXXXX. Note that .<version> is removed and duplicates are removed.
    """
    output:
        os.path.join(RESDIR,"tmp","accession.txt"),
    input:
        get_accession,
    shell:
        "cat {input} | sort | uniq > {output}"


# GET ACCESSION FROM TAXID

rule dm_ncbi_aggregate_accession_from_ncbi_cli_taxid:
    output:
        os.path.join(RESDIR,"tmp","accession.cli.txt"),
    input:
        expand(os.path.join(RESDIR,"tmp","taxids", "taxid_{taxid}",'accession_from_taxid_cmdline.tsv'),
            taxid=config["taxid"]),
    shell:
        "cat {input} > {output} "

rule dm_ncbi_accession_from_ncbi_cli_taxid:
    output:
        os.path.join(RESDIR,"tmp","taxids", "taxid_{taxid}",'accession_from_taxid_cmdline.tsv'),
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "datasets summary genome taxon {wildcards.taxid}  --as-json-lines |" 
        "dataformat tsv genome --fields accession | "
        "tail -n +2 > {output}"

rule dm_ncbi_filter_reports:
    '''
        filter assembly based on request taxid/phyla
    '''
    output:
        os.path.join(RESDIR,"tmp",'accession.report.txt'),
    input:
        reports = os.path.join(RESDIR,'tmp',"reports",'{}.report.tsv'.format(date.today().strftime('%Y-%m'))),
        taxids  = os.path.join(RESDIR,"tmp","taxid.txt"),     
    conda:
        "../envs/biopython_1.79.yaml"
    script:
        "../scripts/get_assemblies.py"


# rule dm_ncbi_aggregate_taxid_accessions:
#     " aggregate accession from both sections and for all taxid in a file."
#     output:
#         os.path.join(RESDIR,"tmp","accession.taxid.txt"),
#     input:
#         expand(os.path.join(RESDIR,"tmp",'accession.taxid.{section}.txt'),
#             section=config['section']),
#     shell:
#         "tail -q -n +2 {input} | cut -f 1  | sort | uniq >> {output}"


# rule dm_ncbi_taxid_to_accession:
#     '''
#         get genome accession based on taxid using ncbi-genome-download utility  
#     '''
#     output:
#         os.path.join(RESDIR,"tmp",'accession.taxid.{section}.txt'),
#     input:
#         os.path.join(RESDIR,"tmp","taxid.txt"),
#     params:
#         group = config['group'],
#     conda:
#         "../envs/ngd.yaml"
#     shell:
#         "ncbi-genome-download {params.group} -s {wildcards.section} -T {input} -n -P > {output} || touch {output}"



rule dm_ncbi_expand_taxid:
    output:
        os.path.join(RESDIR,"tmp","taxid.txt"),
    input:
        expand(os.path.join(RESDIR,"tmp","taxids", "taxid_{taxid}",'taxid.txt'),
            taxid=config["taxid"]),
    shell:
        "cat {input} | cut -f 2 | sort | uniq  > {output}"
 

rule dm_ncbi_get_children_taxid:
    '''
        resolve children taxid
    '''
    output:
        taxid = os.path.join(RESDIR,"tmp","taxids", "taxid_{taxid}",'taxid.txt'),
        taxidname = os.path.join(RESDIR,"tmp","taxids", "taxid_{taxid}",'taxid.name.txt'),
    log:
        os.path.join(RESDIR, "logs","taxid_{taxid}.log"),
    params:           
        taxid = lambda wildcards: wildcards.taxid,
    conda:
        "../envs/ngd.yaml"
    script:
        "../scripts/taxid_children.py"   

# DOWNLOAD NCBI REPORTS
rule dm_ncbi_get_reports:
    '''
        merge ncbi reports
    '''
    output:
        os.path.join(RESDIR,"tmp","reports",'{}.report.tsv'.format(date.today().strftime('%Y-%m'))),
    log:
        os.path.join(RESDIR,"tmp","reports",'report.log'),  
    conda:
        "../envs/biopython_1.79.yaml"
    script:
        "../scripts/get_ncbi_reports.py"