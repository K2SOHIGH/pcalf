localrules: dm_ncbi_dataset_table

rule dm_ncbi_dataset_table:
    """
        Retrieve a subset of metadatas from NCBI if name is suitable for this ... (i.e GCF_903937765.1)
    """
    output:
        # naming the output file as `out` is needed , see ncbi_metadatas.py
        out = os.path.join( RESDIR,"ncbi-datas","ncbi_metadatas.tsv" )         
    input:
        # naming the input file as `accession` is needed , see ncbi_metadatas.py
        accession = os.path.join( RESDIR,"ncbi-datas","accession.txt")         
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    script:
        "../scripts/ncbi_metadatas.py"

rule dm_accession:
    """
        Store genome's names in a file
    """
    output:
        os.path.join( RESDIR,"ncbi-datas","accession.txt" )     
    params:
        INPUT    
    run:        
        with open(str(output), 'w' ) as fh:
            for k in INPUT.keys():                
                fh.write("{}\n".format(k))                