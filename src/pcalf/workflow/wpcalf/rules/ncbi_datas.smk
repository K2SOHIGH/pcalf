DATASCOLS = ['Assembly Accession',
 'ANI Best ANI match ANI',
 'ANI Best ANI match Assembly',
 'ANI Best ANI match Assembly Coverage',
 'ANI Best ANI match Type Category',
 'ANI Best ANI match Organism',
 'ANI Best ANI match Type Assembly Coverage',
 'ANI Best match status',
 'ANI Category',
 'ANI Check status',
 'ANI Comment',
 'ANI Declared ANI match ANI',
 'ANI Declared ANI match Assembly',
 'ANI Declared ANI match Assembly Coverage',
 'ANI Declared ANI match Type Category',
 'ANI Declared ANI match Organism',
 'ANI Declared ANI match Type Assembly Coverage',
 'ANI Submitted organism',
 'ANI Submitted species',
 'Annotation BUSCO Complete ',
 'Annotation BUSCO Duplicated ',
 'Annotation BUSCO Fragmented ',
 'Annotation BUSCO Lineage ',
 'Annotation BUSCO Missing ',
 'Annotation BUSCO Single Copy ',
 'Annotation BUSCO Total Count ',
 'Annotation BUSCO Version ',
 'Annotation Count Gene Non-coding',
 'Annotation Count Gene Other',
 'Annotation Count Gene Protein-coding',
 'Annotation Count Gene Pseudogene',
 'Annotation Count Gene Total',
 'Annotation Method',
 'Annotation Name',
 'Annotation Pipeline',
 'Annotation Provider',
 'Annotation Release Date',
 'Annotation Release Version',
 'Annotation Report URL',
 'Annotation Software Version',
 'Annotation Status',
 'Assembly Atypical Is Atypical',
 'Assembly Atypical Warnings',
 'Assembly BioProject Accession',
 'Assembly BioProject Lineage Accession',
 'Assembly BioProject Lineage Parent Accession',
 'Assembly BioProject Lineage Parent Accessions',
 'Assembly BioProject Lineage Title',
 'Assembly BioSample Accession',
 'Assembly BioSample Attribute Name',
 'Assembly BioSample Attribute Value',
 'Assembly BioSample BioProject Accession',
 'Assembly BioSample BioProject Parent Accession',
 'Assembly BioSample BioProject Parent Accessions',
 'Assembly BioSample BioProject Title',
 'Assembly BioSample Description Comment',
 'Assembly BioSample Description Organism Common Name',
 'Assembly BioSample Description Organism Infraspecific Names Breed',
 'Assembly BioSample Description Organism Infraspecific Names Cultivar',
 'Assembly BioSample Description Organism Infraspecific Names Ecotype',
 'Assembly BioSample Description Organism Infraspecific Names Isolate',
 'Assembly BioSample Description Organism Infraspecific Names Sex',
 'Assembly BioSample Description Organism Infraspecific Names Strain',
 'Assembly BioSample Description Organism Name',
 'Assembly BioSample Description Organism Pangolin Classification',
 'Assembly BioSample Description Organism Taxonomic ID',
 'Assembly BioSample Description Title',
 'Assembly BioSample Sample Identifiers Database',
 'Assembly BioSample Sample Identifiers Label',
 'Assembly BioSample Sample Identifiers Value',
 'Assembly BioSample Last updated',
 'Assembly BioSample Models',
 'Assembly BioSample Owner Contact Lab',
 'Assembly BioSample Owner Name',
 'Assembly BioSample Package',
 'Assembly BioSample Publication date',
 'Assembly BioSample Status Status',
 'Assembly BioSample Status When',
 'Assembly BioSample Submission date',
 'Assembly Blast URL',
 'Assembly Description',
 'Assembly Level',
 'Assembly Linked Assembly Accession',
 'Assembly Linked Assembly Type',
 'Assembly Name',
 'Assembly Notes ',
 'Assembly Paired Assembly Accession',
 'Assembly Paired Assembly Name',
 'Assembly Paired Assembly Status',
 'Assembly Refseq Category',
 'Assembly Sequencing Tech',
 'Assembly Status',
 'Assembly Submission Date',
 'Assembly Submitter',
 'Assembly Synonym',
 'Assembly Type',
 'Assembly Stats Contig L50',
 'Assembly Stats Contig N50',
 'Assembly Stats Gaps Between Scaffolds Count',
 'Assembly Stats GC Count',
 'Assembly Stats GC Percent',
 'Assembly Stats Number of Component Sequences',
 'Assembly Stats Number of Contigs',
 'Assembly Stats Number of Scaffolds',
 'Assembly Stats Scaffold L50',
 'Assembly Stats Scaffold N50',
 'Assembly Stats Total Number of Chromosomes',
 'Assembly Stats Total Sequence Length',
 'Assembly Stats Total Ungapped Length',
 'CheckM completeness',
 'CheckM completeness percentile',
 'CheckM contamination',
 'CheckM marker set',
 'CheckM marker set rank',
 'CheckM species tax id',
 'CheckM version',
 'Current Accession',
 'Organelle Assembly Name',
 'Organelle BioProject Accessions',
 'Organelle Description',
 'Organelle Infraspecific Name',
 'Organelle Submitter',
 'Organelle Total Seq Length',
 'Organism Common Name',
 'Organism Infraspecific Names Breed',
 'Organism Infraspecific Names Cultivar',
 'Organism Infraspecific Names Ecotype',
 'Organism Infraspecific Names Isolate',
 'Organism Infraspecific Names Sex',
 'Organism Infraspecific Names Strain',
 'Organism Name',
 'Organism Pangolin Classification',
 'Organism Taxonomic ID',
 'Source Database',
 'Type Material Display Text',
 'Type Material Label',
 'WGS contigs URL',
 'WGS project accession',
 'WGS URL',
 'MAG']

rule dm_ncbi_expand_dataset_table:
    """
        Merging datas from cultured and mags datasets.
    """
    output:
        os.path.join( RESDIR,"ncbi-datas",  "assembly_data_report.tsv"),
    input:        
        expand(
            os.path.join( RESDIR,"ncbi-datas","assembly_report_{mag}.tsv"),         
            mag=["mag","culture"]
        ), 
    run:
        import pandas as pd
        ldf = []
        for f in input:
            if os.stat(f).st_size != 0:
                df = pd.read_csv(f,sep="\t",header=0)
                mag = False
                if "mag" in f:
                    mag = True
                df["is_mag"] = mag
                ldf.append(df)
        if ldf:
            cdf = pd.concat(ldf,axis=0)
        else:
            cds = pd.DataFrame(columns=DATASCOLS)
        cdf.set_index("Assembly Accession",inplace=True)        
        cdf.to_csv(str(output),sep="\t",header=True)
            




rule dm_ncbi_dataset_table:
    """
        Convert assembly report to table
    """
    output:
        os.path.join( RESDIR,"ncbi-datas","assembly_report_{mag}.tsv" )         
    input:
        os.path.join( RESDIR,"ncbi-datas","accession.txt" )
    params:
        mag = lambda wildcards: "exclude" if wildcards.mag == "culture" else "only",
    conda:
        "../envs/ncbi_dataset_14.14.0.yaml"
    shell:
        "datasets summary genome accession --mag {params.mag} --as-json-lines --inputfile {input} | "
        "dataformat tsv genome > {output} || touch {output}"

rule dm_accession:
    """
        Convert assembly report to table
    """
    output:
        os.path.join( RESDIR,"ncbi-datas","accession.txt" )     
    params:
        INPUT    
    run:        
        with open(str(output), 'w' ) as fh:
            for k in INPUT.keys():
                fh.write("{}\n".format(k))

