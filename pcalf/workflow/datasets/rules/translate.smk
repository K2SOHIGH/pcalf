rule tm_make_genome_file:
    output:
        os.path.join(RESDIR,"genomes.yaml"),
    input:
        os.path.join(RESDIR,"translation.done"),
    params:
        genome_dir = os.path.join(RESDIR, "ncbi_dataset","data")
    run:
        import yaml
        assemblies = {}
        for g in os.listdir(str(params.genome_dir)):
            if g.startswith("GCA_") or g.startswith("GCF_"):
                f = os.path.join(RESDIR, "ncbi_dataset","data" , g ,"cds_from_genomic.faa.gz" )
                assemblies[g] = {
                    "genome" : glob.glob(
                        os.path.join(RESDIR, "ncbi_dataset","data" , g, "GC*_genomic.fna.gz" )
                    )[0],
                    "cds_fna": os.path.join(RESDIR, "ncbi_dataset","data" , g ,"cds_from_genomic.fna.gz" ),
                    "cds_faa": os.path.join(RESDIR, "ncbi_dataset","data" , g ,"cds_from_genomic.faa.gz" ),
                }
                for f in assemblies[g].values():
                    if not os.path.exists(f):
                        raise FileNotFoundError("Error : {} does not exists".format(f))
        with open(str(output),'w') as fh:
            yaml.dump(assemblies, fh )


def aggregate_genome_cds(wildcards):
    cds = [] 
    genome_dir = os.path.join(RESDIR, "ncbi_dataset","data")
    for g in os.listdir(genome_dir):
        if g.startswith("GCA_") or g.startswith("GCF_"):
            f = os.path.join(RESDIR, "ncbi_dataset","data" , g ,"cds_from_genomic.faa.gz" )
            cds.append(f)
    return cds
    
rule tm_expand_cds:
    output:
        touch(
            temp(os.path.join(RESDIR,"translation.done")),
        ),
    input:
        aggregate_genome_cds,


rule tm_ncbi_cds_translate:
    output:
        os.path.join(RESDIR, "ncbi_dataset","data", "{assembly}","cds_from_genomic.faa.gz" )
    input:
        os.path.join(RESDIR, "ncbi_dataset" ,"data" , "{assembly}","cds_from_genomic.fna.gz" )
    params:
        unzipf = os.path.join(RESDIR, "ncbi_dataset" ,"data", "{assembly}","cds_from_genomic.faa" )
    conda:
       "../envs/seqkit.yaml"
    shell:
        "gunzip -c {input} | seqkit translate > {params.unzipf} && gzip {params.unzipf} "
        

def get_genome_file(wildcards):
    return glob.glob(
        os.path.join(RESDIR, "ncbi_dataset","data" , wildcards.assembly, "GC*_genomic.fna.gz" )
        )[0]
    

def infere_prodigal_proc(wildcards,input):
    total_bases = 0
    if str(input).endswith(".gz"):
        fh = gzip.open
        mode = "rt"
    else:
        fh = open
        mode = "r"
    with fh(str(input),mode) as stream:
        for line in stream.readlines():            
            if not line.startswith(">"):
                total_bases += len(line.strip())
    return "single" if total_bases > 100000 else "meta"

rule tm_ncbi_cds_prediction:
    output:        
        cdsfna = os.path.join(RESDIR, "ncbi_dataset","data","{assembly}","cds_from_genomic.fna.gz"),        
    input:
        get_genome_file,
    params:
        mode = lambda wildcards,input : infere_prodigal_proc(wildcards,input),
        cmd = lambda wildcards,input : "gunzip -c" if str(input).endswith(".gz") \
            else "cat",
        unzipfile = os.path.join(RESDIR, "ncbi_dataset","data","{assembly}","cds_from_genomic.fna")
    conda:
        os.path.join("..","envs", "prodigal-2.6.yaml")
    threads:
        10          
    shell:
        '{params.cmd} {input} | prodigal '        
        '-d {params.unzipfile} '
        '-p {params.mode} -q > /dev/null &&'         
        'gzip {params.unzipfile}'
