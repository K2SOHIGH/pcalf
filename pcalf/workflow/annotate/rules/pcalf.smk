import yaml
import gzip
import pandas  as pd
from Bio import SeqIO

rule target_pcalf:     
    output:
        touch(temp(
            os.path.join(RESDIR , "pcalf" , "pcalf.done" )
        ))
    input:
        os.path.join(RESDIR , "pcalf" , "pcalf.summary.tsv" ),
        os.path.join(RESDIR , "pcalf" , "ccyA.summary.tsv" ),
        os.path.join(RESDIR , "pcalf" , "pcalf.features.tsv"),
        os.path.join(RESDIR , "pcalf" , "reformat.done"),

def parse_location(location):
    loca = location.replace("[","").replace("]","").split("=")[-1]
    strand = "+"
    if "complement" in loca:
        loca = loca.replace("complement(","").replace(")","")        
        strand = "-"
    start,stop=loca.split("..")     
    return start,stop , strand

def parse_ncbi_header(header):
    if "[location=" in header:
        region = header.split(" ")[0].split("|")[-1].split("_cds_")[0]
        desc = header.split(" ")[1:]
        start= None
        stop= None
        partial = None
        pseudo = False
        for d in desc:
            if "location=" in d:                
                start,stop,strand = parse_location(d)                
            elif "partial" in d:
                partial = d.replace("[","").replace("]","").split("=")[-1]
            elif "pseudo" in d:
                pseudo = d.replace("[","").replace("]","").split("=")[-1]
        return region,start,stop,strand,partial,pseudo,"NCBI annotation"
    return None,None,None,None,None,None,None

def parse_prodigal_header(header):
    ">MTBS01000157.1_30 # 28299 # 29294 # -1 # ID=157_30;partial=01;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.501"
    seqid, start, stop , strand , desc = header.split(" # ")
    pseudo = "NA"
    partial = None
    for d in desc.split(";"):
        k,v = d.split("=")
        if k == "partial":
            if v == "01":
                partial = "3'"
            elif v=="10":
                partial = "5'"
            elif v=="11":
                partial = "5',3'"
            else:
                pass
    return "_".join(seqid.replace('>','').split("_")[:-1]) , start, stop, strand, partial, pseudo, "Prodigal prediction"


rule ccya:
    output:
        os.path.join(RESDIR , "pcalf" , "ccyA.summary.tsv" ),
    input:
        os.path.join(RESDIR , "pcalf" , "pcalf.summary.tsv" ),
        os.path.join(RESDIR , "cds_fna_input.tsv"),        
    run:
        files = {}
        with open(str(input[1])) as cds:
            for line in cds.readlines():
                fid, fpath = line.split()
                files[fid] = fpath

        with open(str(output),'w') as streamout:
            streamout.write("sequence_id\tccyA_genomic_region\tccyA_start\tccyA_stop\tccyA_frame\tccyA_partial\tccyA_pseudo\tccyA_src\tccyA_seq\n")
            df = pd.read_csv(str(input[0]),sep="\t",header=0,index_col=0)
            for seqid, row in df.iterrows():
                gid = row.sequence_src
                cds_file = str(files[gid])
                cds_handle = gzip.open(cds_file,"rt")  if cds_file.endswith('.gz') else open(cds_file, "r" )                
                fna_records = SeqIO.to_dict(SeqIO.parse(cds_handle,format="fasta"))
                ccya = fna_records[seqid]
                
                region,start,stop,strand,partial,pseudo,src = parse_ncbi_header(ccya.description)
                
                if start is None:
                    region,start,stop,strand,partial,pseudo,src = parse_prodigal_header(ccya.description)
                #"ccyA\tccyA_genomic_region\tccyA_start\tccyA_stop\tccyA_frame\tccyA_partial\tccyA_pseudo\tccyA_seq\n"
                streamout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    ccya.id,
                    region,
                    start,
                    stop,
                    strand,
                    partial,
                    pseudo,
                    src,
                    str(ccya.seq)
                    )
                )
            

rule pm_aggregate_reformat:
    output:
        touch(temp(
            os.path.join(RESDIR,"pcalf","reformat.done")
        ))
    input:
        expand(
            os.path.join(RESDIR , "pcalf" ,"MSA","{msa}.msa.tsv"),
            msa = ["Gly1","Gly2","Gly3","Glyx3"]
        ),
        os.path.join(RESDIR , "pcalf" , "nter.tsv"),

rule pm_format_msa:
    output:
        os.path.join(RESDIR , "pcalf" ,"MSA","{msa}.msa.tsv"),
    input:
        os.path.join(RESDIR , "pcalf" , "pcalf.summary.tsv" ),
    params:
        msa = os.path.join(RESDIR , "pcalf" ,"MSA","{msa}.msa.fa"),
    run:
        rows = []
        row = []
        with open(str(params.msa)) as fh:            
            for line in fh.readlines():
                if line.startswith(">"):
                    if row:
                        rows.append(row)
                    row = [ line.split()[0].replace(">","") , "" ]
                else:                
                    if row:
                        row[-1]+=line.strip()
                    else:
                        continue
        rows.append(row)
        pd.DataFrame(rows,columns=["sequence_id","sequence"]).to_csv(
            str(output),sep='\t',header=True,index=False
        )


rule pm_format_nter:
    output:
        os.path.join(RESDIR , "pcalf" , "nter.tsv"),
    input:
        os.path.join(RESDIR , "pcalf" , "N-ter-DB.tsv"),
    run:
        import pandas as pd
        df = pd.read_csv(str(input),sep="\t",header=None)
        df.columns = ["nter","sequence_id","sequence"]
        df.to_csv(str(output),sep="\t", header=True,index = False)

rule pm_pcalf:
    output:
        os.path.join(RESDIR , "pcalf" , "pcalf.summary.tsv" ),
        os.path.join(RESDIR , "pcalf" , "pcalf.features.tsv"),
        temp(os.path.join(RESDIR , "pcalf" , "N-ter-DB.tsv")),
        expand(
            os.path.join(RESDIR , "pcalf" ,"HMM","{hmm}.hmm"),
            hmm = ["Gly1","Gly2","Gly3","Glyx3"]
        ),
        expand(
            os.path.join(RESDIR , "pcalf" ,"MSA","{msa}.msa.fa"),
            msa = ["Gly1","Gly2","Gly3","Glyx3"]
        ),
    input:        
        os.path.join(RESDIR , "cds_faa_input.tsv"),
    params:
        outdir = os.path.join(RESDIR , "pcalf"),
        glyx3  = "--glyx3-msa {}".format(config["config-ccya"]["glyx3_msa"]) if config["config-ccya"]["glyx3_msa"] else "",
        gly1   = "--gly1-msa {}".format(config["config-ccya"]["gly1_msa"]) if config["config-ccya"]["gly1_msa"] else "",
        gly2   = "--gly2-msa {}".format(config["config-ccya"]["gly2_msa"]) if config["config-ccya"]["gly2_msa"] else "",
        gly3   = "--gly3-msa {}".format(config["config-ccya"]["gly3_msa"]) if config["config-ccya"]["gly3_msa"] else "",
        nter   = "--nterdb {}".format(config["config-ccya"]["nterdb"]) if config["config-ccya"]["nterdb"] else "",
        evalue = "--glyx3-evalue {}".format(config["config-ccya"]["glyx3-evalue"]) if config["config-ccya"]["glyx3-evalue"] else "",
        coverage = "--glyx3-coverage {}".format(config["config-ccya"]["glyx3-coverage"]) if config["config-ccya"]["glyx3-coverage"] else "",
        glyzip_evalue = "--glyzip-evalue {}".format(config["config-ccya"]["glyzip-evalue"]) if config["config-ccya"]["glyzip-evalue"] else "",
        glyzip_coverage = "--glyzip-coverage {}".format(config["config-ccya"]["glyzip-coverage"]) if config["config-ccya"]["glyzip-coverage"] else "",
        maxite = config["config-ccya"]['max-iteration'],
        iterative_search = '--iterative-search' if config["config-ccya"]['iterative-search'] else '',
        iterative_update = '--iterative-search' if config["config-ccya"]['iterative-update'] else '',
        Z = "-Z {}".format(config["config-ccya"]["Z"]) if config["config-ccya"]["Z"] else "",
        domZ = "--domZ {}".format(config["config-ccya"]["domZ"]) if config["config-ccya"]["domZ"] else "",
    log:
        os.path.join(RESDIR , "logs" , "pcalf.log"),
    resources:        
        mem_mb= 1000,        
        nodes= 1, 
        cpus_per_task = 10,  
    threads:
        10
    shell:
        "pcalf -i {input} "
        "-o {params.outdir} "
        "--force --thread {threads} "
        "{params.glyx3} "
        "{params.gly1} "
        "{params.gly2} "
        "{params.gly3} "
        "{params.nter} "
        "--log {log} "
        "{params.evalue} "
        "{params.coverage} "
        "{params.glyzip_evalue} "
        "{params.glyzip_coverage} "
        "{params.iterative_search} "
        "{params.iterative_update} "
        "--max-iteration {params.maxite} "
        "{params.Z} "
        "{params.domZ} " 
    
def get_cds(wildcards):
    cds = []
    for gid, files in INPUT.items():
        cds.append(files["cds_{}".format(wildcards.cds_format)])
    return cds



rule pm_make_cds_files:
    output:
        os.path.join(RESDIR , "cds_{cds_format}_input.tsv"),
    input:
        get_cds,
    params:
        input_datas = INPUT,
    run:        
        with open(str(output),'w') as streamout:
            for gid , files in params.input_datas.items():                
                streamout.write("{}\t{}\n".format(
                    gid,
                    files["cds_{}".format(wildcards.cds_format)]
                ))
