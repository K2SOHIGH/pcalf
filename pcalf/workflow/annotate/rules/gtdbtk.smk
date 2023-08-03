GTDBSUMMARYCOLS = ['user_genome',
 'classification',
 'fastani_reference',
 'fastani_reference_radius',
 'fastani_taxonomy',
 'fastani_ani',
 'fastani_af',
 'closest_placement_reference',
 'closest_placement_radius',
 'closest_placement_taxonomy',
 'closest_placement_ani',
 'closest_placement_af',
 'pplacer_taxonomy',
 'classification_method',
 'note',
 'other_related_references(genome_id,species_name,radius,ANI,AF)',
 'msa_percent',
 'translation_table',
 'red_value',
 'warnings']



rule gm_target_gtdbtk:
    output:
        touch(temp(os.path.join(RESDIR,"gtdbtk-res","gtdbtk.done")))
    input:
        os.path.join(RESDIR, "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv"),
    
rule gm_GTDBTK_clean_table:
    output:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv")
    input:
        expand(
            os.path.join(RESDIR , "gtdbtk-res", "{batch}", "gtdbtk.ar53.bac120.summary.tsv")
            ,batch = GENOMESBATCH.keys())
    run:        
        import os
        import pandas as pd
        dfs = []
        for f in input:
            if os.stat(str(f)).st_size != 0:
                t = pd.read_csv(str(f),sep='\t',header=0)
                dfs.append(t)    
        df = pd.DataFrame()
        if dfs:
            df = pd.concat(dfs,axis=0)            
        df.to_csv(str(output),index=False,header=True,sep='\t')


rule gm_concatenate_archaea_and_bacteria_results:
    output:
        os.path.join(RESDIR, "gtdbtk-res", "{batch}", "gtdbtk.ar53.bac120.summary.tsv")                                            
    input:
        bac = os.path.join(RESDIR, "gtdbtk-res", "{batch}", "gtdbtk.bac120.summary.tsv"),
        ar  = os.path.join(RESDIR, "gtdbtk-res", "{batch}", "gtdbtk.ar53.summary.tsv"  ),                         
    run:
        import os
        import pandas as pd
        dfs = []
        for f in input:
            if os.stat(str(f)).st_size != 0:
                t = pd.read_csv(str(f),sep='\t',header=0)
                dfs.append(t)    
        df = pd.DataFrame()
        if dfs:
            df = pd.concat(dfs,axis=0)
            df.user_genome = df.user_genome.str.replace("USER_","")
        df.to_csv(str(output),index=False,header=True,sep='\t')



rule gm_gtdbtk_classify_wf:
    output:
        os.path.join(RESDIR , "gtdbtk-res", "{batch}", "gtdbtk.ar53.summary.tsv"),        
        os.path.join(RESDIR , "gtdbtk-res", "{batch}", "gtdbtk.bac120.summary.tsv"),
    input:
        os.path.join(RESDIR , "gtdbtk-res", "{batch}", "batchfile.tsv"),
    conda: 
        os.path.join(".." , "envs" , "gtdbtk_2.1.yaml")
    params:
        outdir = os.path.join(RESDIR , "gtdbtk-res", "{batch}"),
        gtdbtk_data = config["config-genomes"]["GTDB"],
    threads:
        15
    shell:
        "export GTDBTK_DATA_PATH={params.gtdbtk_data}; "
        "gtdbtk classify_wf "
        "--mash_db {params.gtdbtk_data} "
        "--batchfile {input} "
        "--out_dir {params.outdir} "
        "--cpus {threads} --debug && "
        "touch {output}"


#def get_batch(wildcards):


rule gm_gtdbtk_batchfile:
    output:
        os.path.join(RESDIR , "gtdbtk-res", "{batch}", "batchfile.tsv")
    params:
        batch = lambda wildcards: GENOMESBATCH[wildcards.batch],
    run: 
        with open(str(output),'w') as fh:
            for gid, file in params.batch.items():
                fh.write("{}\tUSER_{}\n".format(                    
                    file,
                    gid
                )) 
            