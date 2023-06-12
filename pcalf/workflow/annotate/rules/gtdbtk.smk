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



rule cm_empty_gtdbtk_summary:
    output:
        temp(os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv")),
    params:
        genomes = INPUT
    run:
        import pandas as pd
        df = pd.DataFrame(
            columns = GTDBSUMMARYCOLS[1:],
            index = list(params.genomes.keys())            
            )
        df.index.name = GTDBSUMMARYCOLS[0]        
        df.to_csv(str(output),sep="\t",header=True)
                
rule cm_target_quick_gtdbtk:
    output:
        touch(temp(os.path.join(RESDIR,"gtdbtk-res","gtdbtk.quick.done")))
    input:
        rules.cm_empty_gtdbtk_summary.output


rule gm_target_gtdbtk:
    output:
        touch(temp(os.path.join(RESDIR,"gtdbtk-res","gtdbtk.done")))
    input:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv"),
    
rule gm_GTDBTK_clean_table:
    output:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.clean.tsv")
    input:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.tsv")
    run:
        import os
        import pandas as pd
        if os.stat(str(input)).st_size != 0:
            df = pd.read_csv(str(input),sep="\t",header=0,index_col=None)        
            df.user_genome = df.user_genome.str.replace("USER_","")
            df.to_csv(str(output),sep="\t",header=True,index=False)
        open(str(output),'w').close()


rule gm_concatenate_archaea_and_bacteria_results:
    output:
        os.path.join(RESDIR , "gtdbtk-res","gtdbtk.ar53.bac120.summary.tsv")                                            
    input:
        bac = os.path.join(RESDIR , "gtdbtk-res","gtdbtk.bac120.summary.tsv"),
                        #"classify", 
        ar = os.path.join(RESDIR , "gtdbtk-res", "gtdbtk.ar53.summary.tsv"), 
                        #"classify", 
    shell:
        "cat {input} | grep -m 1 user_genome > {output} && " # grep header from input files, either input.bac and input.ar can be empty
        "tail -q -n +2 {input} >> {output} "
        # "cat {input.bac} "                # keep header from first input
        # "<(tail -q -n +2 {input.ar}) "    # remove header from other input
        # "> {output} "


rule gm_gtdbtk_classify_wf:
    output:
        os.path.join(RESDIR , "gtdbtk-res",  "gtdbtk.ar53.summary.tsv"),        
        os.path.join(RESDIR , "gtdbtk-res",  "gtdbtk.bac120.summary.tsv"),
    input:
        os.path.join(RESDIR , "gtdbtk-res", "batchfile.tsv")
    params:
        outdir = os.path.join(RESDIR , "gtdbtk-res"),
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


rule gm_gtdbtk_batchfile:
    output:
        os.path.join(RESDIR , "gtdbtk-res", "batchfile.tsv")
    params:
        input_datas = INPUT
    run: 
        with open(str(output),'w') as fh:
            for gid, files in params.input_datas.items():
                fh.write("{}\tUSER_{}\n".format(                    
                    files["genome"],
                    gid
                )) 
            