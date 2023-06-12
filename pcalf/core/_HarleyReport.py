

import json
import os
import re

import sqlite3
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import igraph
from igraph import Graph, EdgeSeq
from jinja2 import Template


class Report:
    NTER_CMAP = {
        "Z-type":"#7ad5a3",
        "X-type":"#322f26",
        "CoBaHMA-type":"#0babc1",
        "Y-type":"#b980d1",
        "Unknown-type":"#717171"
    }

    DOM_CMAP = {
        "Gly1":"#ecec13",   
        "Gly2":"#ec9d13",
        "Gly3":"#e25d3f",
        "GlyX3":"#d9d5d4"
    }
    DOM_CMAP.update(NTER_CMAP)

    def __init__(self):
        self.report = ""

    def connect_to_db(self,db):    
        return  sqlite3.connect(db)

    def make_trace(self,rid,record,y,**kwargs):
        traces = []
        segments = {}
        segment_id = 0
        for _,f in record["features"].items():
            domid = f["feature_id"]     
            if domid == "N-ter":
                domid = f["feature_src"].split("|")[0]

            x = list(range(f["feature_start"],f["feature_end"],20))
            hover = [domid for _ in x]
            segments[segment_id]= (x,hover)
            segment_id+=1


        # ADD FULL SEQ TRACES
        line = dict(color='black', width=1, dash='dash')
        xlen = [0,len(record['sequence'])] #list(range(0,len(record['sequence']),len(record['sequence'])))
        traces.append(go.Scatter(x=xlen, y=[y for i in range(0,len(xlen))],
                            mode='lines',
                            name="full-seq",
                            hovertext=rid,
                            line=line,#showlegend=False,                                             
                            legendgroup="full-seq",
                            legendgrouptitle_text="full-seq",
                            opacity = 0.4,
                            #fill="none",
                            #fillcolor=None,
                            #name="second legend group",                                                            
                            )
            ) 
                
        for _,segm in segments.items():
                if segm:
                    if segm[0][0] and segm[0][1] :
                        sx,sh = segm
                        sy = [y for i in range(0,len(sx))]
                        dom = sh[0]
                        if dom in self.DOM_CMAP:                
                            line = dict(color=self.DOM_CMAP[dom],width=3)
                        else:
                            #print(dom)
                            line = dict(color='black', width=1)

                        neighbor = "NA"
                        if "nter_neighbor" in record:
                            neighbor = record["nter_neighbor"]
                        htxt = "- {}<br>- {}<br>- {} [neighbor]<br>- {}<br>- {}<br>".format(rid,dom , neighbor , record["flag"] ,"From {} to {}".format(min(sx),max(sx)))
                        for i,j in kwargs.items():
                            htxt+="- {}: {}<br>".format(i,j)

                        traces.append(go.Scatter(
                            x=sx, y=sy,
                            mode='lines',
                            name=dom,
                            line=line,#showlegend=False,                                             
                            legendgroup=dom,
                            legendgrouptitle_text=dom,
                            hovertext=htxt,
                            #opacity = 0.4,
                            #fill="none",
                            #fillcolor=None,
                            #name="second legend group",                                                            
                            )
                        )    
        return traces
        
    def make_report(self,db,templatedir):
        dbname = os.path.basename(db).split(".")[0]
        cnx = self.connect_to_db(db)
        summary_df = pd.read_sql_query("SELECT * FROM genomes", cnx,index_col="Assembly Accession")

        # Reformat organism name
        summary_df["Organism Name"] = summary_df.apply(lambda x: re.sub(r'[^a-zA-Z0-9\s]+', '', x["Organism Name"]),axis=1)

        # Fix date 
        #summary_df["Download date"] = summary_df.apply(lambda x : "2022-07-09" if x.name in d2022_07_09 else x["Download date"] , axis = 1)

        # We count the number of genome with or without ccyA depending on their NBCI state (GCA,GCF, GCA and GCF or strain)
        genomes_metrics = {
            "Strain":{
                "ccyA+":0,
                "ccyA-":0,        
                "strains":[]
            },
            "GenBanK-RefSeq":{
                "ccyA+":0,"ccyA-":0,"ccyA~":0,"genomes":[],"l_ccyA~":[],
            }, 
            "GenBanK":{
                "ccyA+":0,"ccyA-":0,"genomes":[]
            }, 
            "RefSeq":{
                "ccyA+":0,"ccyA-":0,"genomes":[]
            }
        }
        for gid , row in summary_df[["Organism Name","flag"]].iterrows():
            is_calc = "Calcyanin" in str(row.flag)
            if row["Organism Name"] not in genomes_metrics["Strain"]["strains"]:
                genomes_metrics["Strain"]["strains"].append(row["Organism Name"])
                if is_calc:
                    genomes_metrics["Strain"]["ccyA+"]+=1
                else:
                    genomes_metrics["Strain"]["ccyA-"]+=1
            db="RefSeq"
            cr = gid.replace("GCF","GCA") 
            if gid.startswith("GCA"):
                db = "GenBanK"
                cr = gid.replace("GCA","GCF")
            
            if gid not in genomes_metrics[db]["genomes"]:
                genomes_metrics[db]["genomes"].append(gid)
                if is_calc:
                    genomes_metrics[db]["ccyA+"]+=1
                else:
                    genomes_metrics[db]["ccyA-"]+=1
            
            if cr in list(summary_df.index):
                consensus_acc = gid.split("_")[1].split(".")[0]
                if consensus_acc not in genomes_metrics["GenBanK-RefSeq"]["genomes"]:
                    genomes_metrics["GenBanK-RefSeq"]["genomes"].append(consensus_acc)
                    # is_cr_calc : does the cross reference has a calcyanin ?
                    is_cr_calc = "Calcyanin" in str(summary_df.at[cr,"flag"])
                    if is_calc and is_cr_calc:
                        genomes_metrics["GenBanK-RefSeq"]["ccyA+"]+=1
                    elif (is_calc and not is_cr_calc) or (not is_calc and is_cr_calc):
                        genomes_metrics["GenBanK-RefSeq"]["ccyA~"]+=1
                        genomes_metrics["GenBanK-RefSeq"]["l_ccyA~"].append(consensus_acc)
                    else:
                        genomes_metrics["GenBanK-RefSeq"]["ccyA-"]+=1
                    
        # We reformat the dictionnary and make a plotly-friendly dataframe.
        metrics_d = []
        for i,j in genomes_metrics.items():
            metrics_d.append(
                (i,"ccyA+",j["ccyA+"])
            )
            metrics_d.append(
                (i,"ccyA-",j["ccyA-"])
            )
            metrics_d.append(
                (i,"ccyA~",j["ccyA~"] if 'ccyA~' in j else 0)
            )

        # We make a figure with multiple pie-chart [strain, GCA-GCF, GCA, GCF].   
        metrics_df = pd.DataFrame(metrics_d)
        metrics_df.columns = ["is","genotype","count"]
        metrics_fig = px.pie(metrics_df,values="count",names="genotype",facet_col="is",color="genotype",
                color_discrete_map={"ccyA-":"lightgrey","ccyA+":"lightgreen","ccyA~":"orange"})
    
        # We count the increasing number of genomes over the time.
        total_genomes_count = pd.DataFrame(summary_df.reset_index().groupby("Download date").count()["Assembly Accession"])
        total_genomes_count.columns=["count_by_date"]
        total_genomes_count["total_genomes"] = total_genomes_count["count_by_date"].cumsum(axis = 0)
        genome_over_time = px.line(total_genomes_count.reset_index(),x="Download date", y="total_genomes",color_discrete_sequence=["mediumvioletred"],markers=True,title="Number of genomes over time")   

        # we load calcyanin and ccyA related datas into a dataframe + formating
        #if utils.check_if_table_exist(db,"pcalf_summary") and utils.check_if_table_exist(db,"ccya"):
        pcalf_df = pd.read_sql_query("""SELECT T1.sequence_id,   
                                                T1.sequence_src, T1.flag, T1.nter, T1.nter_neighbor, T1.cter, T1.sequence, T1.iteration,
                                                T2.ccyA_genomic_region, T2.ccyA_start, T2.ccyA_stop, T2.ccyA_frame, T2.ccyA_partial, T2.ccyA_pseudo, T2.ccyA_seq                                                              
                                    FROM pcalf_summary AS T1 LEFT JOIN ccya AS T2 ON T1.sequence_id = T2.sequence_id""", cnx )    
        pcalf_df["Download date"] = pcalf_df.sequence_src.map(summary_df["Download date"].to_dict())
        pcalf_df["Organism Name"] = pcalf_df.sequence_src.map(summary_df["Organism Name"].to_dict())    
        # First we check that the dataframe is not empty.
        if not pcalf_df.empty:
            pcalf_df["Organism Name"] = pcalf_df.apply(lambda x: re.sub(r'[^a-zA-Z0-9\s]+', '', x["Organism Name"]),axis=1)
            pcalf_df["iteration"] = pcalf_df.apply(lambda x: "Iteration_{}".format(x["iteration"]),axis=1)        
            pcalf_df.nter.fillna("Unknown-type",inplace=True)
            #pcalf_df = pcalf_df.T.drop_duplicates(keep="first").T#   we drop sequence_id column as its duplicated after the join
        
            # We group sequences based on several attributes and we count the number of sequences for each group 
            pcalf_cpt_df = pd.DataFrame(pcalf_df.groupby(["Download date","flag","nter","cter","iteration"]).count().sequence_id).reset_index()  
            pcalf_cpt_df.columns=["Download date","flag","nter","cter","iteration","count"]
            
            # We make a sunburst chart : 
            sunburst = px.sunburst(pcalf_cpt_df, path=['nter', 'flag', 'cter', "Download date", 'iteration'], values='count', color = "nter",
                        color_discrete_map=self.NTER_CMAP)
            sunburst.update_layout(title="No of Calcyanin.")    
            # And another one with the same kind of information : a treemap:
            treemap = px.treemap(pcalf_cpt_df, path=[px.Constant(""), 'nter', 'flag', 'cter', 'Download date'], values='count',
                            color='nter',color_discrete_map=self.NTER_CMAP)
            treemap.update_traces(root_color="lightgrey")
            treemap.update_layout(title="No of Calcyanin.", margin = dict(t=50, l=25, r=25, b=25))
            del(pcalf_cpt_df)
            
            # We retrieve the number of calcyanin by N-ter for each date.        
            pcalf_cpt_df = pd.DataFrame(pcalf_df.groupby(["Download date","nter",'flag']).count().sequence_id).reset_index()
            pcalf_cpt_df.columns=["Download date","nter",'flag','count']
            # We get the cumulative count for each category.
            cnt = {}
            cnt_data = []
            idxs = []
            for _ , row in pcalf_cpt_df.iterrows():
                idx = ";".join([row.nter, row.flag])
                if idx not in idxs:
                    idxs.append(idx)
                d = row["Download date"]
                if  d not in cnt:
                    cnt[d]={}
                if idx not in cnt[d]:
                    cnt[d][idx]=0
                cnt[d][idx]+=row["count"]
                cnt_data.append(cnt[d][idx])

                # we create null node when no calcyanin was found at a specific date.    
                for d in list(summary_df["Download date"].unique()):
                    if d not in cnt:
                        cnt[d]={}
                    for idx in idxs:                    
                        if idx not in cnt[d]:
                            cnt[d][idx]=0
                # and we build a new plotly-friendly dataframe.
                cpt_d = []
                total = {}    
                for _,d in cnt.items():
                    for k,v in d.items():
                        if k not in total:
                            total[k]=0
                        total[k]+=v
                        cpt_d.append({"Download date":_,"nter":k.split(";")[0],"flag":k.split(";")[1],"count":v, "Total count":total[k] })
            timeseries_df = pd.DataFrame(cpt_d)    
            # We then make a grouped line plot.
            nter_over_time = px.line(timeseries_df,x="Download date", y="Total count",
                        color="nter", 
                        symbol="flag",
                        markers=True,title="Number of calcyanin over time")
        else:
            # There is no calcyanin - so we can't make those charts.
            nter_over_time = None
            sunburst = None
            treemap = None    
        
        figures = []
        if nter_over_time:
            figures.append(nter_over_time)
        figures.append(genome_over_time)
        timeseries = make_subplots(rows=len(figures), cols=1) 
        for i, figure in enumerate(figures):
            for trace in range(len(figure["data"])):
                timeseries.append_trace(figure["data"][trace], row=i+1, col=1)
        if len(figures) == 2:
            title="No of calcyanin [top] and No of assembly [bottom] over time."
        else:
            title="No of assembly over time."
        timeseries.update_layout(title=title)

        # Plotting modular organization
        # First we load the feature 'table'
        features_df = pd.read_sql_query("SELECT * FROM pcalf_features", cnx ).set_index("sequence_id")
        features_df["e-value"] = features_df.apply(lambda x : '{:0.2e}'.format(x["e-value"]),axis=1)
        features_dict = {}
        # dictionnary with sequence identifiers as keys and features as values
        for _ , sdf in features_df.groupby("sequence_id"):
            features_dict[_]=sdf.reset_index().T.to_dict()
        features_dict
        # dictionnary with sequence identifiers as keys and sequence-related information + features as values
        sequence_datas = {}
        for _ , row in pcalf_df.set_index("sequence_id").iterrows():
            sequence_datas[_] = row.to_dict()
            sequence_datas[_]["features"]=features_dict[_]  
        
        # For each type of N-ter we make a modular organization chart.    
        oms_plots = {}
        for nter_type in list(set(pcalf_df.nter)):
            fig = go.Figure()
            y = 1
            cpt = 0
            for rid , record in sequence_datas.items():
                if record["nter"] == nter_type:
                    cpt+=1
                    traces = self.make_trace(rid,
                                        record,
                                        y,
                                        Organism_Name=record["Organism Name"],
                                        Assembly=record["sequence_src"]
                                    )
                    if traces:
                        for t in traces:
                            fig.add_trace(t)
                        y+=1
                    cpt+=1
                    leg = []
                    for trace in fig['data']:
                        leg.append
                        if trace['name'] not in leg:
                            leg.append(trace['name'])
                        else:
                            trace['showlegend'] = False
                fig.update_layout(
                    autosize=False,
                    width=1000,
                    height=10*cpt if cpt > 100 else 700,
                    margin=dict(
                        l=50,
                        r=50,
                        b=100,
                        t=100,
                        pad=4
                    ),
                )   
                oms = fig #.to_html().split("<body>")[1].split("</body>")[0]  
                oms_plots[nter_type] = oms
        
        # DECISION TREE CHART
        nr_vertices = 11
        v_label = list(map(str, range(nr_vertices)))
        G = Graph.Tree(nr_vertices, 2) # 2 stands for children number
        lay = G.layout_reingold_tilford(mode="in", root=[0])
        position = {k: lay[k] for k in range(nr_vertices)}
        Y = [lay[k][1] for k in range(nr_vertices)]
        M = max(Y)
        es = EdgeSeq(G) # sequence of edges
        E = [e.tuple for e in G.es] # list of edges
        L = len(position)
        Xn = [position[k][0] for k in range(L)]
        Yn = [2*M-position[k][1] for k in range(L)]
        Xe = []
        Ye = []
        Ec = {}
        color = "#c6587e"
        txt = "N"
        for edge in E:
            Xe=[position[edge[0]][0],position[edge[1]][0], None]
            Ye=[2*M-position[edge[0]][1],2*M-position[edge[1]][1], None]
            Ec[edge] = [Xe,Ye,color,txt]
            # we alternate the vertice color for positive and negative answer
            if color == "#48d38b":
                color = "#c6587e"
                txt = "N"
            else:
                color = "#48d38b"
                txt = "Y"
        # Nodes values
        node_txt = {
            0:"[Sequence has a significative hit against the GlyX3] Does the sequence have three glycine zipper in the right order (G1|G2|G3) ?",
            1:"Does the sequence have<br>at least a G1 and G3 in this order ?",
            2:"Does the sequence have<br>a Known N-ter ?",
            3:"Does the sequence have<br>Known N-ter ?",
            4:"Does the sequence have<br>a N-ter of type Y ?",
            5:"Calcyanin with<br>new N-ter",
            6:"Calcyanin with<br>known N-ter",
            7:"Atypical gly region<br>with new N-ter",
            8:"Atypical gly region<br>with known N-ter",
            9:"Calcyanin with<br>new N-ter",
            10:"Calcyanin with<br>known N-ter",
        }  

        labels = ["<b>"+txt+"</b>" if _ in [9,10,5,6] else txt for _,txt in node_txt.items()]
        decision_tree = go.Figure()
        for e,E in Ec.items():
            decision_tree.add_trace(go.Scatter(x=E[0],
                            y=E[1],
                        mode='lines',
                        line=dict(color=E[2], width=4),
                        text = E[3],
                        hoverinfo="text",
        #                   hoverinfo='none'
                        ))


        decision_tree.add_trace(go.Scatter(x=Xn,
                        y=Yn,
                        mode='markers+text',
                        name='bla',
                        marker=dict(symbol='diamond',
                                        size=50,
                                        color='black',    #'#DB4551',
                                        line=dict(color='rgb(50,50,50)', width=1)
                                        ),
                        text=labels,
                        textposition='top center',
                        hoverinfo='text',
                        opacity=1
                        ))
        decision_tree.update_xaxes(visible=False)
        decision_tree.update_yaxes(visible=False)
        decision_tree.update_layout({
        'margin':dict(l=40, r=40, b=85, t=100),
        'showlegend':False,
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        'title': "Calcyanin classification decision tree.",
        })

        # We convert the "strain" dataframe into a dictionnary - later we will inject those datas into the report file using jinja2.
        DATAS = {}
        #keys are the name of the organism, value are dictionnaries (gid : datas)
        for _ , sdf in summary_df.groupby("Organism Name"):
            DATAS[_]=sdf.T.to_dict()

        # We add calcyanin related data for each strain-genome.
        for seqid, seqrecord in sequence_datas.items():
            org = seqrecord["Organism Name"]
            src = seqrecord["sequence_src"]
            if "sequences" not in DATAS[org][src]:
                DATAS[org][src]["sequences"]={}
            DATAS[org][src]["sequences"][seqid]=seqrecord
        # we convert the dictionnary into a json object for jinja2 injection.
        json_object = json.dumps(DATAS, indent = 4)     

        # And we fill the template 
        with open(os.path.join(templatedir,'template.html')) as file_:
            template = Template(file_.read())

        report = template.render(
            dbname = dbname,
            datas=json_object,
            css= open( os.path.join(templatedir,"template.css" )).read(),
            js = open( os.path.join(templatedir,"template.js"  )).read(),
            decision_tree = decision_tree.to_html().split("<body>")[1].split("</body>")[0],
            sunburst = sunburst.to_html().split("<body>")[1].split("</body>")[0] if sunburst else "<p style='font-style: italic;'>No calcyanin detected - No sunburst :/ </p>"  ,
            treemap = treemap.to_html().split("<body>")[1].split("</body>")[0] if sunburst else "<p style='font-style: italic;'>No calcyanin detected - No treemap :/ </p>"  ,
            cobahma_oms = oms_plots["CoBaHMA-type"].to_html().split("<body>")[1].split("</body>")[0] if "CoBaHMA-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,
            x_oms = oms_plots["X-type"].to_html().split("<body>")[1].split("</body>")[0] if "X-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,
            y_oms = oms_plots["Y-type"].to_html().split("<body>")[1].split("</body>")[0] if "Y-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,
            z_oms = oms_plots["Z-type"].to_html().split("<body>")[1].split("</body>")[0] if "Z-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,
            unknown_oms = oms_plots["Unknown-type"].to_html().split("<body>")[1].split("</body>")[0]  if "Unknown-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,
            timeserie = timeseries.to_html().split("<body>")[1].split("</body>")[0],
            metrics_fig = metrics_fig.to_html(config= dict(displayModeBar = False)).split("<body>")[1].split("</body>")[0],
            strain_total = len(genomes_metrics["Strain"]["strains"]),    
            cr_total =  len(genomes_metrics["GenBanK-RefSeq"]["genomes"]) , 
            gb_total =  len(genomes_metrics["GenBanK"]["genomes"]) , 
            rs_total =  len(genomes_metrics["RefSeq"]["genomes"]) , 
            total = len(genomes_metrics["RefSeq"]["genomes"]) + len(genomes_metrics["GenBanK"]["genomes"]),
        )
        self.report = report

    def save(self,file):
        dirname = os.path.abspath(os.path.dirname(file))
        os.makedirs(dirname,exist_ok=True)
        with open( file , 'w' ) as stream:
            stream.write(self.report)


    if __name__ == "__main__":
        import sys
        import shutil
        dbfile = sys.argv[1]
        templatedir = sys.argv[2]
        outfile = sys.argv[3]
        report_output = make_report(
            dbfile,
            templatedir
        )
        with open( outfile , 'w' ) as stream:
            stream.write(report_output)
        
        if shutil.which("open"):
            os.system("open {}".format( 
                outfile
            ))