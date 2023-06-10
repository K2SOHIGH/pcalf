#!/usr/bin/env python
# coding: utf-8

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

NTER_CMAP = {
        "Z-type":"#7ad5a3",
        "X-type":"#322f26",
        "CoBaHMA-type":"#0babc1",
        "Y-type":"#b980d1",
        "Unknown-type":"#717171"
}

def define_ccya_genotype(x):
    if x.Accession.startswith("GCA"):
        prefix="GCA" 
    elif x.Accession.startswith("GCF"):
        prefix = "GCF"
    else:
        return None
    
    if x.flag == "Calcyanin with known N-ter" or x.flag == "Calcyanin with new N-ter":
        return prefix + "+"
    return prefix+"-"

def count_genotypes(df , group):
    datas = {"ccyA+":0,"ccyA-":0}
    if group in df.columns:
        df = df.reset_index().set_index(group)
        for k in list(set(df.index)):
            #datas[k]="ccyA-"
            loc = df.loc[k]
            if isinstance(loc,pd.Series):
                loc = pd.DataFrame(loc).T
            if "ccyA+" in loc["(ccyA)"].values:                
                datas["ccyA+"]+=1
                continue
            datas["ccyA-"]+=1
        return datas
    else:
        raise KeyError(group)
        


def make_genome_pie_chart(cnx):
    summary_df = pd.read_sql_query("""
        SELECT * FROM genomes as g JOIN 
        harley as h on h.Accession=g.Accession LEFT JOIN 
        summary as s on G.Accession = s.sequence_src""", 
                                   cnx , index_col="Accession").reset_index()
    #return summary_df
    summary_df.Accession = summary_df.apply(lambda x : x.Accession[0],axis=1)
    summary_df["uid"] = summary_df.apply(lambda x : x.Accession.split("_")[-1].split(".")[0], axis=1 )
    summary_df["uidv"] = summary_df.apply(lambda x : x.Accession.split("_")[-1], axis=1 )
    summary_df["(DBccyA)"] = summary_df.apply(lambda x : define_ccya_genotype(x), axis=1 )
    summary_df["(ccyA)"] = summary_df.apply(lambda x : "ccyA+" if x.flag in ["Calcyanin with known N-ter","Calcyanin with new N-ter"] else "ccyA-", axis=1 )
    summary_df = summary_df[["Organism", "uid", "uidv", "Accession", "Assembly name" , "Date" , "(ccyA)" , "sequence_accession", "flag", "nter"]]

    # We reformat the dictionnary and make a plotly-friendly dataframe.
    metrics_d = []
    for g in ["Organism", "uid" , "uidv" , "Accession"] :
        d = count_genotypes(summary_df , group = g)      
        total = len(list(summary_df[g].unique()))
        metrics_d.append(
            ("{} [{}]".format(g,total) ,"ccyA+",d["ccyA+"])
        )
        metrics_d.append(
            ("{} [{}]".format(g,total),"ccyA-",d["ccyA-"])
        )


    # We make a figure with multiple pie-chart [strain, GCA-GCF, GCA, GCF].   
    metrics_df = pd.DataFrame(metrics_d)
    metrics_df.columns = ["RedLevel","genotype","count"]
    fig = px.pie(metrics_df,values="count",names="genotype",facet_col="RedLevel",color="genotype",
            color_discrete_map=
#                 {"ccyA-":"#ffb4b4","ccyA+":"#439775","ccyA~":"orange"}
                 {"ccyA-":"#D5D5D8","ccyA+":"#88D9E6","ccyA~":"orange"}
                 
                )
    return fig

#make_genome_pie_chart(cnx)

def make_decision_tree_chart():
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
            0:"[Sequence has a significative hit against the GlyX3]<br>Does the sequence have three glycine zipper in the right order (G1|G2|G3) ?",
            1:"Does the sequence have<br>at least a G1 and G3 in this order ?",
            2:"Does the sequence have<br>a Known N-ter ?",
            3:"Does the sequence have<br>Known N-ter ?",
            4:"Does the sequence have<br>a N-ter of type Y ?",
            5:"Calcyanin with<br>new N-ter",
            6:"Calcyanin with<br>known N-ter",
            7:"Atypical gly region<br>with new N-ter",
            8:"Atypical gly region<br>with known N-ter",
            9:"Atypical gly region<br>with new N-ter",#"Calcyanin with<br>new N-ter",
            10:"Calcyanin with<br>known N-ter",
        }  

        labels = ["<b>"+txt+"</b>" if _ in [9,10,5,6] else txt for _,txt in node_txt.items()]
        decision_tree = go.Figure()
        for e,E in Ec.items():
            decision_tree.add_trace(go.Scatter(x=E[0],
                            y=E[1],
                        mode='lines',
                        line=dict(color=E[2], width=4),
                        opacity=0.5,
                        text = E[3],
                        hoverinfo="text",
        #                   hoverinfo='none'
                        ))


        decision_tree.add_trace(go.Scatter(x=Xn,
                        y=Yn,
                        mode='markers+text',
                        name='bla',
                        marker=dict(symbol='diamond',
                                        size=30,
                                        color='white',    #'#DB4551',
                                        #opacity=0.5,
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
        'height':750,
        })
        return decision_tree
#make_decision_tree_chart()    

# # Number of cyanobacteria over time
def count_genome_by_date(df , group , label):
    df.sort_values("Date",inplace=True)
    if group in df.columns:
        df = df.reset_index().set_index(group)
        df = df[~df.index.duplicated(keep="first")]
        df.reset_index(inplace=True)
        df = pd.DataFrame(df.groupby("Date").count()[group])
        df.columns=["count_by_date"]
        df["total"] = df["count_by_date"].cumsum(axis = 0)
        df["label"] = label
        return df
    else:
        raise KeyError(group)

def make_genome_over_time_chart(cnx):
    summary_df = pd.read_sql_query("""
        SELECT * FROM genomes as g JOIN 
        harley as h on h.Accession=g.Accession LEFT JOIN 
        summary as s on G.Accession = s.sequence_src""", 
                                   cnx , index_col="Accession").reset_index()
    summary_df.Accession = summary_df.apply(lambda x : x.Accession[0],axis=1)
    summary_df["uid"] = summary_df.apply(lambda x : x.Accession.split("_")[-1].split(".")[0], axis=1 )
    summary_df["uidv"] = summary_df.apply(lambda x : x.Accession.split("_")[-1], axis=1 )
    summary_df["(DBccyA)"] = summary_df.apply(lambda x : define_ccya_genotype(x), axis=1 )
    summary_df["(ccyA)"] = summary_df.apply(lambda x : "ccyA+" if x.flag in ["Calcyanin with known N-ter","Calcyanin with new N-ter"] else "ccyA-", axis=1 )
    summary_df = summary_df[["Organism", "uid", "uidv", "Accession", "Assembly name" , "Date" , "(ccyA)" , "sequence_accession", "flag", "nter"]]


    assembly_count_df = pd.concat(
        [count_genome_by_date(summary_df,"Organism","Strain<br>[i.e Microcystis aeruginosa PCC 9443]"),
         count_genome_by_date(summary_df,"uid","Assembly<br>[i.e XXX_<assembly>.N]"),
         count_genome_by_date(summary_df,"uidv","Version<br>[i.e XXX_<assembly>.<version>]"),     
         count_genome_by_date(summary_df,"Accession","Entry<br>[i.e <ncbi db>_<assembly>.<version>]"),
        ])

    genome_over_time = px.line(assembly_count_df.reset_index().sort_values(["Date","label"]),
                               hover_data=["count_by_date"],
                               x="Date", y="total",color="label" ,markers=True)   
    genome_over_time.update_layout(
        title="Number of entry over time",
        yaxis_title="#entries",
        xaxis_title="Date",
        legend_title="Level of redundancy",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False),
        yaxis=dict(gridcolor="rgba(0,0,0,0.1)")
    )
    return genome_over_time

# # Number of sequences over time
def count_seq_by_date(df , group , subset = [] , label = None):
    df.sort_values("Date",inplace=True)
    if group in df.columns:
        df = df.reset_index().set_index(group)
        df = df[~df.index.duplicated(keep="first")]

        df.reset_index(inplace=True)
        df.fillna("NA",inplace=True)
        grouping_variables = ["Date"]+subset

        df = pd.DataFrame(df.groupby(grouping_variables).count()[group])
        df.columns=["count_by_date"]
        return df
    else:
        raise KeyError(group)

def make_sequence_over_time_chart(cnx):
    summary_df = pd.read_sql_query("""
        SELECT * FROM genomes as g JOIN 
        harley as h on h.Accession=g.Accession JOIN 
        summary as s on G.Accession = s.sequence_src""", 
                                   cnx , index_col="Accession").reset_index()
    summary_df.Accession = summary_df.apply(lambda x : x.Accession[0],axis=1)

    assembly_count_df = count_seq_by_date(summary_df,"sequence_accession", ["flag","nter"])


    datas = []
    for flag,nter in assembly_count_df.reset_index("Date").index.unique():
        tmp = assembly_count_df.reset_index()
        tmp = tmp[(tmp.flag == flag) & (tmp.nter == nter)]
        tmp["total"] = tmp.count_by_date.cumsum(axis = 0)
        datas.append(tmp)
    assembly_count_df = pd.concat(datas)


    fig = px.line(assembly_count_df.reset_index(),x="Date", y="total",
                               color_discrete_map=NTER_CMAP,
                               hover_data=["count_by_date"],color="nter",line_dash="flag" ,markers=True)   
    fig.update_layout(
        title="Number of calcyanin by N-ter type and flag over time",
        yaxis_title="#sequences",
        xaxis_title="Date",
        legend_title="Calcyanin N-ter type and flag",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False),
        yaxis=dict(gridcolor="rgba(0,0,0,0.1)")
    )
    return fig



# # Modular organization
DOM_CMAP = {
        "Gly1":"#ecec13",   
        "Gly2":"#ec9d13",
        "Gly3":"#e25d3f",
        "GlyX3":"#d9d5d4"
    }
DOM_CMAP.update(NTER_CMAP)




def sequence_modular_orga(rid,record,y,**kwargs):
        """
        rid : Sequence accession
        record : 
        """
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
                        if dom in DOM_CMAP:                
                            line = dict(color=DOM_CMAP[dom],width=3)
                        else:
                            #print(dom)
                            line = dict(color='black', width=1)

                        neighbor = "NA"
                        if "nter_neighbor" in record:
                            neighbor = record["nter_neighbor"]
                        htxt = "- {}<br>- {}<br>- {}<br>- {}<br>- {} [nearest neighbor]<br>".format(
                            rid,
                            dom, 
                            record["flag"],
                            "From {} to {}".format(min(sx),max(sx)),
                            neighbor
                        )
                        customdatas = [{"seqid":rid}]
                        for i,j in kwargs.items():
                            htxt+="- {}: {}<br>".format(i,j)
                            customdatas.append({i:j})

                        traces.append(go.Scatter(
                                x=sx, y=sy,
                                mode='lines',
                                name=dom,
                                line=line,#showlegend=False,                                             
                                legendgroup=dom,
                                legendgrouptitle_text=dom,
                                text="",
                                hoverlabel=None,
                                hovertext=htxt,
                                hoveron='points+fills',
                                customdata = customdatas,                                                            
                            )
                        )    
        return traces


def make_modorg_chart(cnx):
    # Plotting modular organization
    # First we load the feature 'table'
    features_df = pd.read_sql_query("SELECT * FROM features as f JOIN summary as s on s.sequence_accession=f.sequence_id", cnx ).set_index("sequence_id")
    features_df["e-value"] = features_df.apply(lambda x : '{:0.2e}'.format(x["e-value"]),axis=1)
    features_dict = {}
    # dictionnary with sequence identifiers as keys and features as values
    for _ , sdf in features_df.groupby("sequence_id"):
        features_dict[_]=sdf.reset_index().T.to_dict()
    features_dict
    # dictionnary with sequence identifiers as keys and sequence-related information + features as values
    sequence_datas = {}

    sequence_df = pd.read_sql_query("""
        SELECT * FROM summary as s JOIN genomes as g on g.Accession = s.sequence_src
        """, cnx ).set_index("sequence_accession")
    for _ , row in sequence_df.iterrows():
        sequence_datas[_] = row.to_dict()
        sequence_datas[_]["features"]=features_dict[_]  

    # For each type of N-ter we make a modular organization chart.    
    oms_plots = {}
    for nter_type in list(set(sequence_df.nter)):
        if not nter_type:
            nter_type = "Unknown-type"
        fig = go.Figure()
        y = 1
        cpt = 0
        for rid , record in sequence_datas.items():
            if not record["nter"]:
                record["nter"] = "Unknown-type"
            if record["nter"] == nter_type:
                cpt+=1
                traces = sequence_modular_orga(rid,
                                    record,
                                    y,
                                    Organism_Name=record["Organism"],
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
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                yaxis=dict(visible=False),
                xaxis=dict(gridcolor="rgba(0,0,0,0.1)"),
                xaxis_title="Sequence length (aa)",
                title="Modular organization [{}]".format(nter_type)
            )   
            oms = fig #.to_html().split("<body>")[1].split("</body>")[0]      
            oms_plots[nter_type] = oms
    return oms_plots



def make_sunburst(cnx):
    sequence_df = pd.read_sql_query("""
        SELECT * FROM summary as s JOIN 
        genomes as g on g.Accession = s.sequence_src JOIN
        harley as h on h.Accession = g.Accession

        """,cnx).set_index("sequence_accession")
    sequence_df.nter.fillna("Unknown-type",inplace=True)
    sequence_df = sequence_df.groupby(
        ["Date","flag","nter","cter"]).count().reset_index() 
    sequence_df = sequence_df[["Date","flag","nter","cter","sequence_src"]]
    sequence_df.columns = ["Date","flag","nter","cter","No of sequences"]
    # We make a sunburst chart : 
    sunburst = px.sunburst(sequence_df, path=['nter', 'flag', 'cter', "Date"], 
                           values='No of sequences', 
                           color = "nter",
                           color_discrete_map=NTER_CMAP)
    sunburst.update_layout(title="No of Calcyanin.")
    return sunburst

def make_calcyanin_treemap(cnx):
    sequence_df = pd.read_sql_query("""
        SELECT * FROM summary as s JOIN 
        genomes as g on g.Accession = s.sequence_src JOIN
        harley as h on h.Accession = g.Accession

        """,cnx).set_index("sequence_accession")
    sequence_df.nter.fillna("Unknown-type",inplace=True)
    sequence_df = sequence_df[["flag","nter","cter","sequence_src"]]
    sequence_df.columns = ["flag","nter","cter","No of sequences"]
    sequence_df = sequence_df.groupby(
        ["flag","nter","cter"]).count().reset_index() 
    # We make a sunburst chart : 
    # And another one with the same kind of information : a treemap:
    treemap = px.treemap(sequence_df, path=['nter', 'flag', 'cter'], values='No of sequences',
                    color='nter',color_discrete_map = NTER_CMAP)
    treemap.update_layout(title="No of Calcyanin.")
    return treemap

#make_calcyanin_treemap(cnx)


def make_data(cnx):
    # We convert the "strain" dataframe into a dictionnary - later  we will inject those datas into the report file using jinja2.
    DATAS = {}

    genomes_df = pd.read_sql_query("""
        SELECT * FROM harley as h  JOIN 
        genomes as g on g.Accession = h.Accession LEFT JOIN 
        checkm as c on g.Accession=c.`Bin Id` LEFT JOIN
        gtdbtk as t on g.Accession=t.`user_genome`
        """,cnx)

    genomes_df = genomes_df.fillna("NA")
    genomes_df = genomes_df.loc[:,~genomes_df.columns.duplicated()]
    genomes_df = genomes_df[['Accession', 'Assembly name','Date', 'Submitter',
       'Submission date', 'Isolate', 'TaxID', 'Organism', 'Biosample',
       'Isolation source', 'Environment (biome)', 'Geographic location',
       'Culture collection', 'Collection date', 'Sample type',
       'Completeness', 'Contamination', 'Strain heterogeneity',
       'Genome size (bp)', '# scaffolds', '# contigs',
       'N50 (scaffolds)', 'N50 (contigs)','# predicted genes',
       'classification', 'fastani_reference','fastani_ani']]
    
    ccyA_plus = []
    ccyA_minus = []
    for org , sdf in genomes_df.groupby(["Organism"]):
        sdf.index = sdf.Accession
        DATAS[org]=sdf.T.to_dict()
        for acc in sdf.index:
            seqs = pd.read_sql_query("""
                SELECT * FROM summary as s JOIN 
                ccya as c on c.sequence_id=s.sequence_accession WHERE
                s.sequence_src="{}"
                """.format(acc),cnx,index_col="sequence_id")
            seqs = seqs.fillna("NA")
            seqs = seqs.T.to_dict()            
            DATAS[org][acc]["sequences"] = seqs
            
            if seqs:
                ccyA_plus.append(org)
            else:
                ccyA_minus.append(org)
            
            for seq in seqs.keys():
                features = pd.read_sql_query("""
                    SELECT * FROM features as f WHERE
                        f.sequence_id="{}"
                """.format(seq),cnx)                
                features = features.fillna("NA")
                features.sort_values(["feature_id","e-value"],inplace=True)
                DATAS[org][acc]["sequences"][seq]["features"] = features.T.to_dict()
                hits = pd.read_sql_query("""
                    SELECT * FROM hits as f WHERE
                        f.sequence_id="{}"
                """.format(seq),cnx)
                hits = hits.fillna("NA")
                hits.sort_values(["hit_src","hit_e_value"],inplace=True)
                DATAS[org][acc]["sequences"][seq]["hits"] = hits.T.to_dict()
    
    DATAS = {i:DATAS[i] for i in list(set(ccyA_plus + ccyA_minus))}
                
                
    # we convert the dictionnary into a json object for jinja2 injection.
    json_object = json.dumps(DATAS, indent = 4)     
    return json_object , DATAS


# And we fill the template 
def save(file,report):
        dirname = os.path.abspath(os.path.dirname(file))
        os.makedirs(dirname,exist_ok=True)
        with open( file , 'w' ) as stream:
            stream.write(report)

            
            
def render(db,templatedir,outfile):
    with open(os.path.join(templatedir,"pcalf.svg"),"r") as f:
        workflow = f.read().rstrip()

    # dbname = os.path.basename(db).split(".")[0]
    cnx = sqlite3.connect(db)
    with open(os.path.join(templatedir,'template.html')) as file_:
        template = Template(file_.read())

    oms_plots = make_modorg_chart(cnx)

    report = template.render(
        
        datas  = make_data(cnx)[0],

        css= open( os.path.join(templatedir,"template.css" )).read(),

        js = open( os.path.join(templatedir,"template.js"  )).read(),
        
        workflow = workflow,

        decision_tree = make_decision_tree_chart().to_html().split("<body>")[1].split("</body>")[0],

        sunburst = make_sunburst(cnx).to_html().split("<body>")[1].split("</body>")[0],# if sunburst else "<p style='font-style: italic;'>No calcyanin detected - No sunburst :/ </p>"  ,

        treemap = make_calcyanin_treemap(cnx).to_html().split("<body>")[1].split("</body>")[0],# if sunburst else "<p style='font-style: italic;'>No calcyanin detected - No treemap :/ </p>"  ,

        cobahma_oms = oms_plots["CoBaHMA-type"].to_html(
            full_html=False,div_id="cobahma-type-plot",include_plotlyjs=False) if "CoBaHMA-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,

        x_oms = oms_plots["X-type"].to_html(
            full_html=False,div_id="x-type-plot",include_plotlyjs=False) if "X-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,

        y_oms = oms_plots["Y-type"].to_html(
            full_html=False,div_id="y-type-plot",include_plotlyjs=False) if "Y-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,

        z_oms = oms_plots["Z-type"].to_html(
            full_html=False,div_id="z-type-plot",include_plotlyjs=False) if "Z-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,

        unknown_oms = oms_plots["Unknown-type"].to_html(
            full_html=False,div_id="unknown-type-plot",include_plotlyjs=False)  if "Unknown-type" in oms_plots else "<p style='font-style: italic;'>No data for this kind of N-ter</p>"  ,

        genome_over_time = make_genome_over_time_chart(cnx).to_html().split("<body>")[1].split("</body>")[0],

        sequence_over_time = make_sequence_over_time_chart(cnx).to_html().split("<body>")[1].split("</body>")[0],

        metrics_fig = make_genome_pie_chart(cnx).to_html(config= dict(displayModeBar = False)).split("<body>")[1].split("</body>")[0]
    )
        
    save(outfile,report)


