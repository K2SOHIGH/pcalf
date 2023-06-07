# PCALF

pcalf stand for Python CALcyanin Finder. 

Benzerara, K., Duprat, E., Bitard-Feildel, T., Caumes, G., Cassier-Chauvat, C., Chauvat, F., ... & Callebaut, I. (2022). A new gene family diagnostic for intracellular biomineralization of amorphous Ca carbonates by cyanobacteria. Genome Biology and Evolution, 14(3), evac026. [DOI:https://doi.org/10.1093/gbe/evac026](https://doi.org/10.1093/gbe/evac026)

Pcalf is a tool to retrieve Calcyanin protein and ccyA gene from genomes. 

## Calcyanin detection

The ccyA gene is searched at the protein level following a simple decision tree based on the specific composition of the C-ter extremity of this protein.
We use a weighted HMM profile specific of the glycine zipper triplication (aka GlyX3) to detect sequences with a putative glycine triplication. Sequences with at least one hit against this profile are annotated with three HMM profiles specific of each Glycines zipper (Gly1, Gly2 and Gly3). 
A set of known calcyanins is used to infere the type of the N-ter extremity (X, Z, Y, CoBaHMA or ?). Finally, a flag is assign to each sequence depending on 
its N-terminus type and its C-terminus modular organization.

## Decision tree

*picture here*

## Dependencies :
```       
         pyhmmer==0.7.1
         biopython==1.81
         numpy>=1.21
         pyyaml>=6
         snakemake==7.22
         pandas==1.5.3
         tqdm==4.64.1
         plotly==5.11.0
         python-igraph==0.10.4
```

## External dependency :
```
  blast
```

## INSTALLATION

```bash
mamba create -n pcalf -c bioconda blast pyhmmer pandas numpy biopython tqdm python=3.9 && conda activate pcalf;
pip3 install https://github.com/K2SOHIGH/pcalf.git
```

## Usage :

Pcalf is composed of several command :
- pcalf
- pcalf-datasets-workflow
- pcalf-annotate-workflow
- pcalf-report

### pcalf : 

This command can be used to look quickly for the presence of calcyanin in a set of amino acid sequences. It take one or more fasta files as input and output several files including a summary, a list of features and a list of raw hits produced during the search.  In addition, it also output updated HMMs, updateds MSAs with calcyanin tagged as Calcyanin with known N-ter detected if any.


### pcalf-datasets-workflow : 

pcalf-datasets-workflow can be used to retrieve genomes from NCBI databases such as RefSeq and GenBank based on accession (GC*_******.*) or TaxID
Genomes and annotations (CDS and genes ) will be downloaded using the new command line tools from [NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/). If annotations does not exists for a genome, then genes and CDS will be predicted with [Prodigal](https://github.com/hyattpd/Prodigal).

A yaml file is also produced and can be used as input for [pcalf-annotate-workflow](#pcalf-annotate-workflow-).
```
GC*_******.*:
  genome : /path/to/genome
  cds_faa: /path/to/cds_faa
  cds_fna: /path/to/cds_fna
GC#_######.#:
...
```

### pcalf-annotate-workflow : 

pcalf-annotate-workflow is actually the workflow that will help you recover calcyanin proteins and ccyA genes from a set of genomes.
This workflow is composed of multiple steps :
- genome taxonomic classification using [GTDB-TK V2](https://github.com/Ecogenomics/GTDBTk)
- genome quality assessment using [CheckM](https://github.com/Ecogenomics/CheckM/)
- calcyanin detection using [pcalf](#paclf-)
- calcyanin / ccyA linking using a specific python script.
- NCBI metadatas recovery using [NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/)

Note, that GTDB-TK and checkM requires external databases, respectively [GTDB](https://gtdb.ecogenomic.org/downloads) and [CheckM datas](https://data.ace.uq.edu.au/public/CheckM_databases). In addition, it's advised to run GTDB-TK and CheckM on a cluster like slurm. Because pcalf-annotate-workflow rely on snakemake you can easily provide a snakemake profile through the --snakargs option to run it on your favorite cluster. On the other hand, you can skip the genome taxonomic classification and the quality assessment with the --quick flag.

pcalf-annotate-workflow take as input a yaml file with a specific format, see [pcalf-datasets-workflow](#pcalf-datasets-workflow-) for details.

The workflow produced several files for each step but the final output is a sqlite3 database storing multiple table: 
- genome          # NCBI metadatas
- gtdbtk          # GTDB-TK classification results
- checkm          # Checkm Results 
- summary         # PCALF summary table
- features        # PCALF features table
- hits            # PCALF hits table
- ccyA            # ccyA table
- gly1            # Gly1 MSA
- gly2            # Gly2 MSA
- gly3            # Gly3 MSA
- glyx3           # Glyx3 MSA
- nterdb          # N-ter table

You can use the sqlite3 database from another run as a basis for a new one. In this case, MSAs stored in the sqlite3 file will be used to generate HMM profiles for [pcalf](#pcalf-).

### pcalf-report : 

This command produce an HTML report from a sqlite3 database produced by [pcalf annotate workflow](#pcalf-annotate-workflow-).


### Workflow : 

*picture here*





