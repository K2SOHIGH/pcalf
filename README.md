[pcalf](https://k2sohigh.github.io/pcalf/) stand for Python CALcyanin Finder. 

pcalf pydoc is availbable [here](https://k2sohigh.github.io/pcalf/docs/src).


            _ (.".) _    
           '-'/. .\'-'   
             ( o o )     
              `"-"`  

Benzerara, K., Duprat, E., Bitard-Feildel, T., Caumes, G., Cassier-Chauvat, C., Chauvat, F., ... & Callebaut, I. (2022). A new gene family diagnostic for intracellular biomineralization of amorphous Ca carbonates by cyanobacteria. Genome Biology and Evolution, 14(3), evac026. [DOI:https://doi.org/10.1093/gbe/evac026](https://doi.org/10.1093/gbe/evac026)

Calcyanins are searched in a set of amino acid sequences using a weighted HMM profile specific of the glycine zipper triplication (aka GlyX3). 
Sequences with at least one hit against this profile are annotated with three HMM profiles specific of each Glycines zipper (Gly1, Gly2 and Gly3). 
A set of known calcyanins is used to infere the type of the N-ter extremity (X, Z, Y or CoBaHMA). Finally, a flag is assign to each calcyanin depending on 
its N-terminus type and its C-terminus modular organization.

## Dependencies :
```
pyhmmer==0.7.1
biopython==1.81
numpy==1.22.4
pyyaml==6.0
snakemake==7.22
pandas==1.5.3
tqdm==4.64.1
```

## External dependency :
```
blast
```


## Input :
CDS in fasta format (faa) [gzip or not] or a tabular file with  one "file designation\tpath/to/fasta/file" per line.

## Output :
```
.
├── HMM/*.hmm (updated)
├── MSA/*.msa.fa (updated)
├── pcalf.features.tsv
├── pcalf.summary.tsv
└── N-ter-DB.tsv (updated)
```

## INSTALLATION

```bash
mamba create -n pcalf -c bioconda blast pyhmmer pandas numpy biopython tqdm python=3.9 && conda activate pcalf;
pip3 install https://github.com/K2SOHIGH/pcalf.git
```

## USAGE
```bash
pcalf -i fasta.cds.fa.gz -o pcalf_res
```


