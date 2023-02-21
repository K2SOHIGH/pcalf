# pyCALF

            _ (.".) _    
           '-'/. .\'-'   
             ( o o )     
              `"-"`  


pycalf stand for PYthon CALcyanin Finder

Potential calcyanins are searched in a set of amino acid sequences using a weighted HMM profile specific of the glycine zipper triplication (aka GlyX3). Then sequences with at least one hit against this profile are annotated with three HMM profiles specific of each Glycines zipper (Gly1, Gly2 and Gly3). A set of known N-ter sequences is used to infere the type (X, Z, Y or CoBaHMA) of the N-ter extremity. Finally, a flag is assign to each calcyanin depending on its N-terminus type and its C-terminus modular organization.

## Dependencies :
```
- pandas
- pyhmmer
- biopython
- numpy==1.22.4
- tqdm
- pdoc3
- blast
```
 
## Input :
CDS in fasta format (faa) [gzip or not].
If you want to process multiple files at once you can pass them through the --input argument with a tabular file (filename  path/to/cds) as entry.

## Output :
```
.
├── HMM/*.hmm
├── MSA/*.msa.fa
├── pycalf.features.tsv
├── pycalf.summary.tsv
└── N-ter-DB.tsv (update)
```

## INSTALLATION

```bash
mamba create -n pycalf -c bioconda blast pyhmmer pandas numpy pdoc3 biopython tqdm python=3.9 && conda activate pycalf;
pip3 install git@github.com:K2SOHIGH/pycalf2.0.git
```

## USAGE
```bash
pycalf -f fasta.cds.fa.gz -o pycalf_res
```