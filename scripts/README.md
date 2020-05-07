## Design and implementation of sequence capture array for phylogenetics

This directory contains Python scripts that I wrote to design the sequence
capture array reported in this manuscript, and to quantify nucleotide diveristy
in the resulting sequence data

### Contents

1. [README.md](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/scripts/README.md): this file

2. [poly_screen.py](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/scripts/poly_screen.py): Screen SNP and indel polymorphism in a VCF, used here to design a targeted sequence capture array

3. [site_degeneracy.py](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/scripts/site_degeneracy.py): Characterize the degeneracy for each site in a FASTA file using the standard
codon table

4. [calculate_pi.py](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/scripts/calculate_pi.py): Summarize variation in nucleotide sequences between individuals within a
species as Nei's &pi;.
