## Sequence capture probes and targets

This directory contains FASTA files for the probes we designed for sequence
capture, as well as the full-length gene sequences that we used as a target
file for assembling the sequence capture data. Detail about the design and
implementation of the array is available in a
[manuscript pre-print here](https://doi.org/10.1101/2020.05.08.084640).

### Contents

1. [README.md](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/probes/README.md): this file

2. [Salicaceae_phylogeny_probes.fasta](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/probes/Salicaceae_phylogeny_probes.fasta): Sequences of probes designed from the *Salix purpurea* 94006 genome annotation v1. The metadata line includes the *S. purpurea* gene name with the start and stop positions of the probe, as well as the ortholog in the *P. trichocarpa* genome assembly.

3. [Salicaceae_phylogeny_probes_supplement.fasta](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/probes/Salicaceae_phylogeny_probes_supplement.fasta): Probes designed from the
*Idesia polycarpa* genome assembly, to assist in target capture of regions
of the genome that are more highly divergent within the family. The metadata
line includes the *S. purpurea* gene name with the start and stop positions
of the probe, as well as the region of the *I. polycarpa* genome assembly that
was orthologous.

4. [phyloTargets.fasta](https://github.com/BrianSanderson/phylo-seq-cap/blob/master/probes/phyloTargets.fasta): The target file used in HybPiper to assemble
the sequence capture data. This includes all genes targeted by the probes, as
well as untargeted paralogs in the *S. purpurea* genome assembly. This file
includes the full length gene coding sequence, rather than just the targeted 
region, to assist with the recovery of off-target sequences and introns via
HybPiper. **Note**: each sequence has the string 'Sapur-' appended to the
the beginning, which is a required convention for use with HybPiper.
