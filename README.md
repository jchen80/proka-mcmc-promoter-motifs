# proka-mcmc-promoter-motifs

## Overview
This repository contains code to run Markov chain Monte Carlo randomization as described in Itzkowitz et al. 2010. The goal is to generate a null distribution of randomized genomes by swapping codons which share the same amino acid translation and flanking codon context. 

To run, download `*cds_from_genomic.fna`, `*translated_cds.faa` and `*genomic.fna` files from the Refseq assembly for each of the species of interest. `mcmc/cds_randomizer.py` either runs in `kmer` mode where it computes kmer frequencies of kmers of length 6 and 7 in the real vs randomized genomes, or alternatively can be run in `plot` mode where it plots the relative changes in kmer frequencies as a function of the number of iterations. In the `kmer` mode the default number of iterations is `3nlog(n)` where `n` is the total number of nucleotides over all the coding segments.

```
python3 mcmc/cds_randomizer.py cds_from_genomic.fna translated_cds.faa results/ --mode kmer
```

From there, one can then calculate the relative enrichment of a kmer in the real vs randomized genomes. For example, this repository includes some scripts to calculate enrichment scores and statistics on the output of the `mcmc/cds_randomizer.py`:

```
python3 kmer_analysis/calculate_kmer_stats.py <Path to kmer frequency csv file>
                                              <length of kmer>
                                              <domain (archaea or bacteria)>
                                              <filepath to REBASE.txt file containing restriction enzyme sites>
```

Additionally, one can generate scatter plots of enrichment of specific motifs of interest vs covariates. For example:
```
python3 kmer_analysis/extract_motifs <parent directory (sub-directories per species)>
                                     <csv file containing tRNA scan learned motifs>
                                     <csv file containing covariates of interest, e.g. NAP abundance>
                                     -d <domain (either bacteria or archaea)>
                                     --k <length of kmer>

python3 kmer_analysis/plot_scatter_covariate.py <csv file from extract_motifs.py summarizing enrichment scores per motif>
                                                <name of covariate column>
```
This repository also includes scripts to fit phylogenetically adjusted model to the covariate vs enrichment relationship using `phylolm`.

### References 
1. Itzkovitz, S., Hodis, E., & Segal, E. (2010). Overlapping codes within protein-coding sequences. Genome Research, 20(11), 1582–1589. https://doi.org/10.1101/gr.105072.110
2. Ho LST, Ane C (2014). “A linear-time algorithm for Gaussian and non-Gaussian trait evolution models.” Systematic Biology, 63, 397-408.

