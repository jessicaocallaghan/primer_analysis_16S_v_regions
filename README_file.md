## Readme file for reanalysis of public data to determine primer bias in upper gential tract microbiome studies. 

Each individual analysis will describe the methods use and then go into detail on the scripts within each directory.



1. *Taxonomic classification of next generation sequencing of original dataset*
  Sequence clustering and operational taxonomic unit (OTU) selection was performed using a modified version of CD-HIT-OTU-454 which does not remove singleton clusters [1]. Taxonomy was assigned to representative sequences by comparison to the latest build of the [Greengenes database][] using [BLAST][], and OTU tables were constructed from the output using a custom Perl script [2]. 

  [Greengenes database]: https://greengenes.secondgenome.com/
  [BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi

  

2. *Phylogenetic analysis*
  Was completed in GUI software (Geneious and MEGA) not by coding. 
  Full-length 16S rRNA sequences for *Lactobacillus* spp. (accession numbers: __AB680529.1, AB690249.1, AB668940.1.1, AB008203.1.1, AB425941.1.1, AB008206.1, AF243169.1, AF243167.1, CP018809.253324, CP018809.1516019, CP018809.1347636, CP018809.500868, AB547127.1, AB517146.1, AB932527.1, AB008209.1, HZ485829.7, LG085736.7 LF134126.7, LG104504.7__), *Pediococcus pentosaceus* (accession numbers: __AB018215.1 and AB362987.1__), and *Bacillus subtilis* (accession numbers: __AP012496.9810 and AP012496.30276__) were downloaded from the [SILVA database][] using the web interface. Sequences were aligned using ClustalW [3] with the default settings. MEGA7 [4] was used to generate the best-known maximum likelihood (ML) tree using a Jukes-Cantor model and 1000 bootstrapping iterations.

  [SILVA database]: www.arb-silva.de

  To generate the phylogenetic trees for the specific regions,  V5-V8 region phylogenetic tree the same full-length 16S rRNA sequences from above were imported into [Geneious][] along with two *Escherichia coli* sequences downloaded from the [SILVA database][] (accession number: __AB045730.1 and AB045731.1__). Sequences were aligned using standard Geneious alignment and trimmed to include the variable specific variable regions (V1-V2, V1-V3, V3-V4, V3-V4. V4 and V5-V8) using the _E.coli_ sequences. Trimmed sequences were then aligned using ClustalW with default settings and imported into MEGA7 and a trees was were constructed and edited same as above. 

  [Geneious]:  https://www.geneious.com/

3. *Hierarchial clustering*
  A dissimilarity matrix was generated based on the relative abundances of *Lactobacillus* spp. in the pyrosequenced and qPCR analysed samples using the vegdist function in the vegan package in [R][]  with the Bray-Curtis dissimilarity metric [5]. Hierarchical clustering was performed using the hclust function in R with ‘average’ linkage (UPGMA) [6]. Clustering and relative abundances were visualized in a heatmap with associated dendrogram using the heatmap.2 function from the R package ggplots [7]. 

  [R]: https://www.r-project.org/

4. *Reanalysis of publicly available endometrium datasets*
  Four datasets (Table 1) were re-analysed using the standard [QIIME2][] data analysis pipeline to determine Lactobacillus abundance using a standard full-length taxonomic classifier. 

  [QIIME2]: http://qiime.org/

  __Table 1__: Datasets used for reanalysis

  | Dataset          | Samples                                                      |
  | ---------------- | ------------------------------------------------------------ |
  | PRJNA329174 [8]  | Endometrial fluid from fertile women at distinct time points before and during pregnancy |
  | PRJNA546137 [9]  | Vaginal fluid, healthy endometrium and endometrium lesions (endometriosis) |
  | PRJEB18626 [10]  | Vaginal, cervical and endometrial samples of control and women with a history of infertility |
  | PRJNA543861 [11] | Mid endometrial samples from women undergoing a hysterectomy and compared with other parts of the reproductive tract |


  The raw reads were downloaded from NCBI using fastq-dump and, if paired end, reads were split. Raw data was then uploaded into QIIME2 using standard parameters and quality control checks were run. Based on the QC data, denoising and trimming was completed with custom values for each dataset using the DADA2 algorithm within QIIME2. The denoised DNA sequences and corresponding biom table were downloaded from QIIME2 and used in R with a full length taxonomic classifier and the DADA2 software. Custom R scripts were used to identify the number of sequences that had been classified as a Lactobacillus species from the full dataset.   

  Each primer pair was used to create in silico 16S rRNA gene read libraries using custom databases with the five *Lactobacillus* species of interest (downloaded from SILVA) using [Grinder][], a sequencing simulator. 

  [Grinder]:  https://sourceforge.net/projects/grinder/

  Custom scripts in R were used to determine the sequence coverage (the percentages of reads amplified by the primers) between different primer pairs. Each primer pair was used with the custom database for all *Lactobacillus* species (n=1371) to create fastq files with approximately 20,000 sequencing reads. The proportion of reads belonging to each species was identified for each primer and compared to the predicted amplification from the unbiased total dataset. These values correspond to the first column (number of reads) and second column (percentage of reads for each species) within the results table for each primer pair. The proportion of species that was expected to be amplified was compared to the actual amplification predicted by each primer set using custom scripts in R. 

  Primer pairs were also used to *in silico* amplify the human genome to predict low biomass contamination by host DNA using [MFEprimer][] software.

  [MFEprimer]:  https://www.mfeprimer.com/

_References_

1. Liu B, Gibbons T, Ghodsi M, Treangen T, Pop M. Accurate and fast estimation of taxonomic profiles from metagenomic shotgun sequences. BMC Genomics. 2011;12 Suppl 2:S4
2. McDonald D, Price MN, Goodrich J, Nawrocki EP, DeSantis TZ, Probst A, et al. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea. ISME J. 2012 Mar;6(3):610–8. 
3. Thompson JD, Higgins DG, Gibson TJ. CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. Nucleic Acids Res. 1994 Nov 11;22(22):4673–80.
4. Kumar S, Stecher G, Tamura K. MEGA7: Molecular Evolutionary Genetics Analysis Version 7.0 for Bigger Datasets. Mol Biol Evol. 2016;33(7):1870–4. 
5. Oksanen J, Kindt R, Legendre P, Hara B, Simpson G, Solymos P, et al. The vegan Package. 2009 Jan 19; 
6. R Core Team. R: A language and environment for statistical computing. [Internet]. R Foundation for Statistical Computing; 2014. Available from: http://www.R-project.org/
7. Wickham H. ggplot2: Elegant Graphics for Data Analysis. [Internet]. 2016. Available from: https://ggplot2.tidyverse.org/
8. Moreno I, Codoñer FM, Vilella F, Valbuena D, Martinez-Blanch JF, Jimenez-Almazán J, et al. Evidence that the endometrial microbiota has an effect on implantation success or failure. Am J Obstet Gynecol. 2016 Dec;215(6):684–703.
9. Hernandes C, Silveira P, Rodrigues Sereia AF, Christoff AP, Mendes H, Valter de Oliveira LF, et al. Microbiome Profile of Deep Endometriosis Patients: Comparison of Vaginal Fluid, Endometrium and Lesion. Diagnostics (Basel) [Internet]. 2020 Mar 17 [cited 2020 Dec 15];10(3). Available from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7151170/
10. Wee BA, Thomas M, Sweeney EL, Frentiu FD, Samios M, Ravel J, et al. A retrospective pilot study to determine whether the reproductive tract microbiota differs between women with a history of infertility and fertile women. Australian and New Zealand Journal of Obstetrics and Gynaecology. 2018;58(3):341–8.
11. Winters AD, Romero R, Gervasi MT, Gomez-Lopez N, Tran MR, Garcia-Flores V, et al. Does the endometrial cavity have a molecular microbial signature? Sci Rep [Internet]. 2019 Jul 9 [cited 2020 Sep 2];9. Available from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6616349/