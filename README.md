# Phillips_et_al_mtDNA
Code associated with Phillips et al. 2022, https://doi.org/10.1002/ece3.9426


1. phylogeny.R: script to generate all phylogeny figures in text and supplement.

2. all_haps_places_stages_25May22.csv: dataframe of global long-fragment mtDNA haplotypes from all seven extant sea turtle species and associated ocean basins and life stages where each haplotype has been found; for input to phylogeny.R

3. seaturtle_mtDNA_longfragment_raw_alm_clustal.fasta: unedited alignment of haplotypes generated with Clustal Omega

4. seaturtle_mtDNA_longfragment_alm.fasta: edited Clustal Omega alignment of haplotypes

5. seaturtle_mtDNA_longfragment_alm.nex: input nexus file for MrBayes

6. infile.nex.con.tre: MrBayes output consensus tree and input for phylogeny.R

7. run2_mtDNA_Strict_Yule_uniform.xml: beauti xlm file with priors specified; use as input to BEAST2

8. run2_seaturtle_mtDNA_longfragment_alm.trees: treeannotator output, input for phylogeny.R **Available upon request; file too large for GitHub upload*

9. run2_treeannotator_Yule_uniform.txt: output from BEAST2, use as input for TreeAnnotator.

10. seaturtle_mtDNA_longfragment_alm.fasta.treefilee: IQtree output, input for phylogeny.R
