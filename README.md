# SweepVsTraitOptimum

This folder contains all the scripts needed to reproduce the results in Barghi N and Schlötterer C. (2020). Distinct patterns of selective sweep and polygenic adaptation in evolve and re-sequence studies. Genome Biology and Evolution.

## Prepare files for simulations

Using the below script you will be able to produce haplotype and selection files (sel.txt for sweep and effect_size.txt for quantitative trait simulations) to be used for MimicrEE2. Phased haplotypes are available from Dryad (doi:10.5061/dryad.744p394)

```
prep_files_mimicree2.py
```

The new trait optimum needs to be changed if the number of loci and their effect sizes change so that the shift in trait optimum is uniform among all simulations for the sake of comparison. Using this script you will be able to change the new trait optimum.

```
scale_new_trait_optimum.py
```

## Run simulations

MimicrEE2 was used to perform simulations (https://sourceforge.net/projects/mimicree2/). Haplotype file and the file for selected loci can be generated using prep_files_mimicree2.py. For quantitative trait simulations, the fitness function can be generated using scale_new_trait_optimum.py. The D. simulans recombination map (Dsim_recombination_map_LOESS_100kb_1.txt) is available from Dryad (doi:10.5061/dryad.744p394). Because D. simulans males do not recombine, we divided the recombination rate estimates by two. For neutral simulations the procedure is the same as selection simulations except that the effect size/selection coefficient of alleles are set to zero. Sample command lines for running sweep (w) and quantitative trait (qff) simulations are provided below:

```
java -jar mim2-v193.jar w --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --fitness sel.txt --output-sync AF_results.sync --output-gpf fitness_result.txt

java -jar mim2-v193.jar qff --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_4500Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --effect-size effect_size.txt --heritability 0.5 --fitness-function sel.txt --output-sync AF.sync.gz --output-gpf GenoPheno.gpf

```

To store haplotypes for specific time points

```
java -jar mim2-v193.jar w --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap --recombination-rate dsim.rr_LDjump-LOESS.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 50 --snapshots 20,40,80,140 --fitness sel.txt --output-sync AF_results.sync --output-gpf fitness_result.txt --output-dir evolved_haplos

```
## Analyze and plot the results

You can analyze and plot the results of simualtions similar to our manuscript using the below scripts:

Plot the allele freqeuncy trajectories of selected alleles in sweep and trait optimum models (Figures 1, 6, 8, S2, S6, S7, S10 and S12) using:

```
AF_traj_plots.py
```

Plot the fitness (for sweep model) and phenotype (for trait optimum model) of populations (Figures 2, 5, 7, S5, S11) using:

```
fitness_phenotype_plots.py

```
The simulations for storing the haplotype information were run in 50 replicates only and haplotypes at generations 0 20,40,80, 140 were saved. Plot the distribution of selected alleles in haplotypes for sweep and trait optimum models (Figures 4 and 9) using:

```
haplotype_distribution.py
```

Compute and plot the similarity among 10 replicates in terms of the number of selected alleles (Figures 3 and S9) using: 

```
Jaccard_index_plots.py
```

Plot figures S1 (expected allele frequency change under neutrality in populations of 450 and 9000 individuals), S3 (allele frequency changes in a population of 450 until generation 2500 under trait optimum model), S4 (allele freqeuncies of 100 selected alleles in populations of 450 and 9000 individuals under sweep and trait optimum models), and S8 (variance of the fitness and variance of the number of selected alleles per haplotype in a population of 450 individuals with different number of selected alleles under sweep model) using:

```
misc_plots.py
```






