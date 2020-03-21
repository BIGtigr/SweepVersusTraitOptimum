# SweepVsTraitOptimum

This folder contains all the scripts needed to reproduce the results in Barghi N and Schlötterer C. (2020). Distinct patterns of selective sweep and polygenic adaptation in evolve and re-sequence studies. Genome Biology and Evolution.

## 1. Prepare files for simulations

Using the below script you will be able to produce haplotype and selection files (sel.txt for sweep and effect_size.txt for quantitative trait simulations) to be used for MimicrEE2. To use this script you need 189 Drosophila simulans phased haplotypes that are available from Dryad (doi:10.5061/dryad.744p394).

```
prep_files_mimicree2.py
```

The new trait optimum needs to be changed if the number of loci and their effect sizes change so that the shift in trait optimum is uniform among all simulations for the sake of comparison. Using this script you will be able to change the new trait optimum.

```
scale_new_trait_optimum.py
```

## 2. Run simulations

MimicrEE2 was used to perform simulations (https://sourceforge.net/projects/mimicree2/). Haplotype file and the file for selected loci can be generated using prep_files_mimicree2.py. For quantitative trait simulations, the fitness function can be generated using scale_new_trait_optimum.py. The D. simulans recombination map (Dsim_recombination_map_LOESS_100kb_1.txt) is available from Dryad (doi:10.5061/dryad.744p394). Because D. simulans males do not recombine, we divided the recombination rate estimates by two. For neutral simulations the procedure is the same as selection simulations except that the effect size/selection coefficient of alleles are set to zero. Sample command lines for running sweep (w) and quantitative trait (qff) simulations are provided below:

```
java -jar mim2-v193.jar w --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --fitness sel.txt --output-sync AF_results.sync --output-gpf fitness_result.txt

java -jar mim2-v193.jar qff --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_4500Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --effect-size effect_size.txt --heritability 0.5 --fitness-function sel.txt --output-sync AF.sync.gz --output-gpf GenoPheno.gpf

```

To study the distribution of selected alleles on haplotypes, in addition to the allele frequencies of populations you need to store the haplotypes too. We stored haplotypes for specific time points using fewer iterations. A sample command line is provided below:

```
java -jar mim2-v193.jar w --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap --recombination-rate dsim.rr_LDjump-LOESS.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 50 --snapshots 20,40,80,140 --fitness sel.txt --output-sync AF_results.sync --output-gpf fitness_result.txt --output-dir evolved_haplos

```
## 3. Analyze and plot the results

You can analyze and plot the results of simualtions similar to our manuscript using the below scripts:

3.1 You can plot the allele freqeuncy trajectories of selected alleles in sweep and trait optimum models (Figures 1, 6, 8, S2, S6, S7, S10 and S12) using the below script. To run this script you need the output of simulations in sync format for simulations with selection and  neutral simulations. you also need the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations).

```
AF_traj_plots.py
```

3.2 You can plot the fitness (for sweep model) and phenotype (for trait optimum model) of populations (Figures 2, 5, 7, S5, S11) using the below script. You only need the genotype/phenotype/fitness result of simulations.

```
fitness_phenotype_plots.py

```
3.3 In our manuscript, the simulations for storing the haplotype information were run in 50 iterations only and haplotypes at generations 0 20, 40, 80, and 140 were saved. You can plot the distribution of selected alleles in haplotypes for sweep and trait optimum models (Figures 4 and 9) using the below script. To run this script you need the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations) and the path to the folder where haplotypes are stored.

```
haplotype_distribution.py
```

3.4 You can compute and plot the similarity among replicates in terms of the number of selected alleles (Figures 3 and S9) using the below script. We have divided 5000 iterations into 50 experimental evolution with 10 replicates anda similarity among 10 replicates are computed. To run this script you need the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations) and the output of simulations in sync format for simulations with selection and neutral simulations.

```
Jaccard_index_plots.py
```

3.5 Using the below script you can plot:
figures S1 (expected allele frequency change under neutrality in populations of 450 and 9000 individuals)
figure S3 (allele frequency changes in a population of 450 until generation 2500 under trait optimum model)
figure S4 (allele freqeuncies of 100 selected alleles in populations of 450 and 9000 individuals under sweep and trait optimum models) figure S8 (variance of the fitness and variance of the number of selected alleles per haplotype in a population of 450 individuals with different number of selected alleles under sweep model). 

For figure S1 you need the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations) and the output of simulations in sync format for neutral simulations. 

For figure S3 you need to run the simulations in the small population for 2500 generation under quantitative trait model (qff function in MimicrEE2) until all alleles are lost/fixed and you need the list of selected sites (effect_size.txt for quantitative trait simulations) and the output of simulations in sync format. 

For figure S4 you need the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations) and the output of simulations in sync format for small and big populations under sweep and trait optimum models.

For figure S8 you need the genotype/phenotype/fitness result of sweep and trait optimum simulations for both small and large populations,  the list of selected sites (sel.txt for sweep and effect_size.txt for quantitative trait simulations) and the path to the folder where haplotypes are stored.

```
misc_plots.py
```

The script were written in Python 2.7 (I know, I know I have to switch to Python 3). 

For any questions please contact me at barghi.neda@gmail.com. 





