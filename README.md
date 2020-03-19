# SweepVsTraitOptimum

This folder contains all the scripts needed to reproduce the results in Barghi N and Schlötterer C. (2020). Distinct patterns of selective sweep and polygenic adaptation in evolve and re-sequence studies. Genome Biology and Evolution.

## prep_files_mimicree2.py

Using this script you will be able to produce haplotype and selection files (sel.txt for sweep and effect_size.txt for quantitative trait simulations) to be used for MimicrEE2. Phased haplotypes are available from Dryad (doi:10.5061/dryad.744p394)

## scale_new_trait_optimum.py

The new trait optimum needs to be changed if the number of loci and their effect sizes change so that the shift in trait optimum is uniform among all simulations for the sake of comparison. Using this script you will be able to change the new trait optimum.

## Run simulations

MimicrEE2 was used to perform simulations (https://sourceforge.net/projects/mimicree2/). Haplotype file and selected loci can be generated using prep_files_mimicree2.py. For quantitative trait simulations, the fitness function can be generated using scale_new_trait_optimum.py. The D. simulans recombination map (Dsim_recombination_map_LOESS_100kb_1.txt) is available from Dryad (doi:10.5061/dryad.744p394). Because D. simulans males do not recombine, we divided the recombination rate estimates by two. For neutral simulations the procedure is the same as selection simulations except that the effect size/selection coefficient of alleles are set to zero. Sample command lines for running sweep (w) and quantitative trait (qff) simulations are provided below:

```
java -jar mim2-v193.jar w --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --fitness sel.txt --output-sync AF_results.sync --output-gpf fitness_result.txt

java -jar mim2-v193.jar qff --haplotypes-g0 FlLines_FreeBayes_biallelicSNPs_q50_4500Ne_No4X_100loci.mimhap --recombination-rate Dsim_recombination_map_LOESS_100kb_1.txt --chromosome-definition "2=2L+2R,3=3L+3R" --replicate-runs 500 --snapshots 10,20,30,40,50,60,70,80,90,100,110,120,130,140 --effect-size effect_size.txt --heritability 0.5 --fitness-function sel.txt --output-sync AF.sync.gz --output-gpf GenoPheno.gpf

```
