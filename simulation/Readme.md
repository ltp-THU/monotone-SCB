[Toc]

# Overview

In our simulation scenarios, we utilize a high-performance computing cluster to conduct simulations using bootstrap algorithms. These simulations vary in sample size, dimension, dependence structure, and bandwidth settings. This Readme file provides a detailed workflow for reproducing the simulation results presented in the tables and figures of our submission, organized into the following sections.


# Reproduction of Table 2
The Table 2 in the main paper shows our main simulation results. In Table 2, we run 400 simulations for 2 models (model (a) and model (b)) under 3 sample sizes ($n=300,500,1000$) and 3 dimensions ($p=9,18,27$). We allocate the simulation into $2\times 3\times 3=18$ tasks according to the model, sample size and dimension settings. 

All the simulation results are arranged as folders under directory `simulation/table_2&B.1&Figure_B.2/`. For example, the directory `simulation/table_2&B.1&Figure_B.2/a300p9sigma1` contains the simulation results of model (a) with sample size $n=300$ and $p=9$ and bandwidth parameters $h_r=0.05,0.1,0.15,0.2,0.3,0.35$. 

The folder `simulation/table_2&B.1&Figure_B.2/GCV` also contains these $2\times 3\times 3=18$ tasks with bandwidth parameter chosen by GCV strategy. 

The folder `simulation/table_2&B.1&Figure_B.2/summary` contains the R script to summarize all the outputs of these tasks and reproduce the Table 2 in the submitted main paper.


Under directory `simulation/table_2&B.1&Figure_B.2/`, we use single R script `simulation_table2.R` to run all the simulation scenarios. The `.sh` file under each subdirectory helps us use different inputs (model/sample size/dimension/bandwidth...) to run the R script on a high-performance computing cluster. 

In specific, we split the 400 simulations into 10 batches and each batch runs 40 simulations on 1 node with 10 cores. Take the folder `simulation/table_2&B.1&Figure_B.2/a300p9sigma1` for example, the `a300p9.sh` file contains code of the form:
```
#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --array=1-10
#SBATCH --job-name="a300p9"
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --output=%x.%a.out
......
......
# Define parameters within the script
index=$SLURM_ARRAY_TASK_ID          # Batch index (1-10)
# Define other parameters (these values are fixed for all jobs)
model="a"                           # Model type ("a" or "b")
sample_size=300                      # Sample size
dimension=9                        # Dimension
candidate="10:100"                   # Candidate vector (in R format)
bandwidth="0.05 0.1 0.15 0.2 0.25 0.3 0.35"  # Bandwidth vector
GCV="FALSE"                          # GCV parameter (TRUE or FALSE)
simulation=40                       # Simulation for each batch (e.g., 40 for batch 1-10 means 400 simulations in total)

Rscript ../simulation_table2.R $index $model $sample_size $dimension "$candidate" "$bandwidth" "$GCV" $simulation
```
We define all the settings including model (model = "a"), sample size (smaple_size=300), dimension (dimension=9), and tunning parameters in `.sh` file. If users want to lauch task for model b with sample size $n=300$, dimension $p=9$, they only need change the inputtings in `.sh` file to run the same underlying code `simulation_table2.R`. The `#SBATCH --array=1-10`, `#SBATCH --cpus-per-task=10`, and `ARGS=$(printf "%02d" $SLURM_ARRAY_TASK_ID)` means that we launch 10 tasks. Each task uses 1 node with 10 cores and will output a slurm-out file names as `a300p9.X.out`(X=1-10). The `.out` file will report the time usage and SCB's coverage/width results for different bandwidths simultaneously with the execution. The final output of R script would be 10 `.csv` files, each referring to a batch of 40 simulations, resulting in a total of 400 simulations. Finally, the `summary.R` in directory `simulation/table_2&B.1&Figure_B.2/summary/` would summarize all these `.csv` files and reproduce Table 2 as `table2.csv` under directory `simulation/table_2&B.1&Figure_B.2/summary/`.

**Note:** The structures of the subdirectory of `a1000p27sigma1/` and `b1000p27sigma1/` is differen with others. This is not a mistake; rather, it is due to the allocation of computing tasks for faster and more efficient computation. The `a1000p27' scenario represents the largest scenario in our simulation with sample data $n=1000$ and dimension $p=27$, which requires the most time and memory to compute. In this particular case, we split the bandwidth range $h_r = 0.05 - 0.35$ into smaller subranges ($0.05 - 0.1$, $0.15 - 0.2$, and $0.25 - 0.35$) and run three separate computing tasks for each subrange. Such splitting facilitates faster computation by utilizing more nodes. The results from these tasks are then summarized together. In other scenarios, where the sample size and dimension are smaller, we run all the bandwidth choices together.



## Reproduction of Table B.2&Figure B.2
The Table B.2 and Figure B.2 are presented in our supplement to show additional information for the simulation results in Table 2 of the main paper. Based on those `.csv` simulation outputs, the `summary.R` in directory `simulation/table_2&B.2&Figure_B.2/summary/` would reproduce Table B.2 and Figure B.2 as `tableB.2.csv` and `Figure B.2.png`. 


# Reproduction of Table B.3
The codes and results for reproducing Table B.3 are in the directory `simulation/table_B.3/`. The workflow is similar with the reproduction of Table 2. There are $2\times 3\times 3=18$ tasks referring to 2 models (model (a) and model (b)), 3 sample sizes ($n=300,500,1000$) and 3 dimensions ($p=9,18,27$). Each task is put in a task folder under directory `simulation/table_B.3/` e.g. folder `a300p9zhou` refers to model (a), sample $n=300$ and dimension $p=9$. We split the 400 simulations into 10 batches and each batch runs 40 simulations on 1 node with 10 cores. Use the `.sh` files in the 18 task folders to submit R script to the SLURM job scheduler and obtain 10 `.csv` output files for each task. Based on the outputs, the `summary.R` in directory `simulation/table_B.3/` would reproduce Table B.3 as `tableB.3.csv`.

# Reproduction of Table B.4

Table B.4 in our supplement presents the simulation with skewness error innovation. The codes and simulation results of two types of skewed distribution, exponential and log-normal are in the directory `simulation/table_B.4/exp/` and `simulation/table_B.4/lognormal/` respectively. 

We also split the 400 simulations into 10 batches and run them on a high-performance computing cluster. In the directory `simulation/table_B.4/exp/`, one can use the `n500LS_exp.sh` file to submit R script `n500p27exp.R` to run 10 batches and output 10 `.csv` files. In the directory `simulation/table_B.4/lognorm/`, one can use the `n500LS_lognorm.sh` file to submit R script `n500p27lognorm.R` to run these 10 batches and output 10 `.csv` files. The R script `simulation/table_B.4/summary.R` would summarize these `.csv` files and reproduce the Table B.4 as `tableB.4.csv`.

## Reproduction of Figure B.1

The Figure B.1 in our supplement presents the Q-Q plot for our Gaussian approximation under the skewed distribution setting in Table B.3. The codes and results are in directory `simulation/figure_B.1/exp` and `simulation/figure_B.1/lognorm`. One can use the `sub_GA.sh` to run R scripts `GA_exp.R` and `GA_lognorm.R`. The code will output `.csv` file named as `GA_verify_exp.csv` or `GA_verify_lognorm.csv`. Based on these two `.csv` results, one can run the R script `GA_qqplot.R` to reproduce Figure B.1 as `Figure B.1.png`.

# Reproduction of Table B.5
Table B.5 in our supplement presents the simulation for our penalized SCBs when trend functions have flat period. The workflow is similar with the reproduction of Table 2. There are $2\times 4=8$ tasks referring to 2 models (model (a) and model (b)), 4 choices of penalization coefficient $C_1=0,0.1,0.2,0.3$. Each task is put in a task folder under directory `simulation/table_B.5/` e.g. folder `/a300p27penal0.1` refers to model (a), sample $n=500$ and dimension $p=27$ with penalization coefficient $C_1=0.3$. We split the 400 simulations into 10 batches and each batch runs 40 simulations on 1 node with 10 cores. Use the `.sh` files in the 8 task folders to submit the R script `penalization.R` to the SLURM job scheduler and obtain 10 `.csv` output files for each task. Based on the outputs, the `summary.R` in directory `simulation/table_B.5/` would reproduce Table B.5 as `tableB.5.csv`.


# Reproduction of Table B.6

The Table B.6 in the supplement shows our simulation results for regression models. In Table B.6, we run 400 simulations for 2 models (model (c) and model (d)) under 3 sample sizes ($n=300,500,1000$). We allocate the simulation into $2\times 3=6$ tasks according to the model and sample size settings. All the codes and simulation results are arranged as folders under directory `simulation/table_B.6/`. For example, the directory `simulation/table_B.6/c300` contains the simulation task and results of model (c) with sample size $n=300$ and bandwidth parameters $h_r=0.05,0.1,0.15,0.2,0.3,0.35$. 

The folder `simulation/table_B.6/GCV` also contains these $2\times 3=6$ tasks with bandwidth parameter chosen by GCV strategy. For example, the code and results for model (c) with sample size $n=300$ with bandwidth chosen by GCV is in the directory `simulation/table_B.6/GCV/c300gcv`.

We run the simulation task on a high-performance computing cluster. In specific, we split the 400 simulations into 10 batches and each batch runs 40 simulations on 1 node with 10 cores. One can use the `.sh` file in each task folder to submit R script to the SLURM job scheduler and obtain 10 `.csv` output files for each task. Based on the outputs, the `summary.R` in directory `simulation/table_B.6/` would reproduce Table B.6 as `tableB.6.csv`.

