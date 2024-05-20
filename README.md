## IntroAdapt

IntroAdapt makes (*will eventually make*) a joint inference of introgression and selection from longitudinal population genomic data. The inference is based on approximate Bayesian computation via random forest. It takes whole genome data (using SLiM + pyslim/msprime).

### Requirements

IntroAdapt is a collection of scripts in Python and SLiM. They have been tested in an Ubuntu (20.04) machine using a Conda environment and using a Snakemake workflow to run them. The Conda environment was created with the following commands (see file `introadapt.yml` to get version number of each package):

```shell
conda create -n introadapt python==3.8.10 r-base=3.6.3
conda activate introadapt
pip install tskit==0.3.7
pip install msprime==1.0.2
pip install pyslim==0.600
pip install scikit-allel==1.3.5
pip install pandas
pip install scipy
pip install pytest
pip install flake8
conda install slim=3.7
conda install -c conda-forge demesdraw
conda install -c r r-rcarbon=1.2.0
conda install -c r r-ini
conda install -c r r-extraDistr
conda install -c r r-testthat
conda install -c r r-dplyr
conda env export > introadapt.yml
```

Reminder for creating environment from file:
```
conda env create -f introadapt.yml
```

Install `abcrf` in R via `install.packages()` with appropriate path:

```r
install.packages('abcrf','/home/USER/anaconda3/envs/introadapt/lib/R/library').
```

### Usage

Input files:

- *config file*: a ini file with options regarding the Setting, Model, Prior and Statistics to be used in your analysis

- *sample file*: a text file (in the form of a space separated table) with metadata of your samples (ID, age, sequencing coverage, etc)

- *genome file*: a text file with a description of the genome organization (chromosomes, recombination map)

- *data file*: a vcf file with the genetic data

To run different parts of the analysis with snakemake :

```shell
snakemake RULE -C options_file='PATH/TO/YOUR/CONFIG_FILE.ini'
```
Where `RULE` is one of the rules defines in the snakefile. For instance, running `snakemake reftable -C options_file='tests/config_project.ini'` will create small reference table using parameters in file `tests/config_project.ini`. Typically the user will use rule `reftable` to run simulations and create the reference table (parameters, summary statistics and latent variables) and `abc` to grow random forests and perform approximate Bayesian computation inferences (rule `abc` will call `reftable` if there is no reference table; there is no need to do it in two steps but it can be usful in many cases).

Before running your analysis is highly recommended to performs some tests. Unit tests (using pytest) can be run using `snakemake test`.

![Directed acyclic graph for the workflow using the test project configuration (`snakemake --dag | dot -Tsvg > workflow_dag.svg`)](workflow_dag.svg)

#### Cleaning old files

You can remove all files (from a specific project or from all projects) from your results folder:

```shell
snakemake clean_project -C options_file='path/to/your/config_file.ini'
snakemake clean
```


### Input parameters (in *.ini config file)

| Parameter name | type | description |
|---|---|---------------|
|**[Settings]**|||
| *verbose* | integer | Level of information to output on screen (0 = minimal; 1,2,3,... = increasing detailed output; -1 = none).|
| *project_name* | string | name of the analysis project, a folder with that name is created in results folder and all output written inside.|
| *batch* | integer | identifier of a batch of simulations for the project. A folder named with that identifier is created within the project folder and results from the simulations of that bacth written inside.|
| *sample_file* | string |  path + file name containing sample information (see below)|
| *genome_file* | string | path + file name containing genome information (see below)|
| *num_of_sims* | positive integer | number of simulation (in the batch)|
| *seed* | integer | seed for the random number generator |
|**[Model]**|||
| *generations_coalescence* | positive integer | number of generations simulated in coalescence (by msprime : not used for the moment) |
| *generations_split* | positive integer |Generation when the population ancestral split occur in coalescence  (by msprime) |
| *generations_forward* | positive integer | number of generations simulated in forward (by SLiM)|
| *generations_migration* | positive integer | The generation when the migration start (by SLiM|
| *populations* | positive integer | Number of population in the simulation (compt the ancestral population too)|
| *selecvar* | positive integer or float | selection coefficient (s) (by SLiM). The user can indicate two numbers in this case IntroAdapt will generate a uniform randum S (runif()) between the two indicated number.|
| *rmig* | positive integer or float | migration rate (by SLiM)|
| *fmig* | positive integer or float | Number used to calculate the generation when the migration ends (by SLiM)|
| *Mtype* | String | The type of model used for the simulation. MAA: advantageous mutation in donor and recipient population or MNA: neutral mutation in donor population then advantageous in recipient (s = 0.0). |
| *Prandom* | String | The size of the populations of the model are random (sample()): "Yes". Or the size of the ancestral population (pop0) is the addition of the sizes of the daughter populations (pop1 and pop2) : "No", in this case the size of the ancestral population will be equal to pop_size_max in Priors.|
|**[Priors]**|||
| *gen_len_prior_sh1* | | prior for generation length (shape1 of rescaled beta distribution)|
| *gen_len_prior_sh2*| | prior for generation length (shape2 of rescaled beta distribution)|
| *gen_len_prior_min*| | prior for generation length (minimum of rescaled beta distribution)|
| *gen_len_prior_max*| | prior for generation length (maximum of rescaled beta distribution)|
| *pop_size_prior_min* | positive integer | prior for population size (minimum of uniform distribution)|
| *pop_size_prior_max* | positive integer | prior for population size (maximum of uniform distribution)|
| *mut_rate_prior_mean* |  |prior for mutation rate (mean of normal distribution)|
| *mut_rate_prior_sd* |  | prior for mutation rate (sd of normal distribution)|
|**[Vcf]**|||
| *chro* | String| List of chromosomes ("chro" + chromosome number) on which the variants are found. There must be the same number of chromosome names as variants. If the user fills in just "No" then this list will be created automatically from the information (chromosome + positions) saved in the genome.txt file using the chrom_deter() method of write_vcf.py.|
| *id* |  | List of identifiers for each variant. Optional if user specifies "No" instead.|
| *ref* |  | List of reference base(s) for each variant (A,C,G,T,N). Optional if user specifies "No" instead.|
| *alt* |  | List of alternate base(s) for each variant (A,C,G,T,N). Optional if user specifies "No" instead.|
| *qual* |  | List of quality score for the assertion made on the alternative bases. Optional if user specifies "No" instead.|
| *filter* |  | Does not work. Optional if user specifies "No" instead.|
| *info* |  |List of lists of additional information. Optional if user specifies "No" instead.|
| *format* |  | Does not work. Optional if user specifies "No" instead.|
| *isphased* |  | If the genotypes are phased or not. "Yes" or "No" parameter. "Yes" for phased genotypes and "No" for non-phased ones.|
|**[Tools]**|||
| *volcanofinder* |  | Create VolcanoFinder files (allele grequency and SFS) If the user want to this tools."Yes" or "No" parameter.|

Example:

```:tests/config_project.ini

```




### Input files

- *Sample information file*: 

- *Genome information file*:

...






