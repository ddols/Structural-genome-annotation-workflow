# Structural-genome-annotation-workflow
This is a step-by-step guide explaining some of the steps in /sdind's genome_annotation_workflow. All credit goes to they work and be sure to visit their repository for an in-depth dive into the workflow. This guide is meant for people who are new to structural genome annotation.

# Before we get started

This tutorial is meant to showcase one of the multiple ways a genome may be structurally annotated. That is, basically, to physically locate every possible genomic feature such exons, introns, genes and regulatory sequences (whenever it is possible). To do so it is fundamental to not only have a genome to annotate (duh) but most importantly transcriptomic data which will be essential to characterize the previously mentioned features.

In an ideal world, the transcriptome should belong to the same species the genome we are trying to annotate belongs to. However, in times of need a sufficiently close species' transcriptome should do the trick. And in more desperate times, well, any transcriptome more or less related to our taxon is acceptable. Anyway, this workflow is an adaptation of [this](https://github.com/sdind/genome_annotation_workflow?tab=readme-ov-file) workflow by **sdind**. 

Lastly, this workflow is written with `Linux` in mind as the running OS with access to a cluster with a [Slurm](https://slurm.schedmd.com/quickstart.html) system. If you are working on a [SGE](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) cluster or any other kind, you can still run the scripts in this repository, although you will have to modify them accordingly for your queueing system. Needless to say, **all commands** will be executed through the terminal.

# Prerequisites

To run this worklow you will surely need to install `conda` or `mamba`, a tool for package and environment management. When visiting [conda's webpage](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html) you will see that there are different installers available. Anyone will do. When running this workflow, I have been working with `miniconda`. You will find the instructions for installing it [here](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).

Once `conda` is installed (make sure to follow the proper instructions for your OS) and activated your prompt in the terminal should look like this:

```{r terminal_prompt, eval=FALSE}
(base) your_username@host_name: ~$ 
```  

Or this:
```{r terminal_prompt2, eval=FALSE}
(base) [your_username@host_name your_current_location]$
```

The `(base)` part indicates that `conda` is active and running its default environment. If for whatever reason it does not appear in your prompt, typing `source ~/.bashrc` or `source ~/.bash_profile` will do. A quick `ls -la ~` should tell you which type of `.bash` file you have. 

# Step 1 - Setting the proper environment

Okay! Now that you have successfully installed `conda` in your computer, we will focus on creating an environment that will ensure that the workflow runs properly. Based on the [original workflow](https://github.com/sdind/genome_annotation_workflow?tab=readme-ov-file), the following programs are required to run the workflow:

* conda >= v23.7.3 (already installed)
* singularity >= v3.7.3
* snakemake >= v7.32.3

I have prepared a `yaml` that contains all the required programs and dependencies to run the workflow. Get the `snakemake.yaml` file to your home directory and type the following command:

```{r install_snakemake, eval=FALSE}
conda env create -f snakemake.yaml
```

Once the environment is created, type `conda activate snakemake` and your prompt in the terminal should look like this:

```{r prompt_env, eval=FALSE}
(snakemake) [your_username@host_name your_current_location]$
```

# Step 2 - Get the programs

Now that we have our environment ready and activated, we will clone the github repository that contains the workflow. To do so, type: 

```{r git_clone, eval=FALSE}
git clone https://github.com/sdind/genome_annotation_workflow
```

You will see that a `/genome_annotation_workflow` folder is created in your current location. To know more about the directory structure within it, please check the [original repository](https://github.com/sdind/genome_annotation_workflow). If you browse around the inner folders, you will find a `snakemake_annotation.run` file in `/workflow_ann`, which is a `bash` file that sets up the analysis and the file that we will submit to the queue when working on a cluster. 

This workflow does pretty much everything through a series of rules implemented in Snakemake and should work without problems from here. However, I have encountered some incompatibilities that have messed up my analysis, specifically with `multiqc.yaml` in the `workflow_ann/workflow/envs` directory and some versions of certain dependencies included in other `*.yaml`. So, to ensure that everything runs smoothly, you should copy the `/envs` and `/rules` folders from this page and substitute the native one's in your `/workflow` directory. Essentially, I have disabled the `FASTQC` and `multiQC` steps, which recursively crashed my runs. However, it might not be the case for you. 

# Step 3 - Prepare your input data

This workflow performs the structural annotation by integrating RNA-seq and protein data. To run it, you will need the following:

* **Your** genome assembly in `fasta` format
* RNAseq data in `fastq` format
* A protein database in `fasta` format

To create a protein database you could download any of the publicly available at [OrthoDB](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/) that is closest to your taxon's genome or use a concatenated fasta of proteomes from closely related species. Or best, combine both. However, always make of two things: 

1.  Always make sure that the protein database is **uncompressed**, otherwise it will crash `BRAKER`. If you download any of OrthoDB's databases, i.e. the metazoan one, use `wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Metazoa.fa.gz` to obtain it. You can obtain the link by right-clicking on the hyperlink on the database in [OrthoDB's webpage](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/). Then, check that it has downloaded correctly comparing the `md5sum` number of your download with the one provided by the webpage, and unpack it with `tar -xvfz Metazoa.fa.gz`. Afterwards, you can add other proteomes to the metazoan backbone using a simple `cat`command.

2.  To prevent issues with some of the programs in this workflow, it is **best to not have long fasta headers** containing spaces or pipes (`|`). Use the following command to simplify your headers:
```{r simplify_comm, eval=FALSE}
sed '/^>/s/ .*//' input.fasta > output.fasta
```

# Step 4 - Configure the workflow 

We are almost done. Now we need to edit a couple of files to indicate where the input data is laying around within our `/workflow_ann` directory. If you have followed this guide, then there are only 2 files that we have to concern about. 

1.  The `config.yaml` file in `/workflow_ann/config` directory.
2.  The `Snakefile` file in the `/Workflow_ann/worflow` directory.

`config.yaml` contains two blocks of text. And it looks as follows: 

```{r config.yaml, eval=FALSE}

#parameters

asm: '/home/MEGddols/genome_annotation_workflow/workflow_ann/genome/Biomphalaria_glabrata.renamed.def.fasta'
snakemake_dir_path: '/home/MEGddols/genome_annotation_workflow/workflow_ann'
name: 'Biom_glabrata_v1'
busco_phylum: 'metazoa_odb10'
prot: '/home/MEGddols/genome_annotation_workflow/workflow_ann/Metazoa_odb/MetaNudis.fasta'
samples:
  biom_body:
    type: 'paired-end'
    R1: '/home/MEGddols/genome_annotation_workflow/workflow_ann/genome/ERR9682492_1.fastq'
    R2: '/home/MEGddols/genome_annotation_workflow/workflow_ann/genome/ERR9682492_2.fastq'




################################################################################
##################################  RESOURCES  #################################
################################################################################


# Simplified resources specification using presets. Threads are always used
# by Snakemake whenever possible, but 'mem_mb' and 'runtime' are only used
# when executing Snakemake on a cluster (i.e. using an execution profile) or
# a cloud environment (e.g. kubernetes..)
presets:
    threads:                  ###### Number of threads
        small: 1              #
        medium: 8             #
        large: 20              #
        very_large: 40        #
    mem_mb:                   ###### Memory in Mb
        very_small: 500
        small: 5000           #
        medium: 50000          #
        large: 200000          #
        very_large: 500000    #
    runtime:                  ###### Maximum runtime (format: D-HH:MM:SS)
        short: '0-05:00:00'   #
        medium: '0-15:00:00'  #
        long: '2-00:00:00'    #
        very_long: '3-00:00:00'

# File describing resources to use for each rule.
# Only edit if you know what you are doing
resources_info: config/resources.yaml

# Resources can be manually specified for each rule and will override the
# resources file specifications.
#
# resources:
#     busco_run:
#         threads: 32
#         mem_mb: 32000
resources:
```

It is best not to meddle with the `resources` block. However, the `#parameters` block is the one we are truly interested in. In this block we will have to indicate the paths to several input files to run the workflow. The paths can be relative or absolute, I prefer the second approach. The parameters are the following:

* `asm` - Path to the genome we want to annotate.
* `snakemake_dir_path` - Path to the directory containing the `Snakefile` file.
* `name` - Prefix that will be used to name output files generated throughout the workflow.
* `busco_phylum` - Database against which we want to conduct a BUSCO search to assess the completeness of our genome assembly. 
* `prot` - Path to the proteins database file.
* `samples` - Path to the RNA-seq files. In here we need to specify a name under which a series of datasets will be associated. Why? Because maybe we are using transcriptomes from different tissues to annotate our genome and it can be easier to track those files (i.e. transcript maps against the genome) downstream under specific names. We will have to indicate the `type` (often `paired-end`) and the paths to the `*_1.fastq` & `*_2.fastq` files. In the example above we are only using whole-body transcripts, but if we wanted to add data from specific body parts that we had sequenced of retrieved we could say so as follows: 

```{r config.yaml 2, eval=FALSE}
(...)
samples:
  biom_body:
    type: 'paired-end'
    R1: '/home/MEGddols/genome_annotation_workflow/workflow_ann/genome/ERR9682492_1.fastq'
    R2: '/home/MEGddols/genome_annotation_workflow/workflow_ann/genome/ERR9682492_2.fastq'
  rhinophores:
    type: 'paired-end'
    R1: '/path/to/folder/containing/the/reads_1.fastq'
    R2: '/path/to/folder/containing/the/reads_2.fastq'
  (...)
```

As per the `Snakefile`, we will have to make sure that all paths to the `config.yaml` and different rules (i.e. files ended with `.smk`) are correctly indicated. The file looks like this: 

```{r Snakefile, eval=FALSE}

configfile: '/home/MEGddols/genome_annotation_workflow/workflow_ann/config/config.yaml' 
include: '/home/MEGddols/genome_annotation_workflow/workflow_ann/workflow/rules/resources.smk' 
include: '/home/MEGddols/genome_annotation_workflow/workflow_ann/workflow/rules/1_MaskRepeat.smk'    
include: '/home/MEGddols/genome_annotation_workflow/workflow_ann/workflow/rules/2_alignRNA.smk'      
include: '/home/MEGddols/genome_annotation_workflow/workflow_ann/workflow/rules/3_braker.smk'        

rule all:
    input:
        directory(os.path.join(config['snakemake_dir_path'], 'results/2_braker/braker_busco')),
        os.path.join(config['snakemake_dir_path'], "results/2_braker/QC_RNA/multiqc_report.html"),
        expand(os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/hisat2/mapping_stats_samtools/{sample}_mapping_stats_samtools.txt'), sample=config['samples'].keys()),
        expand(os.path.join(config['snakemake_dir_path'], 'results/2_braker/align_RNA/hisat2/mapping_stats_qualimap_bamqc/{sample}_mapping_stats_qualimap_bamqc.pdf'), sample=config['samples'].keys())
```

Once all of this has been taken care of, we can move to the `/workflow_ann` directory where the `snakemake_annotation.run` file is and submit the job with `sbatch snakemake_annotation.run`. The workflow is now running and should take a few days to fully complete without errors. All result files will be at `/workflow_ann/results`. 
