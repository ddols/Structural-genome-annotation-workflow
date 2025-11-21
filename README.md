[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://ddols.github.io/Structural-genome-annotation-workflow/)

Structural genome annotation workflow
================

- [Before we get started](#before-we-get-started)
- [Prerequisites](#prerequisites)
- [Step 1 - Setting the proper
  environment](#step-1---setting-the-proper-environment)
- [Step 2 - Get the programs](#step-2---get-the-programs)
- [Step 3 - Prepare your input data](#step-3---prepare-your-input-data)
- [Step 4 - Configure the workflow](#step-4---configure-the-workflow)
- [Step 5 - Obtaining longest isoforms for downstream analyses](#step-5---obtaining-longest-isoforms-for-downstream-analyses)

This is a step-by-step guide explaining some of the steps in [/sdind’s
genome_annotation_workflow](https://github.com/sdind/genome_annotation_workflow?tab=readme-ov-file). All credit goes to their work. Be sure to
visit their repository for an in-depth dive into the workflow. This
guide is meant for people who are new to structural genome annotation.

# Before we get started

This tutorial is meant to showcase one of the multiple ways a genome may
be structurally annotated. That is, basically, to physically locate
every possible genomic feature such as exons, introns, genes and regulatory
sequences (whenever it is possible). To do so, it is fundamental to not
only have a genome to annotate (duh) but, most importantly, transcriptomic
data, which will be essential to characterise the previously mentioned
features.

In an ideal world, the transcriptome should belong to the same species
whose genome we are trying to annotate. However, in times of
need, a sufficiently close species’ transcriptome should do the trick.
And in more desperate times, well, any transcriptome more or less
related to our taxon is acceptable. Anyway, this workflow is an
adaptation of
[this](https://github.com/sdind/genome_annotation_workflow?tab=readme-ov-file)
workflow by **sdind**.

Lastly, this workflow is written with `Linux` in mind as the running OS
with access to a cluster with a
[Slurm](https://slurm.schedmd.com/quickstart.html) system. If you are
working on a
[SGE](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html) cluster
or any other kind, you can still run the scripts in this repository,
although you will have to modify them accordingly for your queueing
system. Needless to say, **all commands** will be executed through the
terminal.

# Prerequisites

To run this workflow, you will surely need to install `conda` or `mamba`,
a tool for package and environment management. When visiting [conda’s
webpage](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html),
you will see that there are different installers available. Anyone will
do. When running this workflow, I have been working with `miniconda`.
You will find the instructions for installing it
[here](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).

Once `conda` is installed (make sure to follow the proper instructions
for your OS) and activated, your prompt in the terminal should look like
this:

``` r
(base) your_username@host_name: ~$ 
```

Or this:

``` r
(base) [your_username@host_name your_current_location]$
```

The `(base)` part indicates that `conda` is active and running its
default environment. If, for whatever reason, it does not appear in your
prompt, typing `source ~/.bashrc` or `source ~/.bash_profile` will do. A
quick `ls -la ~` should tell you which type of `.bash` file you have.

# Step 1 - Setting the proper environment

Okay! Now that you have successfully installed `conda` on your computer,
let's focus on creating an environment that will ensure that the
workflow runs properly. Based on the [original
workflow](https://github.com/sdind/genome_annotation_workflow?tab=readme-ov-file),
the following programs are required to run the workflow:

- conda \>= v23.7.3 (already installed)
- singularity \>= v3.7.3
- snakemake \>= v7.32.3

I have prepared a `yaml` that contains all the required programs and
dependencies to run the workflow. Get the `snakemake.yaml` file to your
home directory and type the following command:

``` r
conda env create -f snakemake.yaml
```

Once the environment is created, type `conda activate snakemake` and
your prompt in the terminal should look like this:

``` r
(snakemake) [your_username@host_name your_current_location]$
```

# Step 2 - Get the programs

Now that we have our environment ready and activated, we will clone the
github repository that contains the workflow. To do so, type:

``` r
git clone https://github.com/sdind/genome_annotation_workflow
```

You will see that a `/genome_annotation_workflow` folder is created in
your current location. To know more about the directory structure within
it, please check the [original
repository](https://github.com/sdind/genome_annotation_workflow). If you
browse around the inner folders, you will find a
`snakemake_annotation.run` file in `/workflow_ann`, which is a `bash`
file that sets up the analysis and the file that we will submit to the
queue when working on a cluster.

This workflow does pretty much everything through a series of rules
implemented in Snakemake and should work without problems from here.
However, I have encountered some incompatibilities that have messed up
my analysis, specifically with `multiqc.yaml` in the
`workflow_ann/workflow/envs` directory and some versions of certain
dependencies included in other `*.yaml`. So, to ensure that everything
runs smoothly, you should copy the `/envs` and `/rules` folders from
this page and substitute the native ones in your `/workflow` directory.
Essentially, I have disabled the `FASTQC` and `multiQC` steps, which
recursively crashed my runs. However, it might not be the case for you.

# Step 3 - Prepare your input data

This workflow performs the structural annotation by integrating RNA-seq
and protein data. To run it, you will need the following:

- **Your** genome assembly in `fasta` format
- RNAseq data in `fastq` format
- A protein database in `fasta` format

To create a protein database, you could download any of the publicly
available at
[OrthoDB](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/)
that is closest to your taxon’s genome, or use a concatenated fasta of
proteomes from closely related species. Or best, combine both. However,
to avoid unsolicited problems, always check these two things:

1.  Always make sure that the protein database is **uncompressed**,
    otherwise it will crash `BRAKER`. If you download any of OrthoDB’s
    databases, i.e. the metazoan one, use
    `wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Metazoa.fa.gz`
    to obtain it. You can obtain the link by right-clicking on the
    hyperlink on the database in [OrthoDB’s
    webpage](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/).
    Then, check that it has downloaded correctly by comparing the `md5sum`
    numbers of your download with the one provided by the webpage, and
    unpack it with `tar -xvfz Metazoa.fa.gz`. Afterwards, you can add
    other proteomes to the metazoan backbone using a simple
    `cat`command.

2.  To prevent issues with some of the programs in this workflow, it is
    **best to not have long fasta headers** containing spaces or pipes
    (`|`). Use the following command to simplify your headers:

``` r
sed '/^>/s/ .*//' input.fasta > output.fasta
```

# Step 4 - Configure the workflow

We are almost done. Now we need to edit a couple of files to indicate
where the input data is within our `/workflow_ann`
directory. If you have followed this guide, then there are only 2 files
that we have to be concerned about.

1.  The `config.yaml` file in `/workflow_ann/config` directory.
2.  The `Snakefile` file in the `/workflow_ann/workflow` directory.

`config.yaml` contains two blocks of text. And it looks as follows:

``` r
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
# when executing Snakemake on a cluster (i.e., using an execution profile) or
# a cloud environment (e.g., kubernetes..)
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

It is best not to meddle with the `resources` block. However, the
`#parameters` block is the one we are truly interested in. In this block,
we will have to indicate the paths to several input files to run the
workflow. The paths can be relative or absolute. I prefer the second
approach. The parameters are the following:

- `asm` - Path to the genome we want to annotate.
- `snakemake_dir_path` - Path to the directory containing the
  `Snakefile` file.
- `name` - Prefix that will be used to name output files generated
  throughout the workflow.
- `busco_phylum` - Database against which we want to conduct a BUSCO
  search to assess the completeness of our genome assembly.
- `prot` - Path to the proteins database file.
- `samples` - Path to the RNA-seq files. In here, we need to specify a
  name under which a series of datasets will be associated. Why? Because
  maybe we are using transcriptomes from different tissues to annotate
  our genome, and it can be easier to track those files (i.e., transcript
  maps against the genome) downstream under specific names. We will have
  to indicate the `type` (often `paired-end`) and the paths to the
  `*_1.fastq` & `*_2.fastq` files. In the example above, we are only
  using whole-body transcripts, but if we wanted to add data from
  specific body parts that we had sequenced or downloaded, we could say so
  as follows:

``` r
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

As per the `Snakefile`, we will have to make sure that all paths to the
`config.yaml` and different rules (i.e., files ended with `.smk`) are
correctly indicated. The file looks like this:

``` r
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

Once all of this has been taken care of, we can move to the
`/workflow_ann` directory where the `snakemake_annotation.run` file is,
and submit the job with `sbatch snakemake_annotation.run`. The workflow
is now running and should take a few days to fully complete without
errors. All result files will be in `/workflow_ann/results`.

# Step 5 - Obtaining longest isoforms for downstream analyses

Eventually, the workflow will happily end (fingers crossed), and then we will be
able to move on with either the functional annotation of the genome or focus on whatever downstream 
analyses we feel like. If you browse the `workflow_ann/results/2_braker/out_braker/braker`
directory, you will find a compendium of `BRAKER`'s outputs (you can learn more about them 
in [here](https://github.com/Gaius-Augustus/BRAKER), at the section **Output of BRAKER**).

Essentially, both `braker.gtf` and `braker.gff3` are the final gene set of our genome assembly, 
`braker.aa` is a multifasta containing all protein sequences (hence, amino-acid sequences) characterised thanks to
our RNA-seq data, and `braker.condingseq`, which is a multifasta containing the nucleotide sequences of 
the genome's gene set. These files contain all possible isoforms identified during the annotation.
However, in some cases, we are only interested in working with the longest ones to perform [macrosynteny](https://github.com/jtlovell/GENESPACE)
analysis, determine sets of [orthologous genes](https://github.com/OrthoFinder/OrthoFinder), or [phylogenomics](https://github.com/smirarab/ASTRAL)
to name a few examples. 

To extract the longest isoforms, we will need at least the following programs:

- [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA)
- [Augustus](https://github.com/Gaius-Augustus/Augustus)
- `biopython`

### Step 5.1 - Installing the required programs

We will start off by creating a new environment for these programs. You can name it however you like.
For this guide's sake we will just name it **tsebra**. To create it, just type:

```r
conda create -n tsebra
```
Followed by:

```r
conda activate tsebra
```

Then, we will [install TSEBRA](https://bioconda.github.io/recipes/tsebra/README.html) following the instructions, 
and do the same for [Augustus](https://anaconda.org/bioconda/augustus) (you can also try following [this](https://github.com/Gaius-Augustus/Augustus?tab=readme-ov-file#building-augustus-from-source) indications), 
and `biopython` (`pip3 install biopython`). You may have noticed that to install `tsebra` bioconda recommends using `mamba`.
In case you do not know where `mamba` is, fear not, you will find it within the `/bin` directory of your `conda` installer
folder. For example, if you chose `miniconda3` as your installer, to install `tsebra` you should type:

```r
~/path-to/miniconda3/bin/mamba install tsebra
```

### Step 5.2 - Running the scripts

Everything is ready now to extract the longest isoforms from `BRAKER`'s outputs. Within this guide's directory `/longest_isoforms_scripts`,
you will find two scripts numbered according to their intended order of usage. Remember that you will need to modify them
a bit depending on whether you are working locally, on a SGE cluster, or a Slurm one.

The first one `1-get_longest_isoforms.sh` will call a script from `TSEBRA` that will parse the longest isoforms
to a new gtf file ended with `*_longest.gtf`. To run it, you will only need to place the script in `/results/2_braker/out_braker/braker` folder, be
on **tsebra**'s `conda` environment, and indicate the base name of the `braker.gtf` file. Which is precisely `braker`. Of course, you can change that name beforehand with `mv` just like this:

```r
mv braker.gtf biomphalaria.gtf
```

Now the base name is **biomphalaria**. The script should be sent to queue as follows:

```r
sbatch 1-get_longest_isoforms.sh biomphalaria
```
A new file named `biomphalaria_longest.gtf` will have appeared in your directory. Then, we will run `2-get_fastas_of_long_isoforms.sh` in the same directory.
This script uses a script from `Augustus` that will use our masked genome assembly and the newly obtained `*_longest.gtf` file
to extract both the protein set and nucleotide set of the longest isoforms in `fasta` format. You will find the masked assembly in
`/results/1_MaskRepeat/RepeatMasker`, perfectly labelled. Now, to run it type:

```r
sbacth 2-get_fasta_of_long_isoforms.sh /path-to-your/fasta.masked your_longest.gtf
```

And, voilà! Two files will have appeared named `longest_isoforms.aa` and `longest_isoforms.codingseq`. Happy analyzing!







