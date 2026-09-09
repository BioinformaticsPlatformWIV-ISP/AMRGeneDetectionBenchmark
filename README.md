# AMRGeneDetectionBenchmark

This framework can be used to benchmark gene performance directly from long-read datasets. Although it was initially tested for the detection of AMR genes, it can also be used with other databases. It provides an entire codebase dedicated to:
1. Perform gene detection on long-read data using various tools;
2. Generate HTML reports summarizing the main findings, for either a single sample (sample reporter) or across all samples (aggregated reporter).

Version: **1.0**

Test files, as well as the databases, ground truths and genes covered are available in [Zenodo](10.5281/zenodo.18458145). 

## Installation

### Installing the framework

The current codebase runs on python 3.12, and all python dependencies are listed in the associated requirements.txt file, but the environment is installable via pip as follows:

```{bash}
git clone https://github.com/BioinformaticsPlatformWIV-ISP/AMRGeneDetectionBenchmark.git
cd AMRGeneDetectionBenchmark/
python3.12 -m venv amr_gene_detection_venv
source amr_gene_detection_venv/bin/activate
pip install .
```

A config file is required to run the workflow. A sample file is available in `detection/config/config.yml.sample`. Be sure to check and modify the values for your system.

The critical entries to modify are:
- dir_temp and dir_logs, the directories for respectively saving the intermediate and log files
- the database directory (see section [Setting up the database](#setting-up-the-database))
- your conda installation directory
- the name of the conda environments setup for each of the tools (see section [Installing the tools](#installing-the-tools))

```{bash}
cp {PATH_TO_INSTALLATION}/detection/config/config.yml.sample {PATH_TO_INSTALLATION}/detection/config/config.yml
vim {PATH_TO_INSTALLATION}/detection/config/config.yml
```

### Installing the tools

Multiple tools and executables need to be available in your PATH, for which precompiled binaries are available:

- [usearch v11.0.667](https://github.com/rcedgar/usearch_old_binaries/) 
- [seqkit v2.8.2](https://github.com/shenwei356/seqkit/releases/tag/v2.8.2) 
- [KMA v.14.12a](https://github.com/genomicepidemiology/kma), which could also be available as a conda environment

Additionally, the tools implemented in this framework require:
1. the internal NDARO database used for the benchmarking (available on the Zenodo archive);
2. to be available via separate conda environments. Due to dependencies incompatibilities, the tools cannot be bundled in a single environment.

- [Argo v0.2.0](https://github.com/xinehc/argo?tab=readme-ov-file#installation)
- [ARGpore2 v2.1](https://github.com/sustc-xylab/ARGpore2?tab=readme-ov-file#pre-requisites-for-argpore) (see also additional notes)
- [DeepARG v1.0.2](https://github.com/gaarangoa/deeparg?tab=readme-ov-file#use-conda-environment)
- [ShortBRED v0.9.5](https://github.com/biobakery/biobakery/wiki/shortbred#1-install)

### Additional note for ARGpore2

To install ARGpore2, create a conda environment containing the dependencies required as explained in the [associated GitHub repository](https://github.com/sustc-xylab/ARGpore2), then follow-up with the installation. The argpore executable should then be in your PATH, which can be performed with the following command:

```{bash}
export PATH=$PATH:/path/to/argpore_dir/
```

Afterward, you'll have to replace manually the content of the database directory with the one available on the Zenodo repository. The command used would then be:

```{bash}
cp /your/downloads/databases/ARGPORE/* /path/to/argpore_dir/database/
```

### Additional note for ShortBRED

For ShortBRED, the following command was used to create the required conda environment, using mamba:

```{bash}
mamba create -n shortbred_env biobakery::shortbred muscle cd-hit blast
```
Importantly, the usearch version should be below v12.0 as it causes incompatibilities with ShortBRED. The correct usearch version used here (v11.0.667) can be found [here](https://github.com/rcedgar/usearch_old_binaries/). The binary should be renamed `'usearch` and should be added to your PATH variable. To realize this, download the usearch binary and place it in a specific folder. An example command would look like: 

```{bash}
cp /path/to/downloads/usearch11.0.667_i86linux64 /path/to/usearch_dir/usearch
export PATH=$PATH:/path/to/usearch_dir/
```

### Setting up the database

#### Bundled NCBI database

The bundled NCBI database for each tool is available in the Zenodo repository. Copy the directory to a location of your choosing, then this directory should be added to the config file, under the `database_directory` entry.

#### Setting up a custom database

1. ARGpore2
   - To build the custom database, the faster way is to first rename the fasta file of your choice as `SARG_20211207_14210_filter.ffn`;
   - Then, index this fasta file of your choice using the provided `last` executable (i.e., in the `bin` directory of the ARGpore2 project);
   - Finally, the indexed fasta file and index files should be copied in the database directory of the ARGPORE2 project
2. DeepARG
   - First, the sequences of the database of your choice must be translated from nucleotides to amino acid sequences;
   - Then, the procedure to create a new model is explained in the associated wiki [here](https://bitbucket.org/gusphdproj/deeparg-largerepo/src/master/), under the section "Train the DeepARG models".
3. KMA
   - KMA is the easiest to setup, as the nucleotides FASTA file only needs to be indexed by 'kma index'.
4. ShortBRED
   - First, the nucleotide sequences must be translated into amino-acid sequences;
   - Afterwards, the procedure is explained in the ShortBRED wiki [here](https://github.com/biobakery/shortbred?tab=readme-ov-file#running-shortbred-identify-to-create-new-markers). 
5. Argo
   - Finally, for Argo, there are only three databases that can be used, as the setup scripts are fully automated: NDARO, CARD and SARG+ (which is the default).
   - The procedure is documented [here](https://github.com/xinehc/argo-supplementary).

## Running the framework

To run the gene detection framework on a single dataset, the command `amr_detection` can be launched with the following options, where the `--fastq-se` and `--sample-name` arguments are mandatory.

```
amr_detection --help
usage: amr_detection [-h] --fastq-se FASTQ_SE [--sample-name SAMPLE_NAME] [--dir-out DIR_OUT] [--dir-working DIR_WORKING] [--mode {classification,gene_detection}] [--min-identity] [--min-coverage]
                             [--keys-list KEYS_LIST] [--keys KEYS] [--list] [--dataset-list] [--nb-jobs NB_JOBS] [--threads-per-job THREADS_PER_JOB] [--debug]

options:
  -h, --help            show this help message and exit
  --fastq-se FASTQ_SE   Input ONT FASTQ file
  --sample-name SAMPLE_NAME
                        Sample name
  --dir-out DIR_OUT     Directory to save results
  --dir-working DIR_WORKING
                        Working directory
  --min-identity        Minimum identity percentage [0-100]
  --min-coverage        Minimum length coverage [0-100]
  --keys-list KEYS_LIST
                        TXT file containing a list of tools to run
  --keys KEYS           Limit to the tools with these keys (comma separated list)
  --list                Lists the available tools, does NOT run
  --dataset-list        Lists the available datasets, does NOT run
  --nb-jobs NB_JOBS     Number of jobs to start in parallel
  --threads-per-job THREADS_PER_JOB
                        Number of threads per job
  --debug               Run in debug mode, does not delete intermediary files
```
The `--keys` options can be used to run only the tools you are interested in.

The first step in the workflow trims the reads are trimmed using [SeqKit](https://github.com/shenwei356/seqkit) command `seqkit seq` for a minimum length of 1,000 bases and minimum read quality of 7. Statistics are then computed using the `seqkit stats` command.

Afterward, the tools are ran separately on the trimmed reads.

### Sample reporter

The command `reporter_sample` can then be used to generate an HTML report of the results generated previously.

The mandatory options are 
- `--dir-in`, corresponding to the sample directory generated by the amr_detection command
- `--configfile` corresponding to the generated config file, which is available in the same directory
- `--coverage-amr-genes` and `--blast-ground-truth`, available in the Zenodo archive
- `--sample-name`, corresponding to the sample name given previously

The default output directory is named `reporter_out`. Feel free to change this value, but all sample reporters should be stored in the same output directory.

#### Interface & example command

```
reporter_sample --help
options:
  -h, --help            show this help message and exit
  --dir-in DIR_IN       Input directory
  --dir-working DIR_WORKING
                        Working directory
  --dir-out DIR_OUT     Output directory
  --sample-name SAMPLE_NAME
                        Sample name
  --ground-truth GROUND_TRUTH
                        Ground truth
  --configfile CONFIGFILE
                        config file used by the main detection
  --coverage-amr-genes COVERAGE_AMR_GENES
                        Coverage amr genes file
  --blast-ground-truth  tick if you have blast-based ground truth [True]
  --convert-gene-name   convert gene allele to gene name [True]

```

To run a reporter on previously analyzed datasets, the input directory is the sample directory created by the main detection script. 

The HTML report will be saved in `[dir-out]/[sample-name]/reporter_gene/`.

#### Output

```
reporter_gene/
|- [Tool]/                            # One subfolder for each tool, containing the PNG of each plot
|- [Tool]_for_comparison.csv          # Characteristics for each AMR gene detected
|- [Tool]_precision_recall_scores.csv # When available, corresponds to precision and recall values for each coverage and identity threshold evaluated
|- per_sample_report.html             # HTML report containing the main figures and results
|- presence_absence.png               # Figure of the AMR genes present in the ground truth
|- statistics.png                     # Barplot figure of the precision, recall and F1 score values
|- statistics_summary.tsv             # For each tool, summarizes the main metrics and AMR genes detected
```

### Aggregated reporter

Finally, after running the reporter for multiple datasets, the main results can be aggregated into a single report using the command `reporter_aggregated`.

The mandatory options are :
- `--dir-in`, corresponding to the output directory from the reporter_sample
- `--dir-out`, which is the output directory that will contain the output result. 

You can specify the list of tools you'd like to have the aggregated report for as well.

#### Interface & example command

```
reporter_aggregated --help
usage: reporter_aggregated [-h] -i DIR_IN -o DIR_OUT [-t {allele,gene}] [-l {Argo,ARGpore2,KMA,DeepARG,ShortBRED} [{Argo,ARGpore2,KMA,DeepARG,ShortBRED} ...]] [-s {yield,downsampling}]

options:
  -h, --help            show this help message and exit
  -i DIR_IN, --dir-in DIR_IN
                        input directory with reporter results
  -o DIR_OUT, --dir-out DIR_OUT
                        output directory with output figures
  -t {allele,gene}, --type {allele,gene}
                        type of reporter
  -l {Argo,ARGpore2,KMA,DeepARG,ShortBRED} [{Argo,ARGpore2,KMA,DeepARG,ShortBRED} ...], --list-of-tools {Argo,ARGpore2,KMA,DeepARG,ShortBRED} [{Argo,ARGpore2,KMA,DeepARG,ShortBRED} ...]
                        selection of tools to include, default=all
```

To run the reporter for obtaining the aggregated results, the `--dir-in` option corresponds to the output directory of the main script.

Of note, you can generate aggregated reports for only a subset of tools by selecting one or more tools with the option `--list-of-tools` and separating them by commas.

#### Output

```
aggregated_reporter/
|- aggregated_auprc_[Sample].png          # Plots of the aggregated AUPRC values for each sample
|- aggregated_boxplot.png                 # Boxplot of the aggregated values across all samples
|- aggregated_mean_values_[identity%].png # Plots showing the mean metric values across all samples analyzed
|- report.html                            # HTML report containing the main figures and results
|- auprc.tsv                              # Summary of the AUPRC values for each coverage and identtity threshold
|- mean_values.csv                        # Summary of the mean metric values represented in the mean_values plots
```


## Tutorial

This tutorial will show how to run the entire framework using two test datasets and the in-house NDARO database. 
Under "test_data" in the zenodo archive, two downsampled data for the Zymo D6300 and Zymo D6331 R10 datasets are available.

### 1. Run the AMR gene detection

```{bash}
amr_detection --fastq-se test_data/zymo_test.fq.gz --sample-name zymo --keys argo,kma
amr_detection --fastq-se test_data/zymo_gut_test.fq.gz --sample-name zymo_gut_r10 --keys argo,kma
```

After completion, the output directories will contain:

```{bash}
main_out/
|-zymo/
    |- argo/
        |- argo_cleaned.json
    |- config.yml
    |- input/
        |- fastq_filtered.fastq.gz
        |- fastq_filtered_stats.tsv
        |- input.fastq.gz
    |- kma/
        |- kma_cleaned.json
        |- kma_for_pr.tsv
|-zymo_gut_r10/
    |- [same structure as zymo]
```

An overview of the JSON output file will look like (truncated): 

```{bash}
head -n30 argo_cleaned.json
```

```{json}
[
  {
    "lineage": "Bacteria;Bacillota;Bacilli;Lactobacillales;Listeriaceae;Listeria;Listeria monocytogenes_B",
    "type": "lincosamide",
    "reference_accession": "vga(G)",
    "coverage_depth": 21.683,
    "genome": 24.25,
    "abundance": 0.894,
    "sequence_identity": 0.9,
    "coverage_percentage": 0.9,
    "input_sequence_id": "sequence_id"
  },
  {
    "lineage": "Bacteria;Bacillota;Bacilli;Staphylococcales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus",
    "type": "beta-lactam",
    "reference_accession": "blaI",
    "coverage_depth": 26.254,
    "genome": 19.375,
    "abundance": 1.355,
    "sequence_identity": 0.9,
    "coverage_percentage": 0.9,
    "input_sequence_id": "sequence_id"
  },
```

### 2. Run the reporter

```{bash}
reporter_sample --dir-in main_out/zymo --sample-name zymo --configfile main_out/zymo/config.yml
--ground-truth [zenodo_location]/ground_truths/ground_truth_zymo.tsv --coverage-amr-genes [zenodo_location]/amr_genes_covered/zymo_genes_covered.tsv

reporter_sample --dir-in main_out/zymo_gut_r10 --sample-name zymo_gut_r10 --configfile main_out/zymo_gut_r10/config.yml
--ground-truth [zenodo_location]/ground_truths/ground_truth_zymo_gut.tsv --coverage-amr-genes [zenodo_location]/amr_genes_covered/zymo_gut_r10_genes_covered.tsv
```

After completion -> 

```
ls reporter_out
|-zymo/
    |- reporter_gene/
        |- Argo/
        |- KMA/
        |- report/
        |- statistics_summary.tsv
        |- updated_ground_truth.tsv
```

Where the generated HTML report is available under `report/`.

An overview of the `statistics_summary.tsv` file, which contains all statistics used to generate the figures in the report, looks like (truncated):

```{bash}
Argo    identity_filter 90
Argo    length_filter   80
Argo    GT      blaI,vmlR,blaZ,fosX,catB7,tet(L),blaPDC,blaOXA,tet(38),vga(G),aadD1,fosB,emrD,crpP,lsa(A),mexA,aac(6'),mepA,blaR1,blaEC
Argo    TP      blaI,vmlR,blaZ,fosX,catB7,tet(L),blaPDC,blaOXA,tet(38),vga(G),aadD1,fosB,emrD,crpP,mexA,aac(6'),mepA,blaR1,blaEC
Argo    FN      lsa(A)
Argo    FP      dfrC,fosM1,fosM2,mphK,toprJ,sdeY,dfrA51
Argo    Precision       0.7307692307692307
Argo    Recall  0.95
Argo    F1 score        0.8260869565217391
```

### 3. Run the reporter to aggregate the results

```{bash}
reporter_aggregated --dir-in reporter_out --dir-out aggregated_reporter -l Argo,KMA
```

After completion:

```{bash}
ls aggregated_reporter
|- aggregated_boxplot_stats.csv
|- auprc.tsv
|- mean_values.csv
|- report.html
|- figures/
```

## Help & Feedback

Feel free to [open an issue](https://github.com/BioinformaticsPlatformWIV-ISP/AMRGeneDetectionBenchmark/issues) to discuss improvements or ask for help and feedback.

## Citation

In the corresponding publication, the framework was used to analyze defined mock community (DMC) datasets, available in the associated (paper)[].
