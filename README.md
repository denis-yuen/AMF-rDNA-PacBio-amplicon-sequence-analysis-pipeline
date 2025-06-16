## ABOUT
This pipeline analyses PacBio AMF amplicon sequences using a workflow management system.  It has been designed to run with [Snakemake](https://snakemake.readthedocs.io/en/stable/) with [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) environments.


## SUMMARY
The workflow has five main steps: <br>
1) Primer removal, quality filtering and size selection using [Cutadapt](https://github.com/marcelm/cutadapt) (v4.5).
2) Dereplication and chimera removal using [VSEARCH](https://github.com/torognes/vsearch) (v2.24).
3) Clustering of dereplicated sequences using [Swarm](https://github.com/torognes/swarm) (v3.1.4).
4) Filter Swarm clusters based on taxonomy and abundance:
- Taxonomic filtering: the 20 most abundant Swarm clusters are locally aligned using BLAST (v2.16.0) against the core_nt database, and their taxonomy is determined using TaxonKit (v0.17.0). The SSU (~1.2 Kb) and Krüger fragment (end of SSU-ITS-partial LSU, ~1.5 Kb) are queried separately against the core_nt database due to the limited number of 2.7 kb AMF sequences in NCBI. BLAST results based on the Krüger fragments are used to filter out clusters that do not belong to the phylum Mucoromycota.
- Abundance filtering: taxonomically filtered Swarm clusters are filtered based on their abundance to reduce the likelihood that SNPs are due to technical artifacts rather than true biological variation. The filtering process has the following conditions:
  *  If there are fewer than 7 unique abundance values and all abundances are 1, all sequences are removed.
  *  If there are fewer than 7 unique abundance values but not all abundances are 1, only the sequence with the highest abundance is retained.
  *   If there are more than 7 unique abundance values, a 7-cluster K-means algorithm is applied:
       * If the maximum abundance is greater than 150, the top three groups are kept.
       * If the maximum abundance is less than 150, only the top two groups are kept.

5) Alignment: the filtered Swarm clusters are aligned with [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html) (v7.525)

#### Rationale for filtering conditions:

As both the number of Swarm clusters and their associated abundance values can vary considerably between samples, adaptive thresholding based on k-means clustering was used to identify significantly abundant Swarm clusters. The k-means clustering grouped the abundance values into statistical clusters based on their distribution.

The number of k-means clusters was set to seven based on the analysis of 11 PacBio datasets representing six Rhizophagus irregularis strains. For k=7, no more than 11 Swarm clusters per sample were recovered, which is consistent with the known number of rDNA copies in R. irregularis genomes.

For samples with Swarm cluster abundance ≤ 150, the distribution tends to be narrow with low abundance (i.e. noise) Swarm clusters. In this case, only the top 2 k-means clusters are retained to ensure a stricter, more conservative selection of high-confidence Swarm clusters.

For samples with Swarm cluster abundance > 150, the distribution tends to be wider with more high abundance Swarm clusters. In this case, the top 3 k-means clusters are selected to provide a more inclusive but still focused retention of likely biologically significant Swarm clusters.

## REQUIREMENTS

> NOTE 1 - Running on a Linux workstation<br>
You can run this pipeline on a Linux workstation with Conda installed. Before running, edit the ```conda_location``` in the [run_prepare_and_cluster_amplicons.sh](scripts/run_prepare_and_cluster_amplicons.sh) script to match the path to your conda installation. The required environments will be created in the *envs* folder the first time you run the pipeline. You will also need to have 
[snakemake](https://snakemake.readthedocs.io/en/stable/) installed. You can install it into your conda 'base' environment using ```conda install snakemake```.

> NOTE 2 - Running on a slurm cluster<br>
The first time the pipeline is launched, an internet connection is required to build the conda environments. Once these environments have been created, the pipeline can be used on computing nodes without internet.

> NOTE 3 - Input file formats<br>
This pipeline required files in fastq or fasta format  and can use gzipped files (.gz). For BAM files, use [Samtools](http://www.htsalllib.org) to convert them to fastq or fasta. 

1) Create a working directory and copy the files from this repository into it. Your working directory should look like this:

```
work_dir
├── README.md
├── config
    └── config_prepare_and_cluster_amplicons.yaml
└── workflow
	├── fastq_files
	│   ├── bc1013.fq.gz
	│   └── bc1026.fq.gz
	├── envs
	│   ├── blast.yml
	│   ├── mafft.yml
	│   └── snakemake_amplicon_clustering_v2.yml
	├── run_prepare_and_cluster_amplicons.sh
	├── run_prepare_and_cluster_amplicons_using_qsub.sh
	├── sample_list.txt
	├── scripts
	├── sequence.gi
	└── Snakefile_prepare_and_cluster_amplicons
```

2) Copy or create a symbolic link to the  PacBio sequences in the *fastq_files* folder.<br>

```bash
ln -s /path/to/raw/PacBio/*.fq.gz /path/to/fastq_files
```

3) Prepare a [sample_list.txt](workflow/sample_list.txt) file by splitting the sample names from the suffixes for the input files in the fastq_files directory. The first line of the [sample_list.txt](snakemake/sample_list.txt) file should have the word 'Sample' followed by one sample name per line. Examples of these files are available in the repository. 

```bash
cd /path/to/work_dir/workflow
ls fastq_files > sample_list.txt
sed -i 's/\.fastq//g' sample_list.txt
```

4) Provided here as an example, [sequence.gi](workflow/sequence.gi) is a list of GIs of all entries in the core_nt database labeled as uncultured and/or environmental from the Mucoromycota phylum. It is used here with blastn to exclude these GIs from blast results. If you do not need such a list, the line`-negative_gilist {sequence_gi} \` from both blatsn commands can be removed from *Snakefile_prepare_and_cluster_amplicons*.

```bash
blastn \
    -db {BLASTDB}/core_nt \
    -task megablast \
    -query {input} \
    -dust no \
    -evalue 1e-20 \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -negative_gilist {sequence_gi} \ # can be removed
    -query_loc 1-1200 \
    -outfmt "6 qseqid sseqid evalue qcovs qlen slen length pident staxids sscinames sskingdoms" \
    -num_threads {threads} \
    -out {output.blast_ssu}
```

5) The swarm clusters are searched against core_nt database using local BLAST. Instructions for obtaining the DB can be found at [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569850/). In addition, [TaxonKit](https://bioinf.shenwei.me/taxonkit/download/) is used to translate taxIDs from the blast results into full taxonomy ranks. TaxonKit requires a database, which can be installed using [these instructions](https://bioinf.shenwei.me/taxonkit/usage/).


## PARAMETERS
Edit the [config_prepare_and_cluster_amplicons.yaml](scripts/config_prepare_and_cluster_amplicons.yaml) file to match the names in the sample list file and the input file suffix. Also edit the forward and reverse primer sequences and their reverse complements if required as well as the min_length and max_length of amplicons to be retained. Finally, the GI list, BLASTDB and TAXONKIT_DB paths must be edited according to their respective locations.

```bash
### This is the config file for the Snakefile_prepare_amplicons ###
## You can edit it to change sample information and parameters as required ##

min_length: 2400
max_length: 3800

## Primer information ##
forward_primer: "ATCAACTTTCGATGGTAGGATAGA"
reverse_primer: "AACACTCGCAYAYATGYTAGA"
forward_primer_revcom: "TCTATCCTACCATCGAAAGTTGAT"
reverse_primer_revcom: "TCTARCATRTRTGCGAGTGTT"

### Sample information ###
## Split each input filename into sample name plus suffix
## NOTE - filenames must exactly match the sample name plus the suffix.
## Put the list of sample names in a text file
## one sample name per line and with the header 'Sample'

sample_file: "sample_list.txt"

## data for blastn
sequence_gi: "sequence.gi"
BLASTDB: "/path/to/DB/core_nt"
TAXONKIT_DB: "/path/to/DB/core_nt"

## Input file information ##
input_data_directory: "fastq_files"
input_file_suffix: ".fq.gz"


## Processing amplicons using cutadapt and vsearch ##
## Should not need to edit parameters below
processed_data_directory: "amplicons_processed"
processed_file_suffix: ".trimmed.combined.fasta" 

dereplicated_data_directory: "amplicons_dereplicated"
dereplicated_file_suffix: ".vsearch_dereplicated.fasta"

chimeras_removed_directory: "amplicons_chimeras_removed"
nonchimeras_file_suffix: ".nonchimeras.fasta"

## Clustering processed amplicons using swarm ##
swarm_clusters_directory: "clusters_swarm"
swarm_input_directory: "amplicons_chimeras_removed"
swarm_input_file_suffix: ".nonchimeras.fasta"

## TOP20
top20_clusters_directory: "clusters_top20"
top20_file_suffix: ".top20.swarm.seeds.fasta"

## Blast
blast_directory: "clusters_top20/clusters_blasted"
blast_file_suffix: ".taxonomy.tsv"

## Taxonomy results
taxonomy_directory: "clusters_top20/clusters_taxonomy"

## Filter taxonomy swarm clusters
tax_filtered_directory: "clusters_tax_filtered"
tax_filter_suffix: ".tax_filtered.swarm.seeds.fasta"

## Removing small clusters (singletons and Median Absolute Deviation (MAD) based)
filtclusters_directory: "clusters_abd_filtered"
filtclusters_file_suffix: ".abd_filtered.swarm.seeds.fasta"
kmeans_file_suffix: ".swarm.kmeans.tsv"

## Mafft alignment
aligned_directory: "clusters_aligned"
aligned_file_suffix: ".abd_filtered.swarm.seeds.txt"


### Software parameters ###
threads: 4

## Use these parameters
## Cutadapt for trimming, sorting and filtering amplicons 
## vsearch for dereplicating amplicons and removing chimeras
## Swarm for clustering: -d is resolution; -z is for vsearch format IDs;
## -f flag is for fastidious to combine small clusters with larger ones
cutadapt_parameters: "--no-indels --discard-untrimmed --quality-cutoff 50"
vsearch_dereplicate_parameters: "--derep_fulllength"
vsearch_remove_chimeras_parameters: "--abskew 7.0"
swarm_parameters: "-zf -d 1 -t 10"
```

## USAGE

Run the pipeline on a Linux workstation using [```run_prepare_and_cluster_amplicons.sh```](workflow/run_prepare_and_cluster_amplicons.sh). If you run it on a slurm cluster, use [```sbatch run_prepare_and_cluster_amplicons_using_sbtach.sh```](workflow/run_prepare_and_cluster_amplicons_using_sbatch.sh).


## OUTPUT
1) Processed amplicons from step one should be in the directory *amplicons_processed* including the output files named *<sample_name>.trimmed.combined.fasta.gz*

2) Dereplicated amplicons should be in the directories *amplicons_dereplicated* and *amplicons_chimeras_removed*

3) Results of running swarm clustering should be in the directory *clusters_directory* including output files named *<sample_name>.swarm.seeds.fasta*

4) Results after step four are located in three folders:
 - *clusters_top20* & *clusters_tax_filtered* are the output of the first phase where:
   +  *pSSU_ITS_pLSU_taxonomy.tsv* and *SSU_taxonomy.tsv* are taxonomy assignments of all samples
   +  *<sample_name>.top20.swarm.seeds.fasta* are the 20 most abundant clusters
 - *clusters_abd_filtered* is the output after the filtration based on the cluster abundances including:
   +  *<sample_name>.abd_filtered.swarm.seeds.fasta*, the clusters remaining after filtration
   +  *<sample_name>.swarm.kmeans.fasta*, the KMeans cutoff used as threshold

5) Alignment results are available in the *clusters_aligned* folder where *<sample_name>.abd_filtered.swarm.seeds.aln* is the output from MAFFT.


## CREDITS
- [Franck Stefani](franck.stefani@agr.gc.ca) - Ottawa RDC - Project lead
- Wendy Findlay - Ottawa RDC (BICoE - retired)
- [Mario Laterriere](mario.laterriere@agr.gc.ca) - Québec RDC
- [Jackson Eyres](jackson.eyres@agr.gc.ca) - Ottawa RDC (BICoE)

[Go to Top](#top)
