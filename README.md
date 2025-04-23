# NOVOLoci - Haplotype-aware assembly of long-sequencing reads

NOVOLoci is a haplotype aware assembler for targeted assembly or whole genome assembly of small genomes.
Both HiFi and ONT reads can be used seperatly or combined.

Currently it is only available for Nanopore, PacBio and hybrid options will be available soon.

## Cite

Will soon be available on BioRxiv 

## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/NOVOLoci/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="https://gbiomed.kuleuven.be/english/cme/research/laboratories/laboratory-for-cytogenetics-and-genome-research/laboratory-for-cytogenetics-and-genome-research" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/49/KU_Leuven_logo.svg/768px-KU_Leuven_logo.svg.png" width=auto height="70"></a><a href="https://groups.oist.jp/grsu" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/ef/OIST_logo.png/1200px-OIST_logo.png" width=auto height="70"></a>


## Instructions

### With Docker:

<code>git clone https://github.com/ndierckx/NOVOLoci.git</code>

<code>cd NOVOLoci</code>

<code>docker pull ndierckx/novoloci:latest</code>

<code>docker run --rm ndierckx/novoloci perl NOVOLoci0.1.pl -c config.txt</code>

### OR

### With Singularity:

<code>git clone https://github.com/ndierckx/NOVOLoci.git</code>

<code>cd NOVOLoci</code>

<code>singularity build novoloci.sif singularity.def</code>

<code>singularity run novoloci.sif -c config.txt</code>

### OR

### 1. Install dependencies

- Install BLAST
- Install MAFFT
- Install Perl modules: MCE::Child, MCE::Channel and Parallel::ForkManager
  
  <code>cpan install MCE</code>
  
  <code>cpan install Parallel::ForkManager</code>
  
### 2. Find a suitable seed (only Targeted approach)

- Sequence from a reference genome or from a previous assembly
- Length should be at least 500 bp
- Make sure you take a sequence before the complex region that you target, DO NOT start in repetitive or duplicated region!
- The format should be like a standard fasta file (first line: >Id_sequence)

### 3. Create configuration file

You can download the example file (config.txt) and adjust the settings to your liking.  
Every parameter of the configuration file is explained below. 


### 4. Run NOVOLoci

<code>perl NOVOLoci1.0.pl -c config.txt</code>


----------------------------------------------------------------------------------------------------------

## Configuration file

This is an example of a configuration file for the assembly of a chloroplast.
To make the assembler work, your configuration file has to have the exact same structure.
(Make sure there is always a space after the equals sign and every parameter is captured in one single line)

**1. Example of configuration file:**

<pre>

Project:
-----------------------
Project name          = projectname
Assembly length       = 1000000
Save assembled reads  = 
Seed Input            = /path/to/seed_file/Seed.fasta
Genome size           =
Ploidy                = 2
Circular              =
Threads               = 30
Output path           = /path/to/output_folder/
TMP path              = /path/to/temporary_folder/

Nanopore reads:
-----------------------
Nanopore reads        = /path/to/reads/
Local DB and NP reads = /path/to/database/
Sequencing depth NP   = 
R10                   =
Min read length NP    =

PacBio reads:
-----------------------
PacBio reads          = /path/to/reads/
Local DB and PB reads = /path/to/database/
Sequencing depth PB   = 
Min read length PB    =

</pre>

**2. Explanation parameters:**
<pre>

Project:
-----------------------
Project name          = Choose a name for your project, it will be used for the output files.
Assembly length       = If you want the assembly to terminate after a certain length, you can give the desired length; 
                        If you want to assemble the complete dataset write: "WG"
Save assembled reads  = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input            = The path to the file that contains the seed sequence.
Genome size           = Either you give the genome size (in Gbp) or you give the sequencing depth below.
Ploidy                = Give the ploidy of the sample. If it is a very heterozygous diploid species (>2%), you can give ploidy 1
Circular              = "Yes" for when the targeted sequence is circular, make sure to give an assembly length, it will try to circularize after reaching that length
Threads               = It is strongly adviced to use multiple cores for the assembly, give here the available cores
Output path           = /path/to/output_folder/
TMP path              = /path/to/temporary_folder/

Nanopore reads:
-----------------------
Nanopore reads        = Only use this when you run the dataset for the first time. 
Local DB and NP reads = If you ran the dataset before, you can give the path of the previous output folder to reuse the database
Sequencing depth NP   = Give an estimation of the sequencing depth
R10                   = If you are using R10 data, please write "yes" here
Min read length NP    = Give the minimum read length to be used in the assembly, (default: 1000)

PacBio reads:
-----------------------
PacBio reads          = Only use this when you run the dataset for the first time. 
Local DB and NP reads = If you ran the dataset before, you can give the path of the previous output folder to reuse the database
Sequencing depth PB   = 
Min read length PB    = Give the minimum read length to be used in the assembly, (default: 500)

</pre>
