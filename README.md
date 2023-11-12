# NOVOLoci - Haplotype-aware assembly of long-sequencing reads

NOVOLoci is a haplotype aware assembler for targeted assembly or whole genome assembly of small genomes.
Both HiFi and ONT reads can be used seperatly or combined.

Will be available soon.

## Cite

Will soon be available on BioRxiv 

## Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/NOVOLoci/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

If your assembly was unsuccessful, you could already add the log file and the configuration file to the mail, this could help me to identify the problem!

<a href="https://gbiomed.kuleuven.be/english/cme/research/laboratories/laboratory-for-cytogenetics-and-genome-research/laboratory-for-cytogenetics-and-genome-research" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/49/KU_Leuven_logo.svg/768px-KU_Leuven_logo.svg.png" width=auto height="70"></a><a href="https://groups.oist.jp/grsu" target="_blank"><img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/ef/OIST_logo.png/1200px-OIST_logo.png" width=auto height="70"></a>


## Instructions

### Targeted approach

### 1. Find a suitable seed

- Sequence from a reference genome or from a previous assembly
- Length should be at least 500 bp
- Make sure you take a sequence before the complex region that you target, DO NOT start in repetitive or duplicated region!
- The format should be like a standard fasta file (first line: >Id_sequence)


### 2. Create configuration file

You can download the example file (config.txt) and adjust the settings to your liking.  
Every parameter of the configuration file is explained below. 

### 3. Install dependencies

- Install BLAST
- Install MAFFT


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
Project name          = /path/to/seed_file/Seed.fasta
Assembly length       = 1000000
Subsample             = 
Save assembled reads  = 
Seed Input            = /path/to/seed_file/Seed.fasta
Split haplotypes      = yes
Reference sequence    = 
Variance detection    = 
Cores                 = 30
Output path           = /path/to/output_folder/

Nanopore reads:
-----------------------
Nanopore reads        = /path/to/reads/
Local DB and NP reads =
Sequencing depth NP   = 
Use Quality scores    =

PacBio reads:
-----------------------
PacBio reads          = /path/to/reads/
Local DB and PB reads =
Sequencing depth PB   = 

</pre>

**2. Explanation parameters:**
<pre>

Project:
-----------------------
Project name         = Choose a name for your project, it will be used for the output files.
Assembly length      = If you want the assembly to terminate after a certain length, you can use this option
Subsample            = This option is currently not available
Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)
Seed Input           = The path to the file that contains the seed sequence.
Split haplotypes     = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed 
                       originates from the same sample or a close reference (yes/no)
Reference (optional) = This option is currently not available.
Variance detection   = This option is currently not available.
Cores                = It is strongly adviced to use multiple cores for the assembly, give here the available cores
Output path          = /path/to/output_folder/

Nanopore reads:
-----------------------
Nanopore reads        = /path/to/reads/
Local DB and NP reads =
Sequencing depth NP   =
Use Quality scores    =

PacBio reads:
-----------------------
PacBio reads          = /path/to/reads/
Local DB and NP reads =
Sequencing depth PB   = 

</pre>
