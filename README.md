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
