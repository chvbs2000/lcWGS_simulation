# lcWGS_simulation
A tool to simulate lcWGS MGMMs. This repository contains a Python script `simulate_1x.py` designed to process high-coverage whole-genome sequencing VCF files. The script takes an input VCF file and generates outputs in a directory of interest.

## Prerequisites
Before you run the script, ensure you have Python 3.8.3 installed on your system. You can download Python 3 from the official Python website.


## Installation
Clone this repository to your local machine:
```
git clone https://github.com/chvbs2000/lcWGS_simulation.git
cd lcWGS_simulation
```

## Usage
To use the script, you need to provide it with a VCF file and specify an output directory where the results will be stored. Here is the basic command structure:
```
python3 simulate_1x.py --input_vcf $vcf --output_dir $myoutdir
```

## Parameters
- `--input_vcf`: The path to the input VCF file. This should be a high-coverage whole-genome sequencing VCF file, for example, NA12878.chr22.30x.vcf.gz.
- `--output_dir`: The directory where the output files will be saved. If the directory does not exist, the script will create it.

## Example Command
Here is an example command using a specific VCF file and output directory:
```
python3 simulate_1x.py --input_vcf NA12878.chr22.30x.vcf.gz --output_dir ./output
```
This command processes the VCF file NA12878.chr22.30x.vcf.gz and places the output in the ./output directory.

### Example case: simulate chromosome 22 for genome build GRCh38
- Download NA12878 chr22 30x: The New York Genome Center and the 1000 Genomes Project have released 30x coverage data for NA12878. For downloading this data, please refer the official documentation at the following URL:
  `https://www.internationalgenome.org/data-portal/data-collection/30x-grch38`

- Download GRCh38 reference genome from below ftp site, please refer this blog for more details about reference genome: `https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use`
`ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`

- Downsampling 30x to 1x
```
HG38_FASTA="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
cram=""
samtools view -T ${HG38_FASTA} -s 1.0333 -bo NA12878.final.chr22.1x.bam $myoutdir/$sample/cram/$cram;
samtools index NA12878.final.chr22.1x.bam
```

## Output
The script will generate several files in the specified output directory, depending on the processing performed within the script. Ensure you have enough disk space and the necessary permissions to write to the directory.

## Support
For support, create an issue in this repository or contact the repository owner.
