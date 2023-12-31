# Dosage-score
## Table of contents
 - [Introduction of Dosage-score](Introduction-of-Dosage-score)
 - [Installation](#Installation)
   + [Dependencies](#Dependencies)
   + [Installation using bioconda](#Installation-using-bioconda)
 - [Usage](#Usage)
   + [command : dosage-score](#command-:-dosage-score)
   + [command : dosage-score-plot](#command-:-dosage-score-plot)
 - [The input file format](#The-input-file-format)
 - [The example of execution](#The-example-of-format-execution)
 - [The output file format](#The-output-file-format)

## Introduction of DNAMarkMaker
  
  Dosage-score is a tool for estimating dosage of each genomic region from Next-Generation Sequencing (NGS) alignment data.

## Installation
### Dependencies
   - [samtools](https://github.com/samtools/samtools) >= 1.16.1
   - python >= 3
   - pandas
   - matplotlib.pyplot
   - numpy

### Installation using bioconda
  You can install DNAMarkMaker using bioconda.
  ```
  conda install -c bioconda dosage-score
  ```
  Alternatively, if you want to create DNAMarkMaker specific environment.
  ```
  conda create -n Dosage-score -c bioconda dosage-score
  conda activate Dosage-score
  ```
  
## Usage
### command : dosage_score
  Purpose: Calculate the Dosage-score of each genomic region from the NGS alignment data.

  Required options
  ```
  -r1 Fasta                     Full path of reference fasta file
  -g Txt                        Genome information file
  -b Txt                        BAM information file
  -o Directry                   Output directory
  ```

  Additional options
  ```
  -r2 Fasta                     Full path of reference fasta file with repeat region replaced by N [None]
  -l Link                       Link file generated by GetTwoGenomeSyn.pl in NGenomeSyn package [None]
  -t Int                        Thread [1]
  -minMQ Int                    Minimum mapping quality [60]
  -minBQ Int                    Minimum base quality [13]
  -window_size Int              Sliding window size [2000000]
  -step_size Int                Step size [500000]
  ```

### command : dosage_score_plot
  Purpose: Plot the Dosage-score.

  Required options
  ```
  -r1 Fasta                     Full path of reference fasta file used to "dosage-score" command
  -b Txt                        BAM information file used to "dosage-score" command
  -o Directry                   Output directory used to "dosage-score" command
  ```

  Additional options
  ```
  -chr_position Str:Int:Int     Specific chromosomal regions focused for plotting. [Whole genome]
  -window_size Int              Sliding window size used to "dosage-score" command [2000000]
  -c Txt                        Plot color file [None]
  -min_rate Float               Excloud position less than minimum ratio (=analysis/window size) [0.05]
  ```

## The input file format
### Genome information file (-g)
  This is a genome informational file sepalated tab.

  Example:
  ```
  A A01,A02,A03,A04,A05,A06,A07,A08,A09,A10
  C C1,C2,C3,C4,C5,C6,C7,C8,C9
  ```
  Colom 1: Genome name.
  Colom 2: Chromosome names listed in reference fasta separated by ",".

### BAM information file (-b)
  This is a sample informational file sepalated tab.

  Example:
  ```
  Rapa	/Full/path/Rapa.bam	Control	A
  Oleracea	/Full/path/Oleracea.bam	Control	C
  Napus	/Full/path/Napus.bam	Test	A,C
  BC1F1	/Full/path/BC1F1.bam	Test	A
  ```
  Colom 1: Sample name.
  Colom 2: Full path of BAM file.
  Colom 3: Control or Test. "Control" samples are diploid species used ancistrial filter. Positions that exhibit coverage depth =0 in selected genomes and coverage depth >0 are unselected genomes are subsequently filtered out in "Control" samples. 
  Colom 4: The genome name entered in the 1st column of the Genome information file where the dose of the sample is expected to be 2. If the dose is expected to be 2 for multiple genomes such as allotetraploid, you can enter 2 genome separated by ",". Dosage-score calculation uses median values for selected genomic regions for correction.

### Reference fasta file with repeat region replaced by N (-ref2)
  This is a fasta file with repeat region replaced by N used repeat filter. If a file is specified, N regions are excluded from the analysis.
  Example:
  ```
  >A01
  ACATGCTNNNNNNNNCGATCGTACGACT
  GATGCTAGCATACGCATCAGCAGCATGG
  .
  .
  .
  >C1
  GATCTNNNNNNNNNNNNNNNNNNNNNNN
  NNNNNNAAGCTACGATCGATACGGGGCA
  .
  .
  .
  ```
  Fasta file that combines those generated by repeat masker for each fasta file of each subgenome

### Link file (-l)
  This is a link file generated by GetTwoGenomeSyn.pl in [NGenomeSyn](https://github.com/hewm2008/NGenomeSyn) used multi- and non-homologous filter.
  If a file is specified, only the specified area will be used for analysis.

  Example:
  ```
  A01     1296792 1323561 C1      1312469 1341862
  A01     240509  280133  C1      254188  301660
  A01     32967977        32989139        C1      40100444        40122565
  .
  .
  .
  ```

### Plot color file (-c)
  This is a color informational file sepalated tab. You can specify a color for each chromosome when plotting. By default all plots are black.

  Example:
  ```
  A01 #006400
  A02 #ff0000
  A03 #0000ff
  .
  .
  .
  ```
  Colom 1: Chr name.
  Colom 2: Html color code.

## The example of execution

### dosage-score
 ```
 dosage-score -r1 /Full/path/Nupus_ref.fa \
              -r2 /Full/path/Nupus_ref_replaced_N.fa \
              -l /Full/path/A_and_C_genome.link \
              -g /Full/path/genome_info.txt \
              -b /Full/path/bam_info.txt \
              -o output
 ```                           

### dosage-score-plot
 ```
 dosage-score-plot -r1 /lustre7/home/lustre3/segawa-tenta/BnapusXBrapa/230301_alignment/make_refarence_index/Chiifu_V4_OX-heart_923_BVRC.fa \
                   -b bam.txt \
                   -o output \
                   -c chr_color.txt
 ```   

## The output file format
### Output folder by example command
 ```
 output
 |---log
 |   |---info.txt
 |
 |---1_linked_regions
 |   |---analysis_region.txt
 |
 |---2_sliding_window_average
 |   |---region_txt
 |   |   |-A01.1.2000000
 |   |   |-A01.500000.2500000
 |   |   |-A01.1000000.3000000
 |   |   ...
 |   |   
 |   |---sliging_window_total
 |   |   |-A01.1.2000000
 |   |   |-A01.500000.2500000
 |   |   |-A01.1000000.3000000
 |   |   ...
 |   |
 |   |---sliging_window_average.txt
 |
 |---3_dosage_score
 |   |---dosage_score.txt
 |
 |---4_plot
     |---filtered_dosage_score.txt
     |---Rapa.png
     |---Oleracea.png
     |---Napus.png
     |---BC1F1.png
 ```

### info.txt in log
 This file describes the information entered

### analysis_region.txt in 1_linked_regions
 This file lists the chromosome regions to be analyzed based on the link file.

### A01.1.2000000 etc. in 2_sliding_window_average/region_txt
 This file is created for each window based on the reference array. The file name is the chromosome name and the start and end positions of the window. The regions where homolog regions are defined between subgenomes by link files are listed.
 
 Colom 1: Chr name.


 Colom 2: Start position.


 Colom 3: End position.

### A01.1.2000000 etc. in 2_sliding_window_average/sliging_window_total
 This file is created for each window based on the reference array. The file name is the chromosome name and the start and end positions of the window.
 
 Colom 1: Number of positions analyzed after filtering.


 Colom 2~: Total coverage depth for each sample.

### sliging_window_average.txt in 2_sliding_window_average
 Average coverage depth list for each window. Windows with analysis position =0 are excluded by the filter.
 
 Colom 1: Chr name.


 Colom 2: Position.


 Colom 3: Number of positions analyzed after filtering.


 Colom 4: Average coverage depth for each sample.

### dosage_score.txt in 3_dosage_score
 Dosage-score list for each window.

 Colom 1: Chr name.


 Colom 2: Position.


 Colom 3: Number of positions analyzed after filtering.


 Colom 4: Dosage-score for each sample.

### filtered_dosage_score.txt in 4_plot
 A file listing only the windows where the number of analysis positions per window exceeds the min_rate value. Same format as 3_dosage_score/dosage_score.txt

### Rapa.png etc. in 4_plot
 A png file plotting the Dosage-score for each sample.
