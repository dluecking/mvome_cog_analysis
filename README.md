# mvome_cog_analysis
## input needed:
**microbial reads (unzipped):**
- heligoland reads, which are a subsample (1M reads each pair), of the 4 heligoland filter samples (H_022, H_045, H_080, H_300), where we subsampled 250k reads from each.
- 1M subsamples of the tara oceans microbial fraction of each of the three stations


**peDNA reads (zipped):**
- heligoland peDNA reads that ran through mviest and were marked as "true.mvome" from the four samples (the four CsCl ones) combined into on sample.
- tara oceans virome fraction for the 9 stations


**sorted MAGs according to type:**
- run
```R
Rscript scripts/concatenate_MAGs_by_type.R
```

## this pipeline then runs the following steps one after another:
1. subsample reads
2. zip subsampled reads
3. map reads to contigs
4. unzip mapped reads
5. run fraggenescan on mapped reads
6. blast fragments against COG
7. run fraggenescan on microbial reads
8. blast microbial fragments against COG
9. collect the blast data
10. plot

## needed programs:
- seqtk
- bbmap
- fraggenescan
- diamond
- r with:
	- data.table
	- stringr
	- ggplot2
	- dplyr
	- ggpubr
- python with:
	- pandas

## setting up conda env:
Install with:
```bash
conda install -y -c bioconda fraggenescan diamond bbmap seqtk
conda install -y -c conda-forge r-data.table r-ggplot2 r-dplyr r-ggpubr
conda install -y -c r r-stringr
conda install -y -c anaconda pandas
```
Then export with:
```
conda env export | grep -v "^prefix: " > env/mvome_cog_analysis.yml
```
