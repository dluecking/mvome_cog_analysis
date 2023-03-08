# mvome_cog_analysis
## input needed:
**microbial reads (zipped):**
UNZIPPED!
- heligoland reads, which are a subsample (1M reads each pair), of the 4 heligoland filter samples (H_022, H_045, H_080, H_300), where we subsampled 250k reads from each.
- 1M subsamples of the tara oceans microbial fraction of each of the three stations



**peDNA reads (unzipped):**
- heligoland peDNA reads that ran through mviest and were marked as "true.mvome" from the four samples (the four CsCl ones) combined into on sample.
- tara oceans virome fraction for the 9 stations


## runs the following steps one after another:
1. subsample reads
2. zip subsampled reads
3. collect MAGs by type and combine
4. map reads to contigs
5. unzip mapped reads
6. run fraggenescan on mapped reads
7. blast fragments against COG
8. run fraggenescan on microbial reads
9. blast microbial fragments against COG
10. collect the blast data
11. plot

## needed programs:
- seqtk
- bbmap
- fraggenescan
- diamond
- r with:
	- library(data.table)
	- library(stringr)
	- library(seqinr)
	- library(tidyr)
	- library(ggplot2)
	- library(dplyr)
	- library(ggpubr)


## setting up conda env:
```bash
conda install -y -c bioconda fraggenescan diamond bbmap seqtk
conda install -y -c conda-forge r-data.table r-ggplot2 r-dplyr r-ggpubr
conda install -y -c r r-stringr
```

# exporting conda env to one env file
```
conda env export | grep -v "^prefix: " > env/mvome_cog_analysis.yml
```

