# mvome_cog_analysis

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
seqtk
bbmap
fraggenescan
diamond
r with:
	library(data.table)
	library(stringr)
	library(seqinr)
	library(tidyr)
	library(ggplot2)
	library(dplyr)
	library(ggpubr)


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

