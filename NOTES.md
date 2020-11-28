# Notes
- Install develop prerequisites:
``` bash
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n yapc pandas matplotlib matplotlib-venn snakemake-minimal pybigwig idr wget
```
- clone [`develop` branch](https://nvie.com/posts/a-successful-git-branching-model/):
``` bash
$ git clone --branch develop https://github.com/jurgjn/yapc.git
```
- install in editable mode using pip:
``` bash
pip install -e yapc
```
- Broad overview of the IDR workflow -- replicates, pseudoreplicates, self-pseudoreplicates, etc -- in [(Landt et al. 2012)](https://dx.doi.org/10.1101%2Fgr.136184.111) around [Fig 7](http://genome.cshlp.org/content/22/9/1813/F7.expansion.html).  
- There are two IDR implementations: [an initial R version](https://sites.google.com/site/anshulkundaje/projects/idr), and a [v2 python rewrite](https://github.com/nboley/idr/) ([available in bioconda](https://anaconda.org/bioconda/idr)). The latter [has several improvements](https://groups.google.com/forum/#!topic/idr-discuss/A7PaMnzoFwg) and is being used by the current version of this pipeline.
- `idr` script in IDR v2 is the equivalent of `batch-consistency-analysis.r` in the initial R version; [it should be used in a similar fashion](https://groups.google.com/forum/#!topic/idr-discuss/_a_GKfw7kwM) (replicates, pseudoreplicates, self-pseudoreplicates)
- An explanation of [globalIDR vs localIDR](https://groups.google.com/forum/#!topic/idr-discuss/FY2K5VKx8AQ)
    > Use the global IDR for thresholding. Essentially the local IDR is akin to the posterior prob. of a peak belonging to the irreproducible noise component. The global IDR is analogous to a multiple hypothesis correction on a p-value to compute an FDR.
- Replicate-specific peaks/peak scores can be replaced with [a metric calculated from a replicate-specific coverage track at regions defined by the oracle peaks](https://groups.google.com/forum/#!topic/idr-discuss/HuS9RyLQsi4)
- Interpreting IDR v2 [diagnostic plots](https://groups.google.com/forum/#!topic/idr-discuss/KCUPDALSBmI)
- Weak but highly correlated peaks [can be problematic with IDR](https://sites.google.com/site/anshulkundaje/projects/idr/deprecated#TOC-CALL-PEAKS-ON-INDIVIDUAL-REPLICATES)
    > truncate the number of peaks to the top 100k-125k. Using more than this simply increases the running time of the IDR analysis with no advantage. Infact using more peaks with MACS2 can cause problems with the IDR model because MACS2 seems to produce strange highly correlated peak scores for very weak and noisy detections. This can confuse the IDR model.
- The method is sensitive to the choice of the smoothing window. It could, in theory, be further improved by utilising multiple windows (and wavelets), as has been previously done for mass spectrometry data in ([Du et al 2006](https://doi.org/10.1093/bioinformatics/btl355)).
