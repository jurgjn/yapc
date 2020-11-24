# Changelog
All notable changes to [yapc](https://twitter.com/aylwyn_scally/status/1311543573435822080) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.dev0] - Unreleased
### Planned
- A comparison to a published set of mammalian ATAC-seq peak calls.

### Added
- Option to specify input data in a textual table (in addition to the current way of using command line arguments).
    - A single positional input argument will be treated as a file name pointing to a sample sheet.
    - Multiple positional arguments will be treated as previously (output prefix, conditions, and BigWig tracks).
- `--chroms` to run on a manually specified a comma-separated subset of chromosomes, e.g. for quick testing.
- Rename `--min-concave-region-width` to `--min-peaklet-width`.
- Overall tweaks to output messages and documentation.
- [A Changelog](https://keepachangelog.com/en/1.0.0/).

### Removed
- `--pseudoreplicates` to automatically choose between replicates and pseudoreplicates.

## [0.1] - 2020-02-02
### Added
- `yapc` script with the peak calling method described in [Identification of accessible sites](https://elifesciences.org/articles/37344#s3-9) ([Jänes et al., 2018](https://doi.org/10.7554/eLife.37344)).
- [Package](https://anaconda.org/bioconda/yapc)/[recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/yapc/meta.yaml)
on [bioconda](https://bioconda.github.io) ([Grüning et al., 2018](https://doi.org/10.1038/s41592-018-0046-7))
for convenient installation: `conda install -c bioconda yapc`.
