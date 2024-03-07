# LSEA v0.2.4

*This is a new version of the tool that is now being extensively benchmarked. Please note that this is a development repository. Please use the contents of this repository with caution*

## Prerequisites
- plink v1.9
- bedtools
- scipy for Python3
- rpy2 for Python3

## Usage

This version of LSEA substitutes the p-value cutoff model with the new approach to select optimal cutoff for your data. LSEA will create result files for multiplt cutoffs and write statistics file to select the optimal parameter value. 

To run full pipeline:  
```
python3 LSEA_2.3.py -input ./data/in.tsv \
                    -plink_dir <plink_folder_path> \
                    -bfile <.bed file for PLINK> \
                    -universe ./data/universe.json
```

All results will be written to the `lsea_result` directory which will be created for you. If you want to change the location of the result files, please specify a relative path using the `-out` option. By defaul, LSEA will perform enrichment analysis using two p-value cutoffs (`1e-05, 5e-08`). If you want to test more cutoffs, please specify them using the `-p` option sepeated by space.

To generate universe from BED and GMT files, run:  
```
python3 universe_generator.py -variants ./data/positions.tsv \
                              -features ./data/gencode_formatted.bed \
                                        ./data/c2.all.v7.0.symbols.gmt \
                              -interval 100000 -o universe.json
```

To generate universe using a directory with per-set BED files, run:
```
python3 universe_generator.py -variants ./data/positions.tsv \
                              -feature_files_dir <dir with BED files> \
                              -interval 100000 -o universe_dir.json
