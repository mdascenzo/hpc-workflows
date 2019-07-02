

###### Cloning the repository:
```

```

###### Example usage:
```
snakemake --configfile config.yml -s rnaseq.smk --cores 8
```

###### Analysis Workflow Options
![workflow-full](doc/rnaseq/img/dag_options_3.png)

###### Todo:
- Possibly create installer to move R-code to library path
- Update genome and transcriptome annotation files


##### STAR

Genome Indexing:
- Configured to utilize a splice junction database.
- sjdbOverhang option should be set to ReadLength-1 (e.g. 74 for 75bp PE reads)

options:
- sjdbGTFfile
- sjdbOverhang

Example splice junction db file:
```
/Precyte1/stage/refs/genomes/hg38/annotation/hg38wERCC92/gencode.v25.primary_assembly.annotation.wERCC92.gtf
```