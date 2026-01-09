# ***Phacelia hastata*** (silverleaf scorpionweed)

<p align="center" width="100%">
    <img width="70%" src="images/phha_photo5.jpg"> 
</p>

## Ponderosa directory organization

```mermaid
flowchart TD;
    A(**PROJECT DIRECTORY** <br> on ponderosa: <br> /working/romero/PHHA)
    A --> B(**fastq** <br> *individual fastq files*)
    A --> C(**assembly** <br> *reference assemblies and associated files*)
    A --> D(**mapping** <br> *individual bam files, population lists, etc*)
    A --> E(**filtering** <br> *vcf files, filtering logs, hard-GT matrices*)
    A --> F("**entropy** <br> *entropy related input and output (gprob, q, dic files)*")
    C --> G(**final** <br> *chosen assembly/index files*)
    A --> H(**angsd**)
    E --> I(**final** <br> *post-filtering VCF, individual, 012 file used for entropy*)
    A --> J(**scripts** <br> *mostly bash scripts referenced in workflows*)
    H --> K(**diversity**)
    H --> L(**relatedness**)
```

## Git organization

`workflows` will have markdowns related to processses on ponderosa. Things like:
    
+ fastq --> GT calls (assembly, mapping, variant calling, filtering)
+ running entropy
+ running angsd
+ anything that generates intermediate files or output files needed for statistical analyses

All other analyses can be organized into associated directory (e.g. maps, population structure, ancestry, diversity, relatedness, GEA, phylogeny, etc.). These would hold things like:

+ input statistical files (e.g. genotype probability matrices, diversity estimates by population, etc.)
+ Rmarkdown/Python scripts for stats or figure making
+ Output figures - but be conscious about file sizes if they're exploratory. Total repo size shouldn't exceed ~ 1 GB

## Latex organization

Will add to later depending on how we want to organize. We should have seemless integration between Github, Overleaf, and Zotero (or whatever citation manager we prefer).
