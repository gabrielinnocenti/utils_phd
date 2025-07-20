# 1. Introduction
In the present project we are trying to isolate from a vast variety of bacterial genomes belonging (mostly) to human microbiome, proteins which may play a role in cancer implications and immune-system modulation. Our scope, is to find out if among these variety of proteins, we can reclute specific disease markers that could be detected with appropriate methods (PhIP-Seq).
## 1.1 Selection of functionally relevant proteins
We screened the literature and decided to search for proteins which could be important for immune system modulation, biofilm formation, genotoxic activity, potential immunogeneity and providing some functions such as antibiotic resistance.
### 1.1.1 Adhesins
Adhesins are are cell surface proteins which mediate adhesion of bacteria to host cells and participate in biofilm formation. These are often involved in bacterial virulence and pathogenesis such as FimH.
The adhesins differ from Gram positive to Gram negative bacteria. In Gram negative bacteria, these are structured in fimbri or pili form (e.g. FimH) - these are typical for example of Escherichia coli or Salmonella enerica. In Gram positive bacteria, these are often surface proteins or secreted proteins (in capsule or vesicle form, such as in Porphyromonas gingivalis). Also some enzimes can act like adhesins, such as enolases.
Classes of adhesins include fimbriae, pilin, lectin, microbial surface component recognizing adhesive matrix molecules (MSCRAMMs), agglutinins, mucus-binding, collagen-binding, fibronectin-binding, fibrinogen-binding, laminin-binding, and other types of adhesins. We removed pre-pilin from the list, as these are not adhesins but rather precursors of pilin proteins.
### 1.1.2 Flagellins
Flagellins are structural proteins which are present in the flagella of many (gram-negative bacteria), which has a fundamental role in motility and chemotaxis. Flagellns are also important for immune responses in humans, as they activators of a wide range of cell types within the immune system. Specifically they are a ligand to host recognition receptors, such as the Toll like receptor 5 (TLR5), exposed on the cell surface of dendritic cells and T-lympohocytes. Also they can worg as an agonist for citosolic NAIP5/NLRC4 inflammasome. For these charateristics, flagellins have been also used as adjuvants in vaccines and immunotherapeutics. 
Their role has been studied in the context of IBD (especially Crohn's disease) and colorectal cancer, playing a key role in the interplay between immune modulation and cancer development.
Check
https://doi.org/10.1016/j.coviro.2023.101330 - Flagellins biology review
https://doi.org/10.1016/j.fsi.2013.02.029 - Flagellins biology review
https://doi.org/10.1038/s41598-019-47347-6 - Flagellin C and Colorectal cancer (research article)
doi: 10.1172/JCI20295. - Crohn disease
### 1.1.3 Ureases
Ureases are enzymes that catalyze the hydrolysis of urea into ammonia and carbon dioxide. They are metallo enzymes produced by many bacteria and plants as well. 
Check:
https://doi.org/10.1016/j.toxicon.2015.11.020 - review on Ureases as multifunctional toxin proteins
https://doi.org/10.1016/j.jare.2018.05.010 - general review on Ureases

### 1.1.4 Invasion proteins (invasins, intimins etc.)
These proteins are a subcategory of adhesins which often mediate pathogenic functions allowing bacteria to attach to a specific host cell membrane and access the cell. They are present in most pathogenic bacteria, and is one of they characterize the most common infection strategies.

### 1.1.5 Cancer associated bacterial proteins
These proteins were studied in the context of cancer development and progression, which mechanistically directly or indirectly enhanced cancer risk, through activation of certain oncogens, genotoxic activity, immune system modulation.
The specific proteins are: 
- FadA and Fap2 proteins from Fusobacterium nucleatum (subsp. animalis)
- colibactin
- Bacteroides fragilis toxin (BFT) 
- fragilysin
- PCWBR2

### 1.1.6 Antibiotic resistance proteins
[[[To add]]]

## 1.2 Search strategy
We performed a pattern search across all proteins using standard annotations provided by RefSeq files downloaded from NCBI. These annotations are based on PGAP annotation tool which is the standard procedure for storing them in the NCBI database.
### 1.2.1 Problem of hypothetical proteins
To minimize the lack of annotation details concerning some proteins, we also ran Interproscan against all proteins and merged the annotations with the previous ones (not provided in the example in this test set).

# 2. Running intructions
Run the script by using:

```bash
python3 get_relevant_proteins.py --outdir test --config config.yaml
```
Needed files:
- config.yaml
- Overview_library_10k.tsv - A specifically formatted file output from create_protein_overview.py
- taxonomic_ranking_library.csv - A specifically formatted file output from create_protein_overview.py
- Optional: Interproscan output. (not in the test examples)

# 3. Conclusions
