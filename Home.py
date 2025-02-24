import streamlit as st

# Title of the app
st.title("NuoHMMER Web v1b")
st.info('This site is actively being improved, and new features may be added. Core functionalities are stable for research use.')

# Welcome text using Markdown
st.markdown("""
### Overview

This app provides an interactive platform for researchers to explore the results of our comprehensive study on the distribution and diversity of Complex I (NADH:ubiquinone oxidoreductase) across prokaryotic genomes. Our dataset includes over 45,000 genomes, and we have performed extensive analyses to investigate the presence, variation, and phylogenetic relationships of Complex I subunits in bacteria and archaea.

### Key Steps in Our Study:

- **Genome Selection and Curation**:
  We selected **47,278 prokaryotic genomes** from the [NCBI Genome database](https://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/) (September 2024), encompassing **10,608 bacterial species** and **431 archaeal species**. The dataset focused on genomes at the "Complete Genomes" or "Chromosomes" level.

- **Extraction and Analysis of Complex I Subunits**:
  Coding sequences (CDS) were screened to identify genes encoding Complex I subunits, supplemented by additional protein sequences from the [InterPro](https://www.ebi.ac.uk/interpro/) database. This combined dataset was used to build Hidden Markov Model (HMM) profiles specific to each subunit.

- **Construction of HMM Profiles**:
  We employed [MMSeq2](https://github.com/soedinglab/MMseqs2) for sequence clustering (at a 85% identity threshold) and [MAFFT](https://mafft.cbrc.jp/alignment/software/) for alignment. HMM profiles were then constructed using [HMMER3](http://hmmer.org/) for alignment, with subunit profiles developed including bacterial and archaeal sequences.

- **Identification of Complex I Subunits**:
  Proteomes from the genomes were generated using [Prodigal](https://github.com/hyattpd/Prodigal) and searched for Complex I subunits using the developed HMM profiles.

- **Species Classification**:
  We utilized [TaxonKit](https://bioinf.shenwei.me/taxonkit/) for species classification and curated additional physiological metadata from the [BacDive](https://bacdive.dsmz.de/) database and literature.

- **Phylogenetic Analysis of Complex I Variants**:
  A phylogenetic tree of concatenated Complex I subunit sequences was generated using [IQ-Tree](http://www.iqtree.org/), employing the **Le-Gascuel 2008 model** to explore phylogenetic relationships.

### Explore Our Results

Use the navigation menu to explore our results, examine the distribution and diversity of Complex I subunits across different species.
""")

# Footer Section
st.markdown("""
---

### Contact
For inquiries, please contact:
**Akshay Shirsath**
Email: akishirsath@gmail.com

### GitHub Repository
Workflow and Pipeline Code:
[GitHub Repository Link](https://github.com/Anand-Research-Group/complex-I)
""")
