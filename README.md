### "Chronic Marijuana Use Is Associated with Gene Expression Changes in BAL"
 
**Summary:** We found differences in bronchoalveolar lavage (BAL) gene expression profiles from long-term, frequent marijuana smokers without concurrent tobacco use compared to those of tobacco-only smokers and of non-smokers, including differentially expressed genes and pathways suggestive of altered respiratory immunity.

* This research letter is open access at [ATS Journals](https://www.atsjournals.org/doi/10.1165/rcmb.2021-0285LE).
* Deidentified participant metadata and raw RNA-seq data are publicly available at GEO, accession [GSE155213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155213).

### File Structure
```
├─ Data
│   └─ count_table.txt             # RNA-seq counts (gene x sample)
│   └─ participant_metadata.txt    # (sample x variable)
├─ Scripts                         # 22 scripts; data processing & analysis
└─ Results                         
│   └─ ResultsTable1.xlsx          # pairwise diff. expression results (multiple tabs) 
│   └─ ResultsTable2.xlsx          # GO enrichment results (multiple tabs)
```

### Links

* doi: [10.1165/rcmb.2021-0285LE](https://doi.org/10.1165/rcmb.2021-0285LE)
* PMID: [35103554](https://pubmed.ncbi.nlm.nih.gov/35103554/)
* Repo: [chooliu/ajrcmb_marijuana_rnaseq](https://github.com/chooliu/ajrcmb_marijuana_rnaseq) 

### Citation

> Liu C, Gaydos J, Johnson-Paben R, Kechris K, Burnham EL, Sharma S. Chronic Marijuana Use Is Associated with Gene Expression Changes in BAL. American Journal of Respiratory Cell and Molecular Biology. 2022 Feb;66(2):238-9.
