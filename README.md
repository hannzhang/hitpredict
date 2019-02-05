# `BCB420.2019.hitPredict`

#### (hitPredict data annotation of human genes)

&nbsp;

###### [Han Zhang](https://orcid.org/0000-0002-1134-6758), University of Toronto, Canada. &lt;helenzhang0814@outlook.com&gt;
----

# 1 About this package:

This package describes the workflow of download the hitPredict data on human genes from the [hitPredict database](http://hintdb.hgc.jp/htp/index.html), how to map HGNC symbols to hitProdict domains, how to annotate the example gene set and provides examples database statistics.


&nbsp;

#### In this project ...

```text
 --BCB420.2019.hitPredict/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.Pfam.Rproj
   |__DESCRIPTION
   |__dev/
      |__rptTwee.R
      |__toBrowser.R               # display .md files in your browser
   |__inst/
      |__extdata/
         |__sym2enst.RData         # ENSP ID to HGNC symbol mapping tool
      |__img/
         |__[...]                  # image sources for .md document
      |__scripts/
         |__recoverIDs.R           # utility to use biomaRt for ID mapping
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md                    # this file

```

&nbsp;

----

## 2 hitPredict Data

Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are combined, annotated and scored. The reliability score is calculated based on the experimental details of each interaction and the sequence, structure and functional annotations of the interacting proteins.


&nbsp;


#### 2.1 Data semantics

Protein-protein interactions from IntAct, BioGRID, HPRD, MINT and DIP are combined, annotated and scored. The reliability score is calculated based on the experimental details of each interaction and the sequence, structure and functional annotations of the interacting proteins.

HitPredict is a resource of experimentally identified physical protein-protein interactions. It was compiled in the following manner:

All physical interactions were downloaded from IntAct, BioGRID, HPRD, DIP and MINT. These were combined to form a non-redundant dataset. Species with less than 10 interactions were removed.

All non-physical interactions, such as genetic interactions, were excluded.

Interactions of proteins from different species were also excluded.

All proteins were assigned valid UniProt IDs. In cases where UniProt IDs could not be directly assigned, protein sequences were aligned to UniProtKB using BLAST and the ID of the longest hit with 99% sequence identity was used.

Interactions in which one or both of the proteins were no longer present in UniProt, were removed.

Interactions of proteins that did not map to valid UniProt IDs were removed. Interactions were excluded only after manual confirmation.

The interacting proteins were annotated with Entrez gene IDs, Ensembl IDs, Pfam domains and Gene Ontology terms.

&nbsp;

## 3 Data download and cleanup

To download the source data from hitPeodict ... :

1. Navigate to the [**Ensembl**](http://useast.ensembl.org/index.html) and follow the link to the [BioMart](http://useast.ensembl.org/biomart/martview/b13049612725fd3cd2d8e840477651f0).
2. Choose "Ensembl Genes 95" as database and "Human genes (GRCh38.p2)" as dataset.
3. Select Filters -> Expand Protein Domains and Families -> Limit genes to ... select "With Pfam domain ID(s)"
4. Select Attributes -> Expand Gene section -> Select APPRIS annotation, Gene Stable ID, and Transcript Stable ID 
5. Expand External References -> Select HGNC symbols
6. Expand Protein Domains and Families -> Select Pfam domain ID, Pfam start, Pfam end
7. Click "Results" button in top left corner
8. "Export all results to ..." -> Select "Compressed file (.gz)" and "TSV". Press "Go"
9. Save downloaded file in sister directory "data"

* `mart_export.txt.gz` (1.5 Mb)	;

&nbsp;

Also download Pfam domain description via Pfam website.

1.  Navigate to [**hitPredict**](http://hintdb.hgc.jp/htp/index.html)
2.  Click on the download() link at the top, which will direct you to the download page
3.  Download H_sapiens_interactions.txt.tgz and save it in sister directory data

* `H_sapiens_interactions.txt.tgz` (13.2 MB) ;

&nbsp;

## 4 Mapping HGNC symbols to hitPredict domain IDs

#### Preparations: packages, functions, files

To begin, we need to make sure the required packages are installed:

**`readr`** provides functions to read data which are particularly suitable for
large datasets. They are much faster than the built-in read.csv() etc. But caution: these functions return "tibbles", not data frames. ([Know the difference](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html).)
```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```

&nbsp;

**`biomaRt`** biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```

&nbsp;

**`dplyr`** is a convenient package to transform and summarize tabular data with rows and columns. It contains a set of functions that perform common data manipulation operations such as filtering roles, re-ordering rows and summarizing data.;

```R
if (! requireNamespace("igraph")) {
  install.packages("igraph")
}
```

&nbsp;

Finally we fetch the HGNC reference data from GitHub. (Nb. This shows how to load `.RData` files directly from a GitHub repository!)

&nbsp;

```R
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame

```

&nbsp;

#### 4.1 Step 1: Which symbols do we have to map ?

&nbsp;

```R
# Load Pfam data with start and end coordinates
data <- as.data.frame(readr::read_delim(file.path("../downloads", "H_sapiens_interactions.txt"),
                          delim = "\t",
                          skip = 1,
                          col_names = c("uniprot1", "uniprot2", "a",
                                        "b", "Entrez1", "Entrez2", "Ensembl1", "Ensembl2", 
                                        "Taxonomy", "type", "source", "method_s", 
                                        "annotation_s", "interaction_s", "confidence")))
                      
head(data)
#   uniprot1	uniprot2	a	b	Entrez1	Entrez2	
# 1 Q16665	P40337	HIF1A_HUMAN	VHL_HUMAN	3091	7428
# 2 P04637	Q09472	P53_HUMAN	EP300_HUMAN	7157	2033
# 3 Q9Y6K9	O14920	NEMO_HUMAN	IKKB_HUMAN	8517	3551
# 4 Q5S007	Q5S007	LRRK2_HUMAN	LRRK2_HUMAN	120892	120892
# 5 P49427	P62877	UB2R1_HUMAN	RBX1_HUMAN	997	9978
# 6 Q12834	Q13257	CDC20_HUMAN	MD2L1_HUMAN	991	4085
#   Ensembl1	Ensembl2	Taxonomy	type source  method_s annotation_s interaction_s confidence
# 1 ENST00000323441 [Q16665-2],ENST00000337138 [Q16665-1],ENST00000539097 [Q16665-3]	ENST00000256474 [P40337-1],ENST00000345392 [P40337-2]	9606	Small-scale	intact
# 2 ENST00000269305 [P04637-1],ENST00000420246 [P04637-2],ENST00000445888 [P04637-1],ENST00000455263 [P04637-3],ENST00000504290 [P04637-9],ENST00000504937 [P04637-7],ENST00000510385 [P04637-8],ENST00000610292 [P04637-4],ENST00000610538 [P04637-6],ENST00000617185 [P04637-2],ENST00000619485 [P04637-4],ENST00000620739 [P04637-4],ENST00000622645 [P04637-5]	ENST00000263253	9606	Small-scale	intact	0.995550000000000046	1	0.997772519164563998	High
# 3 ENST00000594239 [Q9Y6K9-1],ENST00000611071 [Q9Y6K9-1],ENST00000611176 [Q9Y6K9-3],ENST00000618670 [Q9Y6K9-2]	ENST00000416505 [O14920-4],ENST00000519735 [O14920-3],ENST00000520810 [O14920-1],ENST00000520835 [O14920-2]	9606	Small-scale	intact	0.995099999999999985	1	0.997546991374340997	High
# 4 ENST00000298910	ENST00000298910	9606	Small-scale	intact	0.995059999999999945	1	0.997526941992044947	High
# 5 ENST00000215574	ENST00000216225	9606	Small-scale	intact	0.993829999999999991	1	0.996910226650324027	High
# 6 ENST00000310955,ENST00000372462	ENST00000296509 [Q13257-1],ENST00000333047 [Q13257-2]	9606	Small-scale	intact	0.993580000000000019	1	0.996784831345260969	High


# how many unique IDs do we have to map?
  uENSP <- unique(c(tmp$a, tmp$b))  # 19,354 IDs need to be mapped

&nbsp;

#### 4.2 Step 2: Check if each symbol only maps to one principal transcript

There are some symbols that have multiple principal transcripts (APPRIS=principal1). I have checked the protein sequences for these transcripts and they all get the same sequence. Therefore, we could just remove the additional transcript and keep one principal transcript.

```R
Here is a list of assets provided with `rpt` and why they are included. You can delete everything you don't need, but note: you can't push empty directories to your repository. Make sure you keep at least one file in every directory that you want to keep during development.
 
```
.gitignore                     <- defines files that should not be committed to the repository
.Rbuildignore                  <- defines files that should not be included in the package
DESCRIPTION                    <- the metadata file for your package
dev                            <- optional: see (Note 1)
dev/functionTemplate.R         <- optional: see (Note 1)
dev/rptTwee.R                  <- optional: see (Note 1)
inst/                          <- optional: see (Note 2)
inst/extdata/                  <- optional: see (Note 3)
inst/extdata/test-lseq.dat     <- optional: see (Note 3)
inst/scripts/                  <- optional: see (Note 4)
inst/scripts/scriptTemplate.R  <- optional: see (Note 4)
LICENSE                        <- License(s)
man/        