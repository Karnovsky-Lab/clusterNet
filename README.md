# Identifying Metabolic Modules within Biological networks via Consensus Clustering

<font size = "4"> Gayatri Iyer, Marci Brandenburg, Christopher Patsalis, George Michailidis, Alla Karnovsky1 </font>


<h2 align="center" id="heading"><font size = 5>Special focus on clustering output from the Correlation Calculator java tool</font></h2>


### *NOTE: The output for the example code here can be seen by accessing the vignette within the package. To do so, run browseVignettes("clusterNet") in the console after installing the package.* 

<h3 align = "center" id = "heading"><font size = 4>Abstract</font></h3>

The rise of -omics data generation and the ability to extract meaningful biological insight from subsequent analyses has shone a light on the potential of metabolomics research. However, one such common analyses - enrichment analysis - is limited within metabolomics due to poorly annotated pathway databases. In addition, due to the large number of features and relatively small sample size in -omics studies *(colloquially known as the **n >>> p problem**)* many correlation-based statistical methods are not applicable to this datatype. To circumvent these issues, our group has developed several data-driven network analysis tools that take advantage of regularization techniques to create partial correlation networks on metabolomics and lipidomics data:

- Correlation Calculator for data with a single network <span style="color: red;">(Basu et al. 2017)</span>

- Filigree for data with two networks <span style="color: red;">(Iyer et al. 2020)</span>

As your dataset increases in size, the networks may become quite dense, complicating the task of identifying relavent findings. With a dense network, it is often beneficial to cluster the features within the network to create metabolic modules. Clustering will help identify important feature-feature connections within your network and highlight relavant biological associations.

<h3 align = "center" id = "heading"><font size = 4>Workflow</font></h3>
There are several steps upstream that create the edge list used in this walk-through. For a full description of these methods please see the corresponding paper for each tool, or visit [metscape](http://metscape.ncibi.org/calculator.html) to download the tools and access the user manuals. Finally, a full description of the example data and analysis can be found in <span style="color: red;">Goutman et al. 2023</span>.  

<h4 align = "left" id = "heading"><font size = 3>Required clusterNet Parameters</font></h3>

There are many approaches to clustering data that is suitable for identifying subnetworks. Our preferred method is described in <span style="color: red;">Ma et al. 2016</span> and implemented here in the clusterNet package, as well as the Filigree java tool. We focus on metabolomics and lipidomics data in our discussion, particularly for the output of Correlation Calculator. The method provided in clusterNet has been simplified into one easy-to-use function. It has several required parameters that change depending on the input type. 

<font size = "2">**FROM AN ADJACENCY MATRIX:**</font>

-*data* an adjacency matrix with feature names as column and row headers in the same order (can be a logical matrix if unweighted or numeric if weighted)

-*type*: a string corresponding to the input type. In the case of an adjacency matrix, "adjacency_matrix"

-*weighted*: a logical indicating whether the graph should be weighted or not


<font size = "2">**FROM AN EDGE LIST:**</font>

-*data*: an edge list that contains two columns of feature names, each row indicating a connection between said features. NOTE: additional columns, as is the case with the Correlation Calculator output, will not affect the algorithm

-*type*: a string corresponding to the input type. In the case of an edge list, "edge_list"

-*metaboliteA*: A string indicating the column name of the first column in the edge pair. For the example data and any output from Correlation Calculator, is is "metab1"

-*metaboliteB*: A string indicating the column name of the first column in the edge pair. For the example data and any output from Correlation Calculator, is is "metab2"


<h4 align = "left" id = "heading"><font size = 3>Running clusterNet</font></h3>

Downloading clusterNet from its [Github Repository](https://github.com/Karnovsky-Lab/clusterNet) is quick and easy. The following commands download the devtools package if not already installed, and then download the clusterNet package from GitHub.
```{r setup, eval = FALSE}
#if devtools is not already installed
install.packages("devtools")

#Download clusterNet from GitHub
devtools::install_github("Karnovsky-Lab/clusterNet", build_vignettes = TRUE)
```
Now that the package is installed, we can load it into our environment and get started! It is a known trait of clustering algorithms that limiting the input to meaningful information often improves performance. The output of Correlation Calculator will provide values every possible network connection, or (P*(P-1)/2) connections. Therefore, we need to filter our data before consensus clustering. One way to filter is to limit the edge list to only significant connections, as defined as an adjusted p-value less than or equal to 0.1. First we can access the example ALS data and check the number of connections.
```{r}
library(clusterNet)
data(ALS)
ALS[c(1,2,4,5),1:5]
dim(ALS)
```
The ALS dataset has 640 metabolites, therefore the Correlation Calculator output has (640*(640-1)/2) connections. Before running clusterNet, we need to filter this dataset for only the significant connections.
```{r}
ALS <- ALS[ALS$adj.pval <= 0.1,]
dim(ALS)
```
Filtering  reduced our dataset down to 887 connections. Now we're able to run clusterNet. The required parameters include ALS as the 'data' parameter, and "edge_list" as the 'type' parameter since we're clustering from an edge list. There is one additional parameter required regardless of the 'type', main.seed, that allows you to pick the seed for random number generation. It is set by default and does not need to be changed. Since our input is an edge list, we also need the parameters: 'metaboliteA' which equals the name of the first column or "metab1", and 'metaboliteB' which equals the name of the second column or "metab2".
```{r, results = 'hide', message = FALSE, warning = FALSE, error = FALSE}
ALScluster <- clusterNet(data = ALS, type = "edge_list", main.seed = 417,
                         metaboliteA = "metab1", metaboliteB = "metab2")
```

Similarly, we can use the same function to cluster network features with an adjacency matrix. We will provide an adjacency matrix as the 'data' parameter, "adjacency_matrix" as the 'type' parameter, and a logical value for whether the resulting graph should be weighted. We can create a simulated adjacency matrix and then run clusterNet.
```{r, eval = FALSE}
adjacencyMatrix <- matrix(data = 0, nrow = 25, ncol = 25)
rownames(adjacencyMatrix) <- paste0("metabolite", 1:25)
colnames(adjacencyMatrix) <- paste0("metabolite", 1:25)
adjacencyMatrix[runif(200, min = 1, max = 625)] <- runif(200, min = 0, max = 1)
AdjacencyCluster <- clusterNet(data = adjacencyMatrix, type = "adjacency_matrix", main.seed = 417, weighted = TRUE)
```
Now we can access our subnetwork classifications as well as a summary of the results
```{r}
head(ALScluster$subnetworks)
ALScluster$summary
```
<h3 align = "center" id = "heading"><font size = 4>References</font></h3>

-Basu S, Duren W, Evans CR, Burant C, Michailidis G, Karnovsky A.
Sparse network modeling and Metscape-based visualization methods for the analysis of large-scale metabolomics data. Bioinformatics. 2017 Jan 30 [PMC free article](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5860222/) [PubMed](https://pubmed.ncbi.nlm.nih.gov/28137712/) [Google Scholar](https://scholar.google.com/scholar_lookup?journal=Bioinformatics&title=Sparse+network+modeling+and+metscape-based+visualization+methods+for+the+analysis+of+large-scale+metabolomics+data&author=S.+Basu&author=W.+Duren&author=C.R.+Evans&author=C.F.+Burant&author=G.+Michailidis&volume=33&publication_year=2017&pages=1545-1553&pmid=28137712&doi=10.1093/bioinformatics/btx012&)

-Goutman SA, Boss J, Iyer G, Habra H, Savelieff MG, Karnovsky A, Mukherjee B, Feldman EL. Body mass index associates with amyotrophic lateral sclerosis survival and metabolomic profiles. Muscle Nerve. 2023 Mar;67(3):208-216. doi: 10.1002/mus.27744. Epub 2022 Nov 18. PMID: 36321729; PMCID: PMC9957813.[PubMed](https://pubmed.ncbi.nlm.nih.gov/36321729/) [Google Scholar](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C23&q=Goutman+SA%2C+Boss+J%2C+Iyer+G%2C+Habra+H%2C+Savelieff+MG%2C+Karnovsky+A%2C+Mukherjee+B%2C+Feldman+EL.+Body+mass+index+associates+with+amyotrophic+lateral+sclerosis+survival+and+metabolomic+profiles.+Muscle+Nerve.+2023+Mar%3B67%283%29%3A208-216.+doi%3A+10.1002%2Fmus.27744.+Epub+2022+Nov+18.+PMID%3A+36321729%3B+PMCID%3A+PMC9957813.&btnG=)

-Iyer GR, Wigginton J, Duren W, LaBarre JL, Brandenburg M, Burant C, Michailidis G, Karnovsky A. Application of Differential Network Enrichment Analysis for Deciphering Metabolic Alterations. Metabolites. 2020 Nov 24;10(12):479. doi: 10.3390/metabo10120479. PMID: 33255384; PMCID: PMC7761243.[PubMed](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7761243/) [Google Scholar](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C23&q=Application+of+Differential+Network+Enrichment+Analysis+for+Deciphering+Metabolic+Alterations&btnG=)

-Ma J., et al. (2016) Network-based pathway enrichment analysis with incomplete network information. Bioinformatics, 32, 3165–3174. [PubMed](https://pubmed.ncbi.nlm.nih.gov/27357170/) [Google Scholar](https://scholar.google.com/scholar?hl=en&as_sdt=0%2C23&q=%282016%29+Network-based+pathway+enrichment+analysis+with+incomplete+network+information.+Bioinformatics%2C+32%2C+3165%E2%80%933174.&btnG=)
