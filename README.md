### HNSCC-HPVProject ###
This project aims to characterize key genes which are assocaited with key transcriptome changes in HNSCC HPV-Positive and HPV-negative subtypes. Also, the study aims to address the heterogeneity problem

### Project Overview ######

Head and Neck Squamous Cell Carcinoma (HNSCC) is group of malignancies that results from the head and neck squamous cell linings. It is the seventh most prevalent cancer worldwide accounting for ~4.5% of the global death toll, with 890,000 new cases and 450,000 deaths per annum with India leading the chart and accounting for around 80% of all the cases. One of the main causes of these rising cases is due to chewing of the areca nut, which is carcinogenic in nature. Human papillomavirus (HPV) infection has been seen in the HNSCC, based on the detection of HPV DNA, notably HPV-16 it can be broadly categorized into two subtypes: HPV-positive (HPV-Pos) and HPV-negative (HPV-Neg). These two categories differ significantly in their molecular profile, etiopathogenesis, and clinical outcomes and HPV-Neg subtype accounts alone for ~75% of the cases and has a significantly worse prognosis compared to the HPV-Pos subtype. Treatment options include surgery, chemotherapy, radiotherapy, and targeted agents such as EGFR monoclonal antibody (mAb) inhibitors and two PD-1 inhibitors; however, the overall response rates have been moderate. Hence, there is a need to identify novel drug targets for better health outcomes.
Most transcriptome wide studies for HNSCC HPV subtypes have focused on differentially expressed genes (DEGs) which does not necessarily reveal the genes that mediate the transcriptomic changes. In contrast, our network-based tool, PathExt mines transcriptomic data for differentially expressed paths in a curated gene network to identify central genes mediating the transcriptomic response. We have previously established the superiority of PathExt over the conventional approaches relying on DEGs. (Agrwal et.al., iScience, 2024 Apr 16;27(5):109752, Agrawal et.al., Front. Immunol. 2022 Jul1:13:918817).
Here we apply PathExt to TCGA-HNSCC on 501 TCGA-HNSCC samples (64 HPV positive & 437 HPV negative) to characterize key genes and found that our method outperformed DEGs in both subtypes in recapitulating disease etiology. PathExt genes were enriched for processes like “epithelial cell proliferation” as well as subtype-specific processes (immune and metabolic related processes in subtype 1 and peptide-related processes in subtype 2). Additionally, PathExt key genes show significant overlap with HNSCC-specific external validation datasets compared to DEGs. Unlike DEGs, PathExt genes exhibit significant expression in various cell types, enrichment for cancer hallmarks and mutated protein systems. Support Vector Machine (SVM) model developed using PathExt key gene’s expression outperformed DEGs in classifying responders with AUROC of 0.74. Lastly, top 10 potential therapeutic targets and drugs were proposed.


**##################### Instructions for the Users ###################**

In this study, we implemented netwrok based approach to identify differentially regulated paths (also referred as TopNets). These TopNets are of two types; (i) Activated; and (ii) Repressed. Further, from these TopNets, we identify top central genes.

PathExt requires node weights as an input for each gene in the corresponding sample. For more details, refer PathExt paper at "https://academic.oup.com/bioinformatics/article/37/9/1254/5952670?login=true"

**#################### Node Weight Computation ######################**

For computing node weights, user need to provide the input data in the csv file format. The input data should consist of 3 columns. First column should be gene list. Second column should be the value of those genes in control sample and Third column should be the value of those genes in the case. Order of the column is very important for the code.

It's recommended that user should provide the gene values in the quantile normalize form.

Next, run code **"Activated_node_weight.r"** to compute node weight for the Activated Network. The codes are provided in the **"code"** folder

The provided code compute value for one sample. User can run the code in loop for multiple samples.

**Code Usage:**

**/usr/local/bin/Rscript Activated_node_weight.r**      ##### For computing Activated Node Weight<br>

**########################### Running PathExt for Generating the TopNets ##########**

1. Compute the percentile threshold and q-score at which you user want minimum nubmer of nodes in the topnet. Run the following commands inside the folder where all the python codes are present.

**a. mkdir test_data/results/temp**<br>
**b. python3 node_weight_matrix_colname_Pijs.py test_data/input_data Sample1 test_data/human_PPIN.txt 0.1 2 1000 test_data/results/Activated_response test_data/results/temp/Pij**<br>
**c. python3 fdr_rand_pijs_boxcox.py test_data/results/temp test_data/results/Pij_zscores.txt**<br>
**d. rm -rf test_data/results/temp**<br>
**e. python3 try_different_thresholds_node_weight_matrix.py test_data/input_data Sample1 test_data/human_PPIN.txt 2 test_data/results/Pij_zscores.txt test_data/results/thresh_TopNet_sizes.txt**

2. After running the following commands, select the best values from the output file **"thresh_TopNet_sizes.txt"**. For example, user selected 0.01 as percentile and 0.05 as q-score. Now run the following commands for generating the topnets with user selected pecentile and q-score.

**mkdir test_data/results/temp**<br>

**python node_weight_matrix_colname_Pijs.py test_data/input_data Sample1 test_data/human_PPIN.txt 0.01 2 1000 test_data/results/Activated_response test_data/results/temp/Pij**<br>

**python fdr_rand_pijs_boxcox.py test_data/results/temp test_data/results/Pij_zscores.txt**<br>

**rm -rf test_data/results/temp**<br>

**python benjamini_hochberg_boxcox.py test_data/results/Pij_zscores.txt 0.05 test_data/results/Pij_zscores_fdr.txt**<br>

**python extract_fdr_network.py test_data/results/Activated_response test_data/results/Pij_zscores_fdr.txt 0.05 test_data/results/Activated_Response_TopNet.txt**

Here,

a) "input_data" is a tab seprated microarray data file with all the samples to be studied;<br>
b) "human_PPIN.txt" is the unweighted network file<br>
c) "Sample1" is the name of perturbation sample to study<br>
d) "0.01" is the percentile threshold<br>
e) "2" is the path length threshold<br>
f) "0.05" is the q-score cutoff<br>
g) "1000" is the number of randomizations<br>
h) "results" is the output directory<br>
i) "thresh_TopNet_sizes.txt" is the output file with all the percentile and q-score threshold.<br>
j) "Activated_Response" is the file name for base response network (we'll put it in the output directory)<br>
k) "Activated_Response_TopNet.txt" is the file name for TopNet (we'll put it in the output directory)<br>

**########################## Computing centrality score and top central genes ###########**

After generating the topnet file, compute the centrality score of each gene by running the following command

**python calc_ripple_centrality.py test_data/results/Activated_Response_TopNet.txt test_data/results/Activated_epicenter**

Next step is sorting the output file on the basis of "ripple_centrality" score and selecting the top required genes.<br>
We have provided the example ouput file in the **"test_data/results/"** folder for the reader.

<br>

**########################## Computing DEGs ###########**

DEGs in the current study were computed using DESeq2 by running the following command

**/usr/local/bin/Rscript DESeq2.r**
The input file comprises of raw reads for case and control and metadata file. These files are provided here by the name DESeq2_input and DESeq2_metadata. The output will comprise the output generated by the DESeq2 tool. 
The code and files are present in the code folder.

User should install DESeq2 library before running the command.

<br>

**############################ Creating Venn Diagram #######################################**
Create a comma separated file to generate venn diagram. Sample input file is provided "venn_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript venn_diag.r**


**############################ Creating Boxplot #######################################**
Create a comma separated file to generate boxplot. Sample input file is provided "boxplot_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript boxplot.r**


**############################ Creating Barplot #######################################**
Create a comma separated file to generate barplot. Sample input file is provided "barplot_input.csv"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript barplot.r**

**############################ Creating Heatmap #######################################**

Create a comma separated file to generate heatmap. Sample input file is provided "heatmap_input.csv" and the generated output is a jpg image "Heatmap.jpg"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript heatmap.r**

**############################ Creating DotPlot #######################################**

Create a comma separated file to generate Dotplot. Sample input file is provided "dotplot_input.csv" and the generated output is a jpg image "Dotplot.jpg"<br>
Code is present in the **"code"** folder

**Code usage:**

**/usr/local/bin/Rscript dotplot.r**

**############################ Creating Enrichment Plot among GO terms #######################################**

Clusterprofiler package was used for creating enrichment plots among GO terms. Code used is provided in the folder **"code"**
Run the command as

**/usr/local/bin/Rscript BP_clusterprofiler.r**<br>
Note that above code can be used to identify molecular functions. User need to change **ont="BP"** to **ont="MF"** <br>

**#################### Predicting Non-responder using SVM based model ######################**

SVM model was trained on TCGA-HNSCC dataset to predict non-responder to a cetuximab treatment using top100 gene signature expression as a feature. User can use this model to predict whether a patient will respond to the treatment or not (positive label classify as non-responder). Input file is provided by the name "ml_test.csv". Model is provided by the name **HPV_finalized_model.sav**. All the 3 files (code, test file and model is present in the maiin directory)

**Code usage:**
**python SVM_predict.py**



**Tools and packages used:**<br>
Python 3.6.9<br>
R 3.6<br>
Pandas 0.25.3<br>
Networkx 1.11<br>
Numpy 1.17.4<br>
Statsmodel<br>
Random<br>
Sys<br>
Math<br>
ggplot2<br>
dplyr<br>
tidyverse<br>
pheatmap<br>
loess<br>
Clusterprofiler<br>
