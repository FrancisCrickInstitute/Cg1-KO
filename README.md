# γδ T cell receptor dependencies define a unique immunosurveillance modality

## Nicolas Veland, Bethania Garcia-Cassani, Angela Zarco-Cuadrillero, Annamaria Mavrigiannaki, Ambra Natalini, Josephine Eum, Alejandro Suarez-Bonnet, Pierre Vantourout, Duncan R. McKenzie, Jannik Franken, Ana V. Marin, Anett Jandke, Rosa Andrés-Ejarque, Jessica Strid, Adrian Hayday and Miguel Muñoz-Ruiz.

Scripts used for the processing and analysis of the bulk RNA-sequencing and ATAC-sequencing data for the above study.

**Instructions:**

1.  Start with RNAseq directory, which contains raw counts output from nfcore as rsem.merged.gene_counts.tsv and the coldData.txt files for each tissue analyzed (e.g. Epidermis_RNAseq).

2.  For a complete nfcore downstream analysis from the beginning, start with script named 01_processing.qmd. The follow the scripts in the number sequence (i.e. 02_script, 03_script, etc). This will generate the entire analysis performed for RNAseq.

3.  Alternatively, if running nfcore is desired: the nfcore_RNAseq directory contain the script to run nfcore_rnaseq.sh and the samplesheet.csv files.

4.  For specific plots in a specific figure (e.g. RNAseq/Figure_2), there is a specific directory for each Figure containing the script as a .qmd file with the name of the panel (e.g. Figure_2g_Epidermis.qmd) and the input data to start that script (e.g. Epidermis_DESeq.csv) and reproduce the plot.

5.  The ATACseq directory contains scripts to run nfcore, scripts for DA analysis and outputs .tsv files for each tissue in their specific directories. It also contains the .qmd scripts to reproduce specific figure plots.
