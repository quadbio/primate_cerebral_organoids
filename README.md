# Script archive: Organoid single-cell genomic atlas uncovers human-specific features of brain development
Reusable scripts and functions archive for the paper "Organoid single-cell genomic atlas uncovers human-specific features of brain development" (https://www.nature.com/articles/s41586-019-1654-9). The preprint version is also available in biorxiv ["Single-cell genomic atlas of great ape cerebral organoids uncovers human-specific features of brain development"](https://www.biorxiv.org/content/10.1101/685057v1)

## File description:
 * ```extdata```: the BrainSpan RNA-seq data
 * ```scATAC_analysis```: the scripts for preprocessing and analyzing the scATAC-seq data
   * ```scATAC_1_processing_singlecells.sh```: the BASH commands to process the scATAC-seq data, including mapping, peak calling, liftover between species, etc.
   * ```scATAC_2_merge_singleCells.sh```: the BASH commands to merge scATAC-seq of all cells, and do peak calling based on the merged data
   * ```scATAC_3_merge_species.sh```: the BASH commands to create a single count matrix for scATAC-seq data from species mapped to different genomes
   * ```scATAC_analysis.R```: script to analyze scATAC-seq data in R
 * ```calculateRSS.r```: functions related to RSS (Reference Similarity Spectrum) calculation, including:
   * ```retrieveABARef```: retrieve the BrainSpan reference. It requires the extdata folder to be available
   * ```calculateRSS```: calculate RSS matrix given the input expression and reference expression matrices
   * ```plotFeature```: the plot function (based on R-base) used in the paper for feature plots
 * ```differential_detection_rate.r```: functions related to differential detection rate analysis, including:
   * ```sample_cells_ctrl_hetero```: subsample cells/nuclei by controlling subtype heterogeneity and cell number differences among groups to compare
   * ```differential_detection_rate_test```: the GLM-based ANCOVA with binomial noise to do differential detection rate test
 * ```differential_expression.r```: functions related to differential expression analysis, including:
   * ```sampling_ctrl_lines_pt```: get cells indices by random sampling with line number restricted and cell distribution along pseudotime controlled
   * ```DE_ftest_pt```: function to do DE along the aligned pseudotime trajectory
 * ```pt_alignment.r```: functions related to pseudotime alignment, including:
   * ```align_diff```: dynamic-time-warping-based alignment
   * ```get_truncate_points```: determine start/end of reference trajectory for truncated pseudotime alignment
   * ```align_pt_traj```: the wrapper function to align given expression matrix and pseudotimes of cells in reference trajectory and query trajectories
   * ```trunc_align_pt_traj```: the wrapper function for truncated alignment given expression matrix, alternative representation matrix and pseudotimes of cells in reference and query trajectories
 * ```preprocess_spring_given_rss.py```: the Python script to generate pseudocells and prepare the input folder for [SPRING](https://github.com/AllonKleinLab/SPRING)

## Related data
 * The raw data reported in this work can be found in ArrayExpress:
   * [E-MTAB-7552](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7552/): single-cell RNA-seq data based on 10x Genomics
   * [E-MTAB-8234](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8234/): single-cell RNA-seq data based on Fluidigm C1/Smart-seq2
   * [E-MTAB-8089](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8089/): single-cell ATAC-seq of human organoids
   * [E-MTAB-8043](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8043/): single-cell ATAC-seq of chimpanzee organoids
   * [E-MTAB-8083](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8083/): single-cell ATAC-seq of bonobo organoids
   * [E-MTAB-8087](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8087/): single-cell ATAC-seq of macaque organoids
   * [E-MTAB-8228](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8228/): the bulk ATAC-seq data
   * [E-MTAB-8230](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8230/): snRNA-seq data of the adult brain samples
   * [E-MTAB-8231](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8231/): bulk RNA-seq data of the adult brain samples
 * The processed scRNA-seq data is available in Mendeley Data: https://data.mendeley.com/datasets/z4jyxnx3vp/2
 * The data is available for exploration in [scApeX](https://bioinf.eva.mpg.de/shiny/sample-apps/scApeX/). If there is any issue accessing the app, please try the [mirror](https://scapex.ethz.ch).
 * The data is also browsable in [NeMO](https://nemoanalytics.org//index.html?share_id=37fd35aa&gene_symbol_exact_match=1)

## Related publication
The 10x Genomics based scRNA-seq data of human cerebral organoids in this work was reprocessed with CSS (Cluster Similarity Spectrum, https://github.com/quadbiolab/simspec) in [He et al. 2020 Genome Biol](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02147-4), which is also available in [biorxiv](https://www.biorxiv.org/content/10.1101/2020.02.27.968560v2).

## Contact
If you have any question related to the data and the analysis, please contact Dr. Zhisong He (zhisong.he(at)bsse.ethz.ch), Prof. Barbara Treutlein (barbara.treutlein(at)bsse.ethz.ch), Prof. J. Gray Camp (grayson.camp(at)iob.ch)
