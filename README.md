# ComputationalGenomics_FinalProject
Code for final project in Computational Genomics - Spring 2019


Files:

  MR_pipeline.py
  
    -contains class writen for running MR analysis
  
    -takes in exposure (through exposure file and column flag)
    
    -takes in geneotype data
    
    -selects best genetic predictor of exposure (must be corellated with pvalue < 0.01 to satisfy #1 assumption of MR)
    
    -takes in outcome data
    
    -checks assumption #2 of MR (gene in independent of outcome given exposure)
    
    -runs MR analysis with 2 stage least squares on outcomes where assumption holds
    
    -output results of MR analysis prediction of causalty for all outcomes where assumptions hold
    
    Input:
      -gene_f: file with genotype info, see genotype_data_EX.csv for formatting example
      -exp_f: file with exposure info, see exposure_data_EX.csv for formatting example
      -outcome_f: file with outcome info: see outcome_data_EX.csv for formatting example
    
    Optional Input:
      -pval (default=0.01): specifies how strong p value for relationship between gene and exposure must be to be valid
      -print_results (deault=True): flag to output text
      -print_graphs (defualt=False): flag to output graphs
      -check_missing (default=None): set to number to reject any outcomes with over that number of missing data points
      -exposure_col (default=0): set to index of data column (ignoring first ID column) that desired exposure is given in exp_f
      
  sample_run.py
    
    -code for simple run of MMR_pipeline on sample files

  __init__.py
    
    -empty file included so MR_pipeline can be imported

  genotype_data_EX.csv
    
    -sample file for genotype data in MR analysis
    
    -all used genotype files should be laid out in the same format, with unique IDs in the first column and data for allele count (0,1,2) for all relevent genes in other columns (each individual should have genes in same order)

  exposure_data_EX.csv
    
    -sample file for exposure in MR analysis
    
    -the default way to input exposure is in a file like this with 2 column, ID and exposure (as an int or float). However, if data is has multiple columns with one of them containing exposure, you can use the entire file and specify exposure_col=<appropriate index> in your MR pipeline run. If exposure_col is not specified, it will default to the first column (after ID)

  outcome_data_EX.csv
    
    -sample file for outcome in MR analysis
    
    -like genotype and exposure, all outcome files should start with individual ID, then have all outcomes in following columns
    
    -data should be in ints or floats for meaningful analysis


Requirements:

-Python 3 with modlues:
  -pyreadr
  -math
  -scipy.stats
  -csv
  -matplotlib.pyplot
  -numpy


How to run on example data:
1. download folder
2. run code in sample_run.py


How to run on larger files:
1. download folder
2. set up MR analysis tool in same way as sample_run 
  
  -create an instance with your genotype, exposure and outcome files
  
  -all variables should be converted to ints/floats
  
  -all files should have data in rows with first column as individual identifier, as in example
  
  -set appropriate flags
  
  -run tool.full_analysis()


How to access Alzheimer's files used for my analysis:
  
  -all data I used is available at http://adni.loni.usc.edu/ by creating an account
  
  -I used adnimerge files from the adnimerge data table in https://adni.bitbucket.io/ 
  
