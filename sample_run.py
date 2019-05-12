from MR_pipeline import MR_pipeline

def main():
    tool = MR_pipeline('genotype_data_EX.csv', 'exposure_data_EX.csv', 'outcome_data_EX.csv')
    tool.full_analysis()

main()