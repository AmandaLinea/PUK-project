from All_funtions_code import csv_to_two_heatmaps
from All_funtions_code import csv_to_histogram
from All_funtions_code import sequence_from_csv
from All_funtions_code import clinical_var_csv_to_list_with_stability
from All_funtions_code import scatterplot_clinical_var_ddg_and_dddg

denheatmap1_2 = csv_to_two_heatmaps(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN', 'score_ml',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
dendddgheatmap3_4 = csv_to_two_heatmaps(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN', 'score_ml_ddg_bind',chainA=False, chainB=False, chainC=True, chainD=True, show_A=False, show_B=False, show_C=True, show_D=True)

dendddghistogram = csv_to_histogram(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN', 'score_ml_ddg_bind',2,8)
denseq = sequence_from_csv(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv','5DEN')

den_clinvar_variants = clinical_var_csv_to_list_with_stability(r'ShinyClinVar_PAH_2024-01-04.csv',r'ShinyClinVar_PAH_benign_2024-01-04 (1).csv',r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN')

mutsalpha = scatterplot_clinical_var_ddg_and_dddg(r'ShinyClinVar_MSH2_2024-01-04.csv', r'ShinyClinVar_MSH2_benign_2024-01-04.csv', r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '2O8B',r'ShinyClinVar_MSH6_2024-01-04.csv', r'ShinyClinVar_MSH6_benign_2024-01-04.csv', r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '2O8B')
