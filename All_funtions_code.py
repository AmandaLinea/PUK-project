import seaborn as sns 
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 

# A function to create heatmaps showing the change in thermodynamic stability (ddg or dddg) for a PDB ID
# across all theoretic mutations across all positions in the sequence
def csv_to_two_heatmaps(csvfilepath, pdbid, valuecolumn, chainA=True, chainB=True, chainC=True, chainD=True, show_A =True, show_B=True, show_C=True, show_D=True):
    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(csvfilepath)
    df["pos"] = df["variant"].str[1:-1].astype("int")

    # Determine the aminoacid residues in the variant column of the csv file
    df["aa"] = df["variant"].str[-1]
    
    # Filter the dataframe to only include rows where the PDB id is the right one
    df = df[df['pdbid'] == pdbid]
    if chainA==True: 
        dfA = df[df['chainid'] == 'A'] 
        # Create a pivot table with the desired columns and index
        structure_A = pd.pivot_table(dfA.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_A = structure_A.reindex(dfA["aa"].unique())
    if chainB==True: 
        dfB = df[df['chainid'] == 'B']
        # Create a pivot table with the desired columns and index
        structure_B = pd.pivot_table(dfB.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_B = structure_B.reindex(dfB["aa"].unique())
    if chainC==True: 
        dfC = df[df['chainid'] == 'C']
        # Create a pivot table with the desired columns and index
        structure_C = pd.pivot_table(dfC.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_C = structure_C.reindex(dfC["aa"].unique())
    if chainD==True: 
        dfD = df[df['chainid'] == 'D']
        # Create a pivot table with the desired columns and index
        structure_D = pd.pivot_table(dfD.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_D = structure_D.reindex(dfD["aa"].unique())
    
    # Change the names for the value column so it is more understandable in the final title
    if valuecolumn == 'score_ml':
        valuecolumn = 'DDG'
    elif valuecolumn == 'score_ml_ddg_bind':
        valuecolumn = 'DDDG'

    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(16,8))

    if show_A==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_A , vmin=-2 , vmax=12, cmap='rainbow',ax=axs[0], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain A')
    if show_B==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_B , vmin=-2 , vmax=12 , cmap='rainbow',ax=axs[1], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain B')
    if show_C==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_C , vmin=-2 , vmax=12, cmap='rainbow',ax=axs[0], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain C')
    if show_D==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_D , vmin=-2 , vmax=12 , cmap='rainbow',ax=axs[1], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain D')
    plt.show()
    return fig



# A function to create a histogram of the number of variants with the different changes in stability
def csv_to_histogram(csvfilepath, pdbid, valuecolumn, minval, maxval, chainA=True, chainB=True, chainC=True, chainD=True):
    # Read the CSV file
    df = pd.read_csv(csvfilepath)
    df = df[df['pdbid'] == pdbid] 
    df = df[df[valuecolumn].between(minval, maxval)]
    df["pos"] = df["variant"].str[1:-1].astype("int")
    if chainA == True:
        dfA = df[df['chainid'] == 'A'] 
        # Create the histogram for the chain
        plt.hist(dfA[valuecolumn], bins=30, color='yellow', alpha=0.7, label='Chain A')
    if chainB ==True:
        dfB = df[df['chainid'] == 'B']
        # Create the histogram for the chain
        plt.hist(dfB[valuecolumn], bins=30, color='orange', alpha=0.4, label='Chain B')
    if chainC ==True:
        dfC = df[df['chainid'] == 'C'] 
        # Create the histogram for the chain
        plt.hist(dfC[valuecolumn], bins=30, color='purple', alpha=0.4, label='Chain C')
    if chainD==True:
        dfD = df[df['chainid'] == 'D']
        # Create the histogram for the chain
        plt.hist(dfD[valuecolumn], bins=30, color='red', alpha=0.4, label='Chain D')

    # Change the names for the value column so it is more understandable in the final title
    if valuecolumn == 'score_ml':
        valuecolumn = 'DDG'
    elif valuecolumn == 'score_ml_ddg_bind':
        valuecolumn = 'DDDG'

    # Set the title and axis labels
    plt.xlabel(f'{valuecolumn} values for {pdbid}')
    plt.ylabel('Frequency')
    plt.legend()

    return plt.show()

# functions used to get the utilized sequence for the protein or protein sequence

def every_nth(list, nth):
    # Use list slicing to return elements starting from the (nth-1) index, with a step of 'nth'.
    return list[nth - 1::nth]

def sequence_from_csv(csvfilepath, pdbid): 
    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(csvfilepath)
    df = df[df['pdbid'] == pdbid]
    df = df["variant"].str[0:1].astype("str")
    aa_sequence = every_nth(df, 20)
    aa_list = ''.join(aa_sequence)
    pd.set_option('display.max_rows', None)

    return print(f'The sequence of {pdbid} is: {aa_list}')

# This function flattens a nested list like [[0], [1], [2], [3], [4]] into just [0, 1, 2, 3, 4] 
def flatten(list):
    return [x for xs in list for x in xs]

# Converting the clinvar variant data into the same format as the RaSP csv file data
def variant_and_related_disease(clinvarcsvfile):
    df = pd.read_csv(clinvarcsvfile)
    df = df[[ "ref_aa", "pos_aa", "alt_aa", "ClinicalSignificance", "PhenotypeList",]]
    df = df.sort_values("pos_aa")
    df.replace("Ala", "A", inplace=True)
    df.replace("Arg", "R", inplace=True)
    df.replace("Asn", "N", inplace=True)
    df.replace("Asp", "D", inplace=True)
    df.replace("Cys", "C", inplace=True)
    df.replace("Glu", "E", inplace=True)
    df.replace("Gln", "Q", inplace=True)
    df.replace("Gly", "G", inplace=True)
    df.replace("His", "H", inplace=True)
    df.replace("Ile", "I", inplace=True)
    df.replace("Leu", "L", inplace=True)
    df.replace("Lys", "K", inplace=True) 
    df.replace("Met", "M", inplace=True) 
    df.replace("Phe", "F", inplace=True)
    df.replace("Pro", "P", inplace=True)
    df.replace("Ser", "S", inplace=True) 
    df.replace("Thr", "T", inplace=True) 
    df.replace("Trp", "W", inplace=True) 
    df.replace("Tyr", "Y", inplace=True) 
    df.replace("Val", "V", inplace=True) 
    
    df["variant"] = df["ref_aa"] + df["pos_aa"].astype(str) +df["alt_aa"]
    # If clinical significance and phenotype information is wanted as well, the hashtag can be removed from the next line
    #df = df[[ "variant" , "ClinicalSignificance", "PhenotypeList"]]
    df = df[[ "variant"]]
    list_of_data = df.values.tolist()
    return list_of_data

# A function to spit out a list of the clinical pathogenic variants and their corresponding theoretical ddg or dddg
def clinical_var_csv_to_list_with_stability(clinvarpathogeniccsvpath, clinvarbenigncsvpath, RaSPdatacsvpath, pdbid):
    df = pd.read_csv(RaSPdatacsvpath)
    df = df[df['pdbid'] == pdbid]

    # Empty lists for the pathogenic variants and benign variants
    path_variant = []
    path_ddg = []
    path_dddg = []
    benign_variant = []
    benign_ddg = []
    benign_dddg = []

    # Getting the clinvar variants as a list
    list1 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath))
    list2 = flatten(variant_and_related_disease(clinvarbenigncsvpath))

    for word in list1:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df[df['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid:
                 path_variant.append(row['variant'])
                 path_ddg.append(row['score_ml'])
                 path_dddg.append(row['score_ml_ddg_bind'])
    
    for word in list2:
        filtered = df[df['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid:
                 benign_variant.append(row['variant'])
                 benign_ddg.append(row['score_ml'])
                 benign_dddg.append(row['score_ml_ddg_bind'])

    #create empty dictionaries for the pathogenic and the benign variants
    pathogenic_dict = {}
    benign_dict = {}

    # Creating the pathogenic dictionary
    for i in range(len(path_variant)):
        pathogenic_dict[path_variant[i]] = (path_ddg[i], path_dddg[i])
    # Creating the benign dictionary
    for i in range(len(benign_variant)):
        benign_dict[benign_variant[i]] = (benign_ddg[i], benign_dddg[i])
    
    pathkeys = []
    pathvalues = []
    benignkeys = []
    benignvalues = []
    for key, value in pathogenic_dict.items():
        if value[1] >= 0.1:
            pathkeys.append(key)
            pathvalues.append(value)

    for key, value in benign_dict.items():
        if value[1] >= 0.1:
            benignkeys.append(key)
            benignvalues.append(value)

    dddg01 = pd.DataFrame({'Pathogenic variants': pathkeys, 'Pathogenic variant values' : pathvalues})
    dddgbenign01 = pd.DataFrame({'Benign variants': benignkeys, 'Benign variant values' : benignvalues})
    return print(f"The pathogenic variants of {pdbid} with dddg of more than 0.1: {dddg01} \nThe benign variants of {pdbid} with dddg of more than 0.1: {dddgbenign01}")


# A function to create scatter plots showing the clinical variants, either pathogenic or benign, for given PDB IDs
def scatterplot_clinical_var_ddg_and_dddg(clinvarpathogeniccsvpath1, clinvarbenigncsvpath1, RaSPdatacsvpath1, pdbid1, clinvarpathogeniccsvpath2, clinvarbenigncsvpath2, RaSPdatacsvpath2, pdbid2):
    df1 = pd.read_csv(RaSPdatacsvpath1)
    df1 = df1[df1['pdbid'] == pdbid1]
    df2 = pd.read_csv(RaSPdatacsvpath2)
    df2 = df2[df2['pdbid'] == pdbid2]

    # Empty lists for the pathogenic variants and benign variants
    path_ddg1 = []
    path_dddg1 = []
    benign_ddg1 = []
    benign_dddg1 = []
    path_ddg2 = []
    path_dddg2 = []
    benign_ddg2 = []
    benign_dddg2 = []

    # Getting the clinvar variants as a list
    list1 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath1))
    list2 = flatten(variant_and_related_disease(clinvarbenigncsvpath1))
    list3 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath2))
    list4 = flatten(variant_and_related_disease(clinvarbenigncsvpath2))

    for word in list1:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df1[df1['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid1:
                 path_ddg1.append(row['score_ml'])
                 path_dddg1.append(row['score_ml_ddg_bind'])
    
    for word in list2:
        filtered = df1[df1['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid1:
                 benign_ddg1.append(row['score_ml'])
                 benign_dddg1.append(row['score_ml_ddg_bind'])
    
    for word in list3:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df2[df2['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid2:
                 path_ddg2.append(row['score_ml'])
                 path_dddg2.append(row['score_ml_ddg_bind'])
    
    for word in list4:
        filtered = df2[df2['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid2:
                 benign_ddg2.append(row['score_ml'])
                 benign_dddg2.append(row['score_ml_ddg_bind'])
    
    figure = plt.figure(figsize = [8,8])
    plt.scatter(path_ddg1, path_dddg1, label=f'Pathogenic {pdbid1} chain A', color='red',  alpha=0.7)
    plt.scatter(benign_ddg1, benign_dddg1, label=f'Benign {pdbid1} chain A', color='lightsalmon', alpha=0.7)
    plt.scatter(path_ddg2, path_dddg2, label=f'Pathogenic {pdbid2} chain B', color='darkorchid',  alpha=0.7)
    plt.scatter(benign_ddg2, benign_dddg2, label=f'Benign {pdbid2} chain B', color='orchid', alpha=0.7)
    plt.xlabel('ddg')
    plt.ylabel('dddg')
    plt.title(f'Scatterplot of selected simple clinvar variants for {pdbid1} and {pdbid2} ddg and dddg values')
    plt.legend()
    plt.grid(True, linestyle='-', alpha=0.5)
    plt.show()
    return figure

