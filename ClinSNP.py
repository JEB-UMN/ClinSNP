#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:54:56 2023

@author: administrator
"""
#Reminder: Please rename the folder from ClinSNP-main to ClinSNP after downloading. 
#If you did not do this before opening ClinSNP.py, please remove and reopen it in your interpreter after changing the folder name.

#Reminder: The ClinSNP folder has to be in your documents folder for it to function correctly.

import requests
import pandas as pd
import os
import gzip

# All of the code has been verified as functional, but some of it yields kickback notifications. This turns it off.
pd.options.mode.chained_assignment = None

# Get the user's home directory
user_home = os.path.expanduser("~")

# Change the working directory to the ClinSNP folder
new_working_directory = os.path.join(user_home, "Documents", "ClinSNP")
os.chdir(new_working_directory)

# The 'aligned_SNPs.csv' file contains all of an individual's SNPs that were annotated in the ClinVar database.
quick_true = input("If you already have an 'aligned_SNPs.csv' file from a previous run, and would like to perform your analysis on the same set of SNPs, please input the file here. Otherwise, type no. Please ensure that your desired SNP file is in the ClinSNP folder: ")

# Try to import aligned_SNPs file  
def file_search():
    if ".csv" in quick_true: # If a valid .csv extension was provided:
        if os.path.exists(quick_true): # If the aligned_SNPs file was found in the ClinSNP folder
            aligned_SNPs=pd.read_csv(quick_true,sep="\t",low_memory=False) # Read in the aligned_SNPs file
            return aligned_SNPs
        else: # If the file is not found in ClinSNP, return an error
            print("Error: 'aligned_SNPs' file not found. Please check whether your file is in the ClinSNP folder, then try again: ")
            return pd.DataFrame()
    else: # If the user doesn't input a valid .csv file, return an empty dataframe
        return pd.DataFrame()

while True:    
    aligned_SNPs = file_search()
    if not aligned_SNPs.empty or quick_true.lower() == "no":
        break
    else:
        quick_true = input("Invalid input. Please input a valid 'aligned_SNPs' file or type 'no' to skip: ")

    
if aligned_SNPs.empty:
    # Set the file path for the ClinVar database
    file_path = os.path.join(user_home, "Documents", "ClinSNP", "variant_summary.txt")

    # If the database has not yet been downloaded, download it
    if not os.path.exists(file_path):
        print("Downloading ClinVar database. This may take a few minutes.")
        # Set URL and destination paths
        url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        destination_filename = "variant_summary.txt.gz"
        unzipped_filename = "variant_summary.txt"
        
        # Construct full file paths
        destination_path = os.path.join(user_home, "Documents", "ClinSNP", destination_filename)
        unzipped_path = os.path.join(user_home, "Documents", "ClinSNP", unzipped_filename)
        
        # Retrieve the ClinVar database
        response = requests.get(url)
        
        # Save the downloaded database
        if response.status_code == 200:
            with open(destination_path, 'wb') as file:
                file.write(response.content)
            print("Database successfully downloaded.")
        
            # Unzip the downloaded database
            with gzip.open(destination_path, 'rb') as gz_file, open(unzipped_path, 'wb') as txt_file:
                txt_file.write(gz_file.read())
            print("Database successfully unzipped.")
        
        # If the download is unsuccessful based on status code, return error
        else:
            print(f"Failed to download database. Status code: {response.status_code}")
    
    # Import the user's SNP file
    patient_input = input("Please input a SNP .vcf file. Please ensure that you are in the correct working directory: ")
    print("Initiating alignment. This may take a few minutes.")
    snpfile = open(patient_input,"r")
    
    # Set up lists for the user's data
    sfinfo = []
    header = []
    chrom = []
    pos = []
    rsid = []
    ref = []
    alt = []
    qual = []
    filt = []
    info = []
    form = []
    bld = []
    zygous = []
    
    # Read the user's SNP data into a series of lists
    for line in snpfile:
        line = line.rstrip()
        if line[0:2] == "##": # Lines with SNP file info (sfinfo)
            line = line.lstrip("#")
            sfinfo.append(line)
        elif line[0:2] == "#C": # The first column in the SNP data header is #Chr
            header = line.split("\t")
        else:
            line = line.split("\t")
            line[0]=line[0].lstrip("chr")
            chrom.append(line[0])
            pos.append(int(line[1]))
            line[2]=line[2].replace("rs","")
            if line[2] == ".":
                line[2] = "-1"
            rsid.append(line[2])
            ref.append(line[3])
            alt.append(line[4])
            qual.append(line[5])
            filt.append(line[6])
            info.append(line[7])
            form.append(line[8])
            bld.append(line[9])
            if line[9][0:3] == "1/1": #1:1 indicates homozygous, 0:1 indicates heterozygous
                zygous.append("Homozygous")
            else:
                zygous.append("Heterozygous")
    # Identify the reference genome for SNP alignment        
    ref_pref = 0
    # Search the header for reference genome identifiers
    for index, entry in enumerate(sfinfo):
        if "grch37" in entry.lower() or "hg19" in entry.lower(): # If there's a GRCh37-associated identifier:
            ref_pref = 37        
            break
        elif "grch38" in entry.lower(): # If there's a GRCh38-associated identifier:
            ref_pref = 38
            break
    if ref_pref == 0: # If the reference couldn't be identified based on the header, ask for it directly from the user
        ask_ref = input("Reference genome not specified. Please define reference genome (type GRCh37 or GRCh38): ")
        while True:
            if "37" in ask_ref:
                ref_pref = 37
                break
            elif "38" in ask_ref:
                ref_pref = 38
                break
            else:
                ask_ref = input("Invalid reference genome. Please choose GRCh37 or GRCh38: ")
    
    # Create a dataframe from the ClinVar data
    df = pd.read_csv("variant_summary.txt", sep="\t",low_memory=False)
    
    # Create a more refined version of the ClinVar df that contains only the relevant information
    cvdf = df[["GeneSymbol","RS# (dbSNP)","Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","ClinicalSignificance","PhenotypeList","Assembly"]]
    
    # Add a column to cvdf indicating its source, and another for the zygosity (to be defined later)
    cvdf.loc[:, "Source"] = "ClinVar"
    cvdf.loc[:, "Zygosity"] = "N/A"
    
    # Create a new clean_cvdf containing only the SNPs with the correct reference genome alignment
    right_grch = []
    for index, row in cvdf.iterrows():
        if row["Assembly"] == "GRCh"+str(ref_pref):
            right_grch.append(index)
    clean_cvdf = cvdf.loc[right_grch].copy().reset_index(drop=True)
    
    #Remove the assembly column, since it adds to the computational requirement and isn't really relevant to the rest of the analysis
    #Also put the columns in the desired order
    cleaner_cvdf = clean_cvdf[["GeneSymbol","RS# (dbSNP)","Zygosity","Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","ClinicalSignificance","PhenotypeList","Source"]]
    
    # Assemble user data into a format that matches the clinvar dataframe
    user_df = pd.DataFrame({
        "GeneSymbol": "",
        "RS# (dbSNP)": rsid,
        "Zygosity": zygous,
        "Chromosome": chrom,
        "PositionVCF": pos,
        "ReferenceAlleleVCF": ref,
        "AlternateAlleleVCF": alt,
        "ClinicalSignificance": "",
        "PhenotypeList": "",
        "Source": "User"
    })
    
    # Combine the dataframes into a single datagrame containing all of the data
    user_cvdf = pd.concat([cleaner_cvdf,user_df],ignore_index=True)
    
    # Remove any entries that lack position data
    bad_pos = []
    for index, row in user_cvdf.iterrows():
        if row["PositionVCF"] == -1:
            bad_pos.append(index)
    clean_user_cvdf = user_cvdf.drop(bad_pos).reset_index(drop=True)
    
    # Sort the SNPs by chromosome, position, and nucleotide data
    user_cvdf_sorted_pos = clean_user_cvdf.sort_values(by=["Chromosome","PositionVCF","AlternateAlleleVCF","ReferenceAlleleVCF"],ignore_index=True)
    
    # Create a new dataframe containing aligned pairs of SNPs from ClinVar and the user
    SNP_pairs = []
    for index, row in user_cvdf_sorted_pos.iterrows():
        if index + 1 < len(user_cvdf_sorted_pos) and (
            (row["PositionVCF"] == user_cvdf_sorted_pos.at[index + 1, "PositionVCF"]) and
            (row["AlternateAlleleVCF"] == user_cvdf_sorted_pos.at[index + 1, "AlternateAlleleVCF"]) and
            (row["ReferenceAlleleVCF"] == user_cvdf_sorted_pos.at[index + 1, "ReferenceAlleleVCF"]) and
            (row["Source"] != user_cvdf_sorted_pos.at[index + 1, "Source"])
        ):
            SNP_pairs.append(index)
            SNP_pairs.append(index+1)
    pre_aligned_SNPs = user_cvdf_sorted_pos.loc[SNP_pairs].copy().reset_index(drop=True)
    
    # Assimilate the ClinVar and user data into a single row by writing zygosity data into a ClinVar row, then selecting it
    clin_index = []
    for index, row in pre_aligned_SNPs.iterrows():
        if index + 1 < len(user_cvdf_sorted_pos) and row["Source"] == "ClinVar":
            pre_aligned_SNPs.loc[index, "Zygosity"] = pre_aligned_SNPs.at[index+1, "Zygosity"]
            clin_index.append(index)
    
    # If an rsID is unknown, switch the value to N/A, which makes more sense conceptually than -1
    for index, row in pre_aligned_SNPs.iterrows():
        if index + 1 < len(user_cvdf_sorted_pos) and str(row["RS# (dbSNP)"]) == "-1":
            pre_aligned_SNPs.loc[index, "RS# (dbSNP)"] = "N/A"    
    
    # Create the aligned_SNPs folder from the ClinVar rows in the pre_aligned_SNPs dataframe
    aligned_SNPs = pre_aligned_SNPs.loc[clin_index].copy().reset_index(drop=True)
    
    # Export the aligned_SNPs file to the ClinSNP folder
    aligned_SNPs.to_csv("aligned_SNPs.csv",sep="\t",index=False)
    print("An 'aligned_SNPs.csv' file has been saved to the ClinSNP folder.")
    snpfile.close()

# Specify the path to create a ClinSNP Analysis folder
analysis_folder_path = os.path.join(user_home, "Documents", "ClinSNP", "ClinSNP Analysis")

# Check if the ClinSNP Analysis folder already exists
if not os.path.exists(analysis_folder_path):
    os.makedirs(analysis_folder_path) # Create the new folder

# Create a dataframe containing all of the SNPs that are not annotated as benign
non_benign_index = []
for index, row in aligned_SNPs.iterrows():
    if row["Source"] == "ClinVar" and "benign"  not in str(row["ClinicalSignificance"]).lower():
        non_benign_index.append(index)
non_benign_alignment = aligned_SNPs.loc[non_benign_index].copy().reset_index(drop=True)

# Establish whether the user wants a table of all non-benign SNPs
print("SNP alignment successful. Follow the remaining prompts to perform a phenotypic analysis. Output files will be saved in a 'ClinSNP Analysis' folder.")
non_benign_request = input("To rapidly identify phenotypes of interest, you can request a table containing all non-benign SNPs. Would you like to perform this analysis? Please input yes or no: ")
if non_benign_request.lower() == "no":
    print("Proceeding to individual phenotype analysis.")
else:
    non_benign_alignment.to_csv(os.path.join(user_home, "Documents", "ClinSNP", "ClinSNP Analysis","non_benign_SNPs.csv"),sep="\t",index=False)
    print("Analysis successful.\nProceeding to individual phenotype analysis")

# Establish whether the user wants all relevant SNPs, or just non-benign ones
function_prompt = input("Would you like to print out all SNPs for your phenotypes of interest, or focus only on non-benign ones? Please input 'all', 'non-benign', or 'both': ")
while True:
    if function_prompt.lower() == "all":
        prompt_all = aligned_SNPs.copy() # Fill prompt_all with all SNPs
        prompt_non_benign = pd.DataFrame() # Make prompt_non_benign empty
        break
    elif function_prompt.lower() == "non-benign":
        prompt_non_benign = non_benign_alignment.copy() # Fill prompt_non_benign with non_benign SNPs
        prompt_all = pd.DataFrame() # Make prompt_all empty
        break
    elif function_prompt.lower() == "both":
        prompt_all = aligned_SNPs.copy() # Fill prompt_all with all SNPs
        prompt_non_benign = non_benign_alignment.copy() # Fill prompt_non_benign with non_benign SNPs
        break    
    else:
        function_prompt = input("Invalid input. Please input 'all' or 'non-benign': ")

# Define a function to pull out SNP tables for individual phenotypes
def phenotype_analysis():
    phen = input("Please input a disease, tissue, or phenotype. Type 'no' to stop the analysis: ")
    if phen.lower() == "no":
        return phen # Skips the analysis, since the user said "no"
    else:
        if not prompt_all.empty: # If the user wanted a table containing all relevant SNPs
            phen_index = []
            for index, row in prompt_all.iterrows(): # Iterate through all SNPs
                if phen.lower() in str(row["PhenotypeList"]).lower(): # Search the phenotypic annotations for each SNP to identify ones that match the query
                    phen_index.append(index)
            # Download a table containing all SNPs pertaining to the desired phenotype
            filename = phen+"_all"+".csv"
            phen_snps = prompt_all.loc[phen_index].copy().reset_index(drop=True)
            phen_snps.to_csv(os.path.join(user_home, "Documents", "ClinSNP", "ClinSNP Analysis",filename),sep="\t",index=False)
        if not prompt_non_benign.empty: # If the user wanted a table containing only non-benign SNPs
            phen_index = []
            for index, row in prompt_non_benign.iterrows(): # iterate through the non-benign SNPs
                if phen.lower() in str(row["PhenotypeList"]).lower(): # Search the phenotypic annotations for each SNP to identify ones that match the query
                    phen_index.append(index)
            # Download a table containing all non_benign SNPs pertaining to the desired phenotype
            filename = phen+"_non-benign"+".csv"
            phen_snps = prompt_non_benign.loc[phen_index].copy().reset_index(drop=True)
            phen_snps.to_csv(os.path.join(user_home, "Documents", "ClinSNP", "ClinSNP Analysis",filename),sep="\t",index=False)
        print("Analysis successful.")
        return phen

# Run the individual phenotype analysis until the user inputs "no"
while True:
    run = phenotype_analysis()
    if str.lower(run) == "no":
        print("Thank you for using ClinSNP!")
        break

# Set the kickback notifications back to the standard settings
pd.options.mode.chained_assignment = "warn"
