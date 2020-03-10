#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that collects RAMD data and determines what extra Replicas/ensembles need to be run/rerun.
Created on Fri Jan 31 14:44:55 2020

@author: andrewpotterton
"""


###Libraries
import numpy as np
from os import path
import pandas as pd

###Functions
def get_end_point(filename : str) -> int:
    '''
    Args: filename = is the RAMDLOG file.
    Returns: Endpoint of simulations or np.nan if did not finish. 
    '''
    end_point = np.nan #sets default returned value as np.nan value (so if the end of the file is not EXIT line then np.nan is returned)
    with open(filename, 'r') as f: #gets the last line of the RAMDLOG file
        lines = f.read().splitlines()
        last_line = lines[-1] 
    last_line = last_line.split(' ') #split last line by spaces
    if last_line[0] == 'EXIT:': #if line starts with EXIT, must be the exit line. Therefore simulation finished. 
        end_point = int(last_line[1]) #assigns the end point as the frame number.
    return end_point

def get_pandas_DF(ligand_name : str):
    '''
    Args: is the ligand name as a string (name of RAMDLOG_[CGS].log)
    Returns: pandas DataFrame of Ensemble and replica data for each ligand
    '''
    globals()['column_names'] = ['Ensemble1', 'Ensemble2', 'Ensemble3', 'Ensemble4', 'Ensemble5', 'Ensemble6', 'Ensemble7', 'Ensemble8']
    globals()['row_names'] = ['Rep1','Rep2','Rep3','Rep4','Rep5','Rep6','Rep7','Rep8','Rep9','Rep10','Rep11','Rep12','Rep13','Rep14','Rep15','Rep16','Rep17','Rep18','Rep19','Rep20','Rep21','Rep22','Rep23','Rep24','Rep25']
    Ensemble_values =[]
    for Ensemble_no in range(1,9):
        Replica_values =[]
        for Rep in range(1,26):
            if path.exists('Ensemble'+str(Ensemble_no)+'/Rep'+str(Rep)+'/RAMDLOG_'+ligand_name+'.log'):
                rep_value = get_end_point('Ensemble'+str(Ensemble_no)+'/Rep'+str(Rep)+'/RAMDLOG_'+ligand_name+'.log')
                Replica_values.append(rep_value)
            else:
                Replica_values.append(np.nan)
        Ensemble_values.append(Replica_values)
    #Transforms Ensemble_values list of lists to DataFrame. columns and index are inversed as in the next step, the dataframe is transposed
    DataFrame = pd.DataFrame(Ensemble_values,columns=row_names, index=column_names)
    DataFrame = DataFrame.T #Transposes the dataframe
    #Converting DataFrame values to ps from no. steps
    DataFrame=DataFrame.transform(lambda x: 2*x/1000)
    return DataFrame 

def bootstrapping(results : list) -> list:
    #RAMD definition of bootstrapping
    #Requires random and numpy packages
    #results is a list of numbers
    import random
    import numpy as np
    mean_results = []
    while len(mean_results) < 200:
        Results_list=[]
        #getting bootstrap samples 80% of sample length e.g. for 10 trajs bootstraps of 8 samples
        length = round((((len(results)+1)/10)*8), 0)
        length = int(length)
        for i in range(1, length):
            Results_list.append(random.choice(results))
        mean_results.append(float(sum(Results_list))/float(len(Results_list)))
    bootstrapped_mean= (float(sum(mean_results))/float(len(mean_results)))
    bootstrapped_sd = np.std(mean_results)
    #returns bootstrapped mean and standard deviation
    return [bootstrapped_mean, bootstrapped_sd]


###Main Code
#Uncomment and Run this code when setting up for the first time:
#noReplicas=pd.DataFrame(np.nan, index=column_names, columns=Lig_list)
#noReplicas.to_pickle('DataFrames/noReplicas.pkl')
#BootError=pd.DataFrame(np.nan, index=column_names, columns=Lig_list)
#BootError.to_pickle('DataFrames/Bootstrapped_errors.pkl')
if __name__ == "__main__":  
    median_list=[]
    Lig_list=['CGS','L3','L4','L7','Theo','XAC','ZM','CGS216','LUF5448','LUF5549','LUF5550','LUF5631','LUF5833','LUF5834','LUF5835','NECA','UK']
    #Opening noReplicas DF (opens file that describes the number of replicas that have been performed for each ligand in each ensemble)
    noReplicas = pd.read_pickle('DataFrames/noReplicas.pkl')
    #Opening BootError DF (holds current bootstrapped error values)
    BootError = pd.read_pickle('DataFrames/Bootstrapped_errors.pkl')
    
    for Lig in Lig_list:
        #Obtaining endpoints
        globals()[Lig]=get_pandas_DF(Lig) #Setting name of pandas DataFrame as string ligand name. Results in dataframes: CGS, L3, L4... with No. of frames completed in each replica and ensemble.
        median_list_Lig=[]
        for i in range(1,9): #looping through each ensemble 1-8
            #Working out medians
            median_list_Lig.append(globals()[Lig]['Ensemble'+str(i)].median()) #working out the median for each ensemble/lig
            #Working out no. replicas completed for each ligand and ensemble
            new_replica_count = 25 - globals()[Lig]['Ensemble'+str(i)].isna().sum() #working out no. of replicas with a value that is not np.nan (e.g. has a run value). This is the number of replicas performed so far.
            #if no. replicas has increased, work out a new bootstrapped error.
            if new_replica_count > noReplicas[Lig]['Ensemble'+str(i)]:
                error_list = globals()[Lig]['Ensemble'+str(i)].values.tolist() #Taking the ensemble which has increased no. of replicas and coverting values to a list
                error_list=[x for x in error_list if ~np.isnan(x)] #removing np.nan values from list
                boot_error = bootstrapping(error_list)[1] #calculating the bootstrapped error in no.steps
                #Saving bootstrapped error to DataFrame with other bootstrapped values
                BootError[Lig]['Ensemble'+str(i)] = boot_error
                #Saving new DataFrame with new replica count
            noReplicas[Lig]['Ensemble'+str(i)] = new_replica_count
        median_list.append(median_list_Lig)
        #Saving ligand dataframes to files for backup
        globals()[Lig].to_pickle('DataFrames/'+Lig+'.pkl')
    
    #Saving DFs to pickle files
    #Saving median data to DataFrame and then to file
    medianDF= pd.DataFrame(median_list, columns=column_names, index=Lig_list) #Setting median data to dataframe
    medianDF= medianDF.T #Transposes the dataframe
    medianDF.to_pickle('DataFrames/median.pkl') #Saving median pickle file
    #Saving bootstrapped error DF to pickle file
    BootError.to_pickle('DataFrames/Bootstrapped_errors.pkl')
    #Saving replica count to file
    noReplicas.to_pickle('DataFrames/noReplicas.pkl')
    
    #Print out which ligand ensemble need 5 more replicas performed. e.g. % SDreplicas > 50%
    SDreplicasAbove50 = BootError/medianDF*100 > 50 #Makes value True if SDreplica is above 50% of the median value
    print ("The following jobs need to be submitted as SDreplicas is above 50%:")
    for i in range(1,9): #looping through ensemble numbers
        print ("\nEnsemble" + str(i)+ ":")
        for row in SDreplicasAbove50: #looping through ligand names
            no_rep = noReplicas[row]['Ensemble'+str(i)]
            if (SDreplicasAbove50[row]['Ensemble'+str(i)]) and no_rep < 25: #if value is above 50% and no. replicas <25 
                print (f'{row} - Replicas {no_rep+1:.0f} to {no_rep+5:.0f}') #Print ligand name and name of replicas
    
    #Working out final tau_comp and final error.
    SD_comp=[]
    SD_replicas=[]
    tau_comp = medianDF.mean().values.tolist() #calculates mean value of the ensembles for each ligand and saves as a list
    SD_r = BootError.mean().values.tolist()
    for Lig in Lig_list: #looping through ligands to get SDreplicas for each ligand
        SD_replicas.append(np.nanstd(medianDF[Lig].values.tolist())) #finds the standard deviation ignoring np.nan values. 
    for SD_r_value, SD_replicas_value in zip(SD_r,SD_replicas): #looping through SD_r and SD_Replicas to find the maximum of the two values = SD_comp
        SD_comp.append(max(SD_r_value, SD_replicas_value))
    
    #Prints out ligands that SDr is above 50% so more replicas need to be performed (first checking SDreplicas is ok for that ligand)
    print ("\nEnsembles that need to be submitted: ")
    for Lig, err, mean in zip(Lig_list, SD_replicas, tau_comp): #looping through ligands associated SD_replicas and tau_comp
        Ensembles_finished = True #Sets the default to True (Determines whether anymore replicas needs to be completed and therefore it cannot yet be determined if more ensembles need to be run for that ligand)
        no_ensembles_peformed = 8 - medianDF[Lig].isna().sum()
        if err/mean > 0.5 and no_ensembles_peformed < 8: #if SD replicas/ tau comp is greater than 50% and number of ensembles performed is less than 8.
            for i in range(1,9): #looping through ensembles
                if (SDreplicasAbove50[Lig]['Ensemble'+str(i)]) or no_rep == 25:
                    Ensembles_finished = False
            if Ensembles_finished == True:#If this is true then continue
                print (f'{Lig} - Ensemble {no_ensembles_peformed+1:.0f}')
    
    
    #Making DF with final values and saving that DF
    Finalresults = pd.DataFrame([tau_comp, SD_comp], columns=Lig_list, index=['Tau_comp', 'SD_comp'])
    Finalresults.to_pickle('DataFrames/FinalResults.pkl')
