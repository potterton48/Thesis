#Always run this section first before running any other section
###Libraries###
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Descriptors
import openbabel
import os
import os.path
from os import path
import collections
import subprocess
import glob
from itertools import combinations
import shutil
#Import my functions
from KineticDataFunctions import get_descriptors
from KineticDataFunctions import generateconformations
from KineticDataFunctions import get_atom_distance
###Main Code###
#Loading in data, each sheet is a separate dataframe
kinetic_data = pd.ExcelFile('Kinetic Data.xlsx')
FiveHT2A = pd.read_excel(kinetic_data, '5-HT2A')
FiveHT2B = pd.read_excel(kinetic_data, '5-HT2B')
A1_Adenosine = pd.read_excel(kinetic_data, 'A1 Adenosine')
A2A_Adenosine = pd.read_excel(kinetic_data, 'A2A Adenosine')
A2B_Adenosine = pd.read_excel(kinetic_data, 'A2B Adenosine')
A3_Adenosine = pd.read_excel(kinetic_data, 'A3 Adenosine')
Alpha2A_Adrenoceptor = pd.read_excel(kinetic_data, 'Alpha2A-Adrenoceptor')
Alpha2C_Adrenoceptor = pd.read_excel(kinetic_data, 'Alpha2C-Adrenoceptor')
Beta1_Adrenoceptor = pd.read_excel(kinetic_data, 'Beta1-adrenoceptor')
Beta2_Adrenoceptor = pd.read_excel(kinetic_data, 'Beta2-adrenoceptor')
Cannabinoid_CB1 = pd.read_excel(kinetic_data, 'Cannabinoid CB1')
Cannabinoid_CB2 = pd.read_excel(kinetic_data, 'Cannabinoid CB2')
CCR2 = pd.read_excel(kinetic_data, 'CCR2')
CRF_1 = pd.read_excel(kinetic_data, 'CRF-1')
H1_Histamine = pd.read_excel(kinetic_data, 'H1-Histamine')
M2_Muscarinic_acetylcholine = pd.read_excel(kinetic_data, 'M2 Muscarinic acetylcholine')
M3_Muscarinic_acetylcholine = pd.read_excel(kinetic_data, 'M3 Muscarinic acetylcholine')
mGlu2 = pd.read_excel(kinetic_data, 'mGlu2')
Mu_Opioid_Receptor = pd.read_excel(kinetic_data, 'Mu-Opioid Receptor')
NK1 = pd.read_excel(kinetic_data, 'NK1')
Orexin_2 = pd.read_excel(kinetic_data, 'Orexin 2')

#Getting chemical properties of data from SMILES strings and appending to original dataframe
FiveHT2A = get_descriptors(FiveHT2A)
FiveHT2B = get_descriptors(FiveHT2B)
A1_Adenosine = get_descriptors(A1_Adenosine)
A2A_Adenosine = get_descriptors(A2A_Adenosine)
A2B_Adenosine = get_descriptors(A2B_Adenosine)
A3_Adenosine = get_descriptors(A3_Adenosine)
Alpha2A_Adrenoceptor = get_descriptors(Alpha2A_Adrenoceptor)
Alpha2C_Adrenoceptor = get_descriptors(Alpha2C_Adrenoceptor)
Beta1_Adrenoceptor = get_descriptors(Beta1_Adrenoceptor)
Beta2_Adrenoceptor = get_descriptors(Beta2_Adrenoceptor)
Cannabinoid_CB1 = get_descriptors(Cannabinoid_CB1)
Cannabinoid_CB2 = get_descriptors(Cannabinoid_CB2)
CCR2 = get_descriptors(CCR2)
CRF_1 = get_descriptors(CRF_1)
H1_Histamine = get_descriptors(H1_Histamine)
M2_Muscarinic_acetylcholine = get_descriptors(M2_Muscarinic_acetylcholine)
M3_Muscarinic_acetylcholine = get_descriptors(M3_Muscarinic_acetylcholine)
mGlu2 = get_descriptors(mGlu2)
Mu_Opioid_Receptor = get_descriptors(Mu_Opioid_Receptor)
NK1 = get_descriptors(NK1)
Orexin_2 = get_descriptors(Orexin_2)

#%%
##########################################
####         Make Receptor Models     ####
##########################################

#Rotating the Receptor PDB to X-Y plane for membrane insertion using Chimera. This is for creation of the membrane.
#Also adds Hydrogen to the receptor in preparation for packmol-memgen
Chimerapath='/Applications/Chimera.app/Contents/MacOS/chimera' #path for chimera executable
Reduce_command = "/Users/andrewpotterton/Documents/amber18/bin/reduce -NOFLIP -Quiet Receptors/Rotated/filename.pdb > Receptors/Rotated/filename_H.pdb"
for filename in os.listdir('Receptors/'): #loops through files in the Receptor folder
    if filename.endswith(".pdb"): #only selects PDB files
        subprocess.call([Chimerapath, '--nogui', "Receptors/"+filename, 'chimera_rotate.py']) #runs rotate
        subprocess.call(['mv','Receptors/Rotated/temp_xy.pdb', 'Receptors/Rotated/'+filename]) #renames the output of that rotation (temporary name) to the filename
        Reduce_command_file = Reduce_command.replace('filename', filename[:-4]) #replaces run command with actual names of PDB files for input and output
        subprocess.call(Reduce_command_file, shell=True) #Runs add Hydrogen command. Has to be shell=True for some strange reason
#%%
#Building membrane around receptor
Location_to_run_packmol = '/Users/andrewpotterton/Documents/amber18/lib/python2.7/site-packages/' #location where to run packmol-memgen to get the correct scripts to load
amber_location_relative = '../../../bin/' #path of amber relative to the "Location_to_run_packmol"
location_of_kinetic_data= '../../../../University/PhD/KineticData/'
os.chdir(Location_to_run_packmol)
for filename in os.listdir(location_of_kinetic_data+'Receptors/Rotated/'): #loops through files in the Rotated folder
    if filename.endswith(".pdb") and filename != 'XYplaneReceptor.pdb' and filename.endswith("_H.pdb") == False: #only selects PDB files XYplaneReceptor is the example receptor doesn't need to be in a membrane. also avoids protonated files
        subprocess.call(['../../../miniconda/bin/python2.7',amber_location_relative+'packmol-memgen', '--lipids', 'DPPC', '--pdb', location_of_kinetic_data+'Receptors/Rotated/'+filename[:-4]+'_H.pdb', '--preoriented', '--output', location_of_kinetic_data+'Receptors/Membrane/'+filename[:-4]])

os.chdir('../../../../University/Phd/KineticData/')

#%% Copying Receptor-membrane complex to MD folder to run EQ.
#Also adds CONECT lines for disulphide bonds.
#Creates Structure folder for .pdb and .prmtop files, EQ folder for equilibration runs and Rep1-10 for production MD runs
if os.path.exists('MD/') == False: #if MD folder doesn't exist:
    os.mkdir('MD/') #Make MD folder
Receptor_names =['5HT2A', '5HT2B','A1', 'A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'Beta2','CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2', 'M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
listpdbs = ['5HT2A.pdb','5HT2B.pdb','A1_FullyActive_charged.pdb','A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active_charged.pdb', 'B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive_charged.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb','Mu_opoid_Inactive_charged.pdb','Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
for Receptor, filename in zip(Receptor_names, listpdbs):
    #Making Folders for MD
    if os.path.exists('MD/'+Receptor) == False: #if Ligand subfolder doesn't exist:
        os.mkdir('MD/'+Receptor) #Make ligand subfolder in MD folder
        os.mkdir('MD/'+Receptor+'/Structures/')
        os.mkdir('MD/'+Receptor+'/EQ/')
        for i in range(1, 11):
            os.mkdir('MD/'+Receptor+'/Rep'+str(i))
    #Copying over PDB of Receptor-membrane complex - correcting HIS and CYX residues
    savereceptorlines=[] #list of lines to save
    CYXatomslist = [] #list of CYX atoms
    first_ter = True #Hasn't encoured TER line in PDB file
    #Opening file and saving lines. Adding the CYX CONECT lines and removing HIS HD1 atoms
    with open("Receptors/Membrane/"+filename) as fname:
        for line in fname: #loops through PDB lines
            #Adds CYX CONECT LINEs
            if line[17:20] == 'CYX' and line[13:15] == 'SG': #checking if residue name equals CYX and atom name equals SG (Sulphur atom)
                CYX = [line[6:11].strip(), line[30:39].strip(), line[39:47].strip(),line[47:54].strip()] #saves atom number, x, y and z coordinates in a list
                CYXatomslist.append(CYX) #adds details to a list of a list
            if line.startswith('TER') and first_ter == True: #checks for the first TER value
                comb = combinations(CYXatomslist, 2) #does combination values of the all atoms in the list
                for atom1, atom2 in comb:
                    if get_atom_distance(atom1, atom2) < 2.4: #if atom distance between two S atoms is under 2.4 angstroms must be disulphide bond
                        #making CONECT line record for PDB
                        if len(atom1[0]) == 3: #if atom number is in the hundreds,  need an extra space in CONECT line
                            CONECTline = 'CONECT '+atom1[0]+' '+atom2[0] + '\n'
                        else:
                            CONECTline = 'CONECT'+atom1[0]+' '+atom2[0] + '\n'
                        savereceptorlines.append(CONECTline)
                        first_ter = False
            if line[17:20] == 'HIS' and line[13:16] == 'HD1': #avoiding writing HIS HD1 Atom as this breaks tleap
                pass
            else:
                savereceptorlines.append(line) #saves the PDB line
    with open('MD/'+Receptor+'/EQ/'+filename, 'w') as fname:
        fname.writelines(savereceptorlines)
#%% Paramerisation of Membrane-receptor complex
Amber_location = '/Users/andrewpotterton/Documents/amber18/bin/'
Receptor_names =['A1', 'A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'Beta2','CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2', 'M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
listpdbs = ['A1_FullyActive_charged.pdb','A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active_charged.pdb', 'B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive_charged.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb','Mu_opoid_Inactive_charged.pdb','Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
for Receptor, filename in zip(Receptor_names, listpdbs):
    savelines = []
    with open("Scripts/default_Protein.in") as fname: #opening file to read
        for line in fname:
            #replaces $$$ with Receptor location and filename
            newline = line.replace("$$$", 'MD/'+Receptor+'/EQ/'+filename[:-4])
            savelines.append(newline)
    #Writing modified tleap script
    with open("Scripts/Protein.in", 'w') as outfile:
        outfile.writelines(savelines) #writing all lines that were saved
    subprocess.call([Amber_location+'tleap', '-f', 'Scripts/Protein.in'])

#%% Determining PBCs for each receptor
vmd_path = '/Applications/VMD 1.9.3.app/Contents/vmd/vmd_MACOSXX86'
listpdbs = ['A1_FullyoActive_charged.pdb','A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active_charged.pdb', 'B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive_charged.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb','Mu_opoid_Inactive_charged.pdb','Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
Receptor_names =['A1', 'A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'Beta2','CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2', 'M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
for Receptor, filename in zip(Receptor_names, listpdbs):
    #carrying out PBC script, which writes minmax XYZ coordinates to txt file called PBC.txt in the Scripts folder
    subprocess.call([vmd_path, '-e', 'Scripts/PBC_box.tcl', '-pdb', 'MD/'+Receptor+'/EQ/'+filename, '-dispdev', 'text'])
    with open('Scripts/PBC.txt') as fname: #open that PBC.txt file
        minmax = fname.readline() #read data from file
    minmax = minmax.strip() #removing newline character
    minmax = minmax.replace('{', '') #removing {}
    minmax = minmax.replace('}', '')
    minmax = minmax.replace(' ', ',') #replacing spaces with commas
    minmaxlist = minmax.split(',') #splitting string by comma to get the individual values in a list
    x = int(abs(float(minmaxlist[0])-float(minmaxlist[3]))-9) #calculating PBC from min max values.
    y = int(abs(float(minmaxlist[1])-float(minmaxlist[4]))-9) #Takes max from the minimum and gets the absolute value. Rounds to the nearest integer.
    z = int(abs(float(minmaxlist[2])-float(minmaxlist[5]))-5) #removes 9 Angstroms from X, Y coordinates and 5 from the Z coordinates
    #Writing PBC conditions to a file called PDBnamePBC.txt
    with open('MD/'+Receptor+'/Structures/'+filename[:-4]+'PBC.txt', 'w') as fname:
        fname.write(str(x)+' '+ str(y) + ' '+ str(z))
#%% Generation of EQ scripts and EQ restraints files
vmd_path = '/Applications/VMD 1.9.3.app/Contents/vmd/vmd_MACOSXX86'
listpdbs = ['A1_FullyActive_charged.pdb','A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active_charged.pdb', 'B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive_charged.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb','Mu_opoid_Inactive_charged.pdb','Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
Receptor_names =['A1', 'A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'Beta2','CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2', 'M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
for Receptor, filename in zip(Receptor_names, listpdbs):
    #Making the first restraint files for EQ
    subprocess.call([vmd_path, '-e', 'Scripts/step1_con.tcl', '-parm7', 'MD/'+Receptor+'/EQ/'+filename[:-4]+'_PM.prmtop','-pdb', 'MD/'+Receptor+'/EQ/'+filename[:-4]+'_PM.pdb', '-dispdev', 'text'])
    subprocess.call(['mv','s6.1_con.ref', 'MD/'+Receptor+'/EQ/'+filename[:-4]+'s6.1_con.ref']) #moving constraints file
    #Obtaining PBC conditions to a file called PDBnamePBC.txt
    with open('MD/'+Receptor+'/Structures/'+filename[:-4]+'PBC.txt') as fname:
        PBC = fname.readline()
        PBClist = PBC.split(' ') #splitting string by comma to get the individual values in a list
    for i in range(1, 7):
        #Modifying the EQ scripts
        savelines =[]
        with open('Scripts/step6.'+str(i)+'_production.inp') as fname:
            for line in fname:
                newline = line.replace('£££', filename[:-4])
                newline = newline.replace('%X%', PBClist[0])
                newline = newline.replace('%Y%', PBClist[1])
                newline = newline.replace('%Z%', PBClist[2])
                savelines.append(newline)
        with open('MD/'+Receptor+'/EQ/'+filename[:-4]+'step6.'+str(i)+'_production.inp','w') as fname:
            fname.writelines(savelines)
#%% #Making .sh run scripts for HPC
#Copying files to HPC
listpdbs = ['A1_FullyActive_charged.pdb','A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active_charged.pdb', 'B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive_charged.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb','Mu_opoid_Inactive_charged.pdb','Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
Receptor_names =['A1', 'A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'Beta2','CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2', 'M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
for Receptor, filename in zip(Receptor_names, listpdbs):
    savelines =[]
    #modifying GRACE script with receptor name
    with open('Scripts/default_GRACE.sh') as fname:
        for line in fname:
            newline = line.replace('$NAME$', filename[:-4])
            newline = newline.replace('$location$', Receptor)
            savelines.append(newline)
    #Saving that script
    with open('Scripts/HPC/'+filename[:-4]+'Grace.sh', 'w') as fname:
        fname.writelines(savelines)
############################################
###Copy HPC folder and MD folder to GRACE###
###              Run HPC                 ###
############################################
#Make a PDB snapshot of EQ version Save in Receptors/PM_EQ/ folder - automate this
#Make charged version of PDB snapshot if necessary
#%%Convert PDB to Amber format (Adding TER lines)
#Add TER after OXT atom
#Change Residue numbering for Membrane - PAPCPA
#Add TER line after every Cl- atom
#Add TER line after every H2 atom in WAT
directory = 'Receptors/PM_EQ/'
for filename in os.listdir(directory):
    savelines = [] #saves line
    if filename.endswith("Snapshot.pdb"):
        membrane_count =0 #counts the number of membrane residue numbers
        membrane_res = 0 #the membrane residue number
        letters = ['B','C','D','E','F','G','H','I'] #Letters for Chain IDs for water
        letter = 'A' #First Chain ID for Water
        water_count = 0
        with open(directory+filename) as fname: #opening file
            for line in fname: #looping through lines
                #Checking if ATOM is in a membrane atom
                if line[17:19] == 'PA' or line[17:19] =='PC': #After the protein
                    if membrane_res == 0: #setting the first residue number of the first membrane atom
                        membrane_res = line[23:26]
                        newmembrane_res = membrane_res #newmembrane_res is the new numbering for the membrane, only goes up every three membrane molecules
                    elif line[23:26] != membrane_res: #if the residue number changes (e.g. goes up by one)
                        membrane_res =line[23:26] #change the numbering to go up
                        membrane_count = membrane_count +1
                    if membrane_count == 3: #reached the next PAPCPA molecule
                        membrane_count = 0 #reset the counter
                        newmembrane_res = str(int(newmembrane_res) +1) #increase this by one, convert back to string
                        savelines.append('TER\n') #add a terline
                    newline = line[0:23]+newmembrane_res+line[26:] #saving line with new residue number
                    savelines.append(newline)
                elif membrane_count != 0: #not membrane atoms, adds first TER after membrane
                    membrane_count = 0 #resets membrane count
                    savelines.append('TER\n')
                    savelines.append(line)
                elif line[17:20] == 'WAT':
                    if line[13:20] == 'O   WAT':
                        water_count = water_count +1
                    if water_count == 9999:
                        water_count = 1
                        letter = letters.pop(0)
                    newline = line[:21]+letter+'{:4d}'.format(water_count)+line[26:]
                    savelines.append(newline)
                else:
                    savelines.append(line)
                #adding TER lines
                if line[13:16] == 'OXT': #After the protein
                    savelines.append('TER\n')
                elif line[13:20] == 'H2  WAT': #after each water molecule
                    savelines.append('TER\n')
                elif line[13:20] == 'Cl- Cl-': #after each chlorine ion
                    savelines.append('TER\n')
        #saving lines as new file
        f= open(directory+'AMBER/'+filename,"w+") #opening file to write
        f.writelines(savelines) #Writes saved lines (protein PDB lines) to a file
        f.close() #closing file
#%% Making a receptor only version of the Membrane-receptor PDB for docking
"""
Takes the protein-membrane PDB and saves another copy with just the protein for use in docking
"""
#Looping through receptors .pdb files in receptor directory
directory = 'Receptors/'
for filename in os.listdir(directory+'PM_EQ/Amber/'):
    write = True #sets saving lines on
    savelines = [] #saves line
    if filename.endswith(".pdb"):
        with open(directory+'PM_EQ/Amber/'+filename) as fname: #opening file
            for line in fname: #looping through lines
                if line.startswith('TER'): #the first TER should be after the protein after which point writing lines will be turned off
                    write = False #turns off writing lines
                if write: #if write is true, it will save the lines
                    savelines.append(line)
        #saving lines as new file
        f= open(directory+'Membrane/WithoutMembrane/'+filename,"w+") #opening file to write
        f.writelines(savelines) #Writes saved lines (protein PDB lines) to a file
        f.close() #closing file
#%%
#Generation of 3D confomations - done for 5-HT2A and 5-HT2B
for index, row in FiveHT2B.iterrows(): #This is code for just H1 histamine change this line to another receptor pandas dataframe
    #Generating conformers of every ligand for each receptor
    m = rdkit.Chem.MolFromSmiles(row['SMILES'])
    generateconformations(m, row['Ligand'], 'H1') #H1 is the folder name of the receptor for which its ligands are searching conformations

#%% Find COM of ligand in rotated/without membrane Receptor file
"""
Issue is the ligand is not present in the model receptor. So:
    Take the crystal structure, matchmaker onto the receptor model Done
    Find COM of ligand in rotated crystal structure
"""
receptor_name='5HT2B_EQ_Snapshot' #name of receptor model
crystalname='4ib4' #name of crystal structure
Ligand_name = 'ERM' #name of ligand in PDB crystal structure
Chain = 'A' #leave as blank string if no chain ID (need a . before Chain ID)
Chimerapath='/Applications/Chimera.app/Contents/MacOS/chimera' #path for chimera executable
savelines =[] #
with open("chimera_rotate_def.py") as fname: #opening file to read
        for line in fname:
            #replaces ??? with ligand name
            newline = line.replace("???", Ligand_name)
            newline = newline.replace("XXX", Chain) #replaces XXX with chain ID
            savelines.append(newline)
    #Writing modified tleap script
with open("chimera_rotate2.py", 'w') as outfile:
    outfile.writelines(savelines) #writing all lines that were saved
subprocess.call([Chimerapath, '--nogui', 'Receptors/Membrane/WithoutMembrane/'+receptor_name+'.pdb', 'Receptors/Crystal_structures/'+crystalname+'.pdb','chimera_rotate2.py']) #runs matchmaker and finds COM
os.remove('Receptors/Crystal_structures/temp_MM.pdb') #deletes temperary file made during production
##############################################################
###                PERFORM DOCKING IN GOLD                 ###
###           Use EQ without membrane snapshot             ###
#### Receptors/Membrane/WithoutMembrane/*EQ_Snapshot.pdb  ####
##############################################################

#%% Moving lowest scoring docking results into dedicated folder
"""
1. Finds the location of the best scoring docking result
2. Converts this confomer to PDB file and saves to new folder "Parameterisation"
3. Preps this PDB file for parameterisation (removes CONECT lines with lone pair atoms)
4. Add ATOM names to PDB file automatically
5. Opens file up with Chimera and saves to ensure CONECT lines are present
"""
Chimerapath='/Applications/Chimera.app/Contents/MacOS/chimera'
Receptors = [A1_Adenosine, A2A_Adenosine, Beta1_Adrenoceptor, Beta2_Adrenoceptor, Cannabinoid_CB1, Cannabinoid_CB2, CCR2, CRF_1, H1_Histamine, M2_Muscarinic_acetylcholine, M3_Muscarinic_acetylcholine, Mu_Opioid_Receptor, Orexin_2] #setting up list of panda dataframes
Receptor_names =['A1', 'A2A', 'Beta1', 'Beta2', 'CB1','CB2','CCR2', 'CRF1', 'H1','M2','M3', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
obConversion = openbabel.OBConversion() #setting up conversion
obConversion.SetInAndOutFormats("mol2", "pdb") #setting input as mol2 and output as PDB file
mol = openbabel.OBMol() #setting up conversion
for Receptor, Receptor_name in zip(Receptors, Receptor_names): #loops through receptors
    for index, row in Receptor.iterrows(): #Loops through non-repeated values of the database
        active=False
        inactive=False
        ligand = row['Ligand'] #Takes the ligand name
        #Finding folder where docking scores are located (this depends whether there is an active/inactive structure)
        folder_name="Conformers/"+Receptor_name+"/Docking/GoldScore/" #folder name where docking is located
        if path.exists(folder_name+ligand+"_m1.rnk") == False: #if can't find file, then it must either be active or inactive
            new_folder_name=folder_name+'Active/'#setting new to active
            active = True
            if path.exists(new_folder_name+ligand+"_m1.rnk") == False: #testing active, if false the ligand must be in the inactive folder
                new_folder_name=folder_name+'Inactive/' #setting folder to inactive
                inactive=True
                active =False
            folder_name=new_folder_name #saving the new_folder_name to overwrite folder_name
        #Determing the docking number which has the lowest energy value by opening rank value
        with open(folder_name+ligand+"_m1.rnk") as fname: #opening rank file with ranking information
            fname.readline() #removing header lines
            fname.readline()#removing header lines
            fname.readline()#removing header lines
            fname.readline()#removing header lines
            toprank=fname.readline() #read the ranking line, top line has the highest scoring molecule
            toprank=toprank.strip() #removing whitespace
            toprank = toprank[0:2].strip() #getting string number of the top ranking molecule - can be from 1-10
        lowest_file = folder_name+"gold_soln_"+ligand+"_m1_"+toprank+".mol2" #name and path of .mol2 file - best scoring docking
        #Getting name of output folder
        dest_folder = "Conformers/"+Receptor_name+"/Parameterisation/"
        if active == True:
            dest_folder = dest_folder+"Active/"
        elif inactive == True:
            dest_folder = dest_folder+"Inactive/"
        #Converting mol2 to pdb file and saving in appropriate parameterisation folder
        obConversion.ReadFile(mol, lowest_file)
        obConversion.WriteFile(mol, dest_folder+ligand+".pdb") #renaming PDB as ligand name + pdb
        #Removing lone pair atoms Atom name "Xx" and CONECT/MASTER lines and header lines
        #Also fills in ATOM names into PDB file
        write_lines = [] #list of lines to be saved
        element_namelist = [] #List of atom names
        lonepairs = []
        with open(dest_folder+ligand+".pdb") as fname: #opening file
            for line in fname: #looping through lines
                if line.startswith("ATOM") == True: #if line starts with ATOM then it will check if atom name is Xx = lone pair
                    if line[76:78].strip() != "Xx": #if ATOM record is not "Xx" e.g. is not lone pair then it will save that line
                        """
                        PDB formatting states that ATOM name is right justified in the 13-14th column
                        and atom numbering is left justified 15-16th column
                        """
                        element_name = line[76:78].strip() #gets element name
                        element_namelist.append(element_name) #adds element name to list, to store all elements
                        element_count = int(collections.Counter(element_namelist)[element_name])  #provides the count of that individual element
                        new_atom_name =element_name+str(element_count) #Atom name = element name + element count
                        if len(element_name) == 1: #if element name is a single number
                            if len(str(element_count)) == 1: #if element count is a single digit
                                new_atom_name = ' ' + new_atom_name + ' ' #adds a space before and after
                            else: #must be a two digit number
                                new_atom_name = ' ' + new_atom_name #adds space before only
                        elif len(str(element_count)) == 1: #else element name must be two letters, checks if number is a single digit
                            new_atom_name = new_atom_name + ' ' #if so add a space at the end
                        #Need to write code that justifies ATOM NAME different for two letter atom names and different lengths of string
                        newline=line[0:12]+new_atom_name+line[16:17]+"LIG  "+line[22:] #removing chain ID and replacing residue name to LIG to make it consistant
                        #Change line to newline when code is complete
                        write_lines.append(newline) #if starts with ATOM and atom name does not equal Xx then line will be saved
                    else: #if element is Xx = lone pair
                        lonepairs.append(line[8:11].strip()) #Adds atom number to list. This lonepairs list contains the atom numbers of the lone pair atoms
                elif line.startswith("CONECT") == True:
                    newline=line
                    if (any(x in lonepairs for x in line.split())) == False: #checks if CONECT line has number from lone pair list
                        #write code that removes CONECT line if it includes atom id that will be removed
                        write_lines.append(newline) #if starts with CONECT line will be saved
        with open(dest_folder+ligand+".pdb", 'w') as outfile: #rewriting file
            outfile.writelines(write_lines) #writing all lines that were saved
            outfile.write("END") #adding END line
        subprocess.call([Chimerapath, '--nogui', dest_folder+ligand+".pdb", 'chimera.py'])
        subprocess.call(['mv', 'temp.pdb', dest_folder+ligand+".pdb"])

###########################################
# Replace crystal structure ligands poses #
# Matchmaker crystal structure to Membrane-receptor PDB structure
        #Select Ligand
        #Save PDB of ligand only relative to membrane-receptor
        #Add Hydrogens to PDB
        #Rename PDB to LIG, remove chain ID and chain ResID if nessessary.
###########################################
#%% Parameterisation of the ligands
Amber_location = '/Users/andrewpotterton/Documents/amber18/bin/'
Receptor = '5HT2A'
for index, row in FiveHT2A.iterrows():
    ligand_name = row['Ligand'] #Takes the ligand name
    #Determing charge of ligand
    charge = '0' #setting charge to 0
    if ('+' in row['SMILES']) == True: #if there a + charge in the SMILES string
        charge = '+1' #charge of the ligand is +1
    #Determines whether ligand is in the Active, Inactive folder or neither
    Activity = ''
    if path.exists('Conformers/'+Receptor+'/Parameterisation/Active/'+ligand_name+'.pdb') == True:
        Activity = 'Active/'
    elif path.exists('Conformers/'+Receptor+'/Parameterisation/Inactive/'+ligand_name+'.pdb') == True:
        Activity = 'Inactive/'
    subprocess.call([Amber_location+'antechamber', '-i', 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.pdb', '-fi', 'pdb', '-o', 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.mol2', '-fo', 'mol2','-c', 'bcc', '-nc', charge, '-s', '0']) #Parameterisation of the ligand
    #Remove intermediate files
    for fl in glob.glob('ANTECHAMBER_*'): #getting list of files with wildcard ANTECHAMBER_*
        os.remove(fl) #removes those files
    os.remove('ATOMTYPE.INF')
    os.remove('sqm.out')
    os.remove('sqm.in')
    os.remove('sqm.pdb')
    #Carry out parmchk
    subprocess.call([Amber_location+'parmchk2','-i', 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.mol2', '-f', 'mol2', '-o', 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.frcmod'])
    #Modify the tleap script
    savelines = []
    with open("Scripts/default_ligand.in") as fname: #opening file to read
        for line in fname:
            #replaces $$$ with ligand location and name
            newline = line.replace("$$$", 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name)
            savelines.append(newline)
    #Writing modified tleap script
    with open("Scripts/ligand.in", 'w') as outfile:
        outfile.writelines(savelines) #writing all lines that were saved
    #run tleap for the ligand
    subprocess.call([Amber_location+'tleap', '-f', 'Scripts/ligand.in'])

#%% Making a ligand-receptor complex (adding the ligand to the receptor membrane PDB file)
    #Add the ligand paramersiation pdb file to the receptor file
    #Adding disulphide CONECT bond lines of the CYX bonds
receptor_pdb_filenames = ['A1_FullyActive.pdb','A1_Inactive.pdb','A2a_Active.pdb','A2a_Inactive.pdb','B1_Active.pdb','B2_Active.pdb','B2_Inactive.pdb','CB1_Active.pdb','CB1_Inactive.pdb', 'CB2.pdb', 'CCR2.pdb', 'CRF-1.pdb', 'H1.pdb', 'M2_Inactive.pdb', 'M3.pdb', 'Mu_opoid_FullyActive.pdb', 'Mu_opoid_Inactive.pdb', 'Orexin2.pdb']
Receptor_names =['A1', 'A1', 'A2A', 'A2A', 'Beta1', 'Beta2', 'Beta2', 'CB1','CB1', 'CB2','CCR2', 'CRF1', 'H1','M2','M3', 'Mu_Opioid', 'Mu_Opioid', 'Orexin2'] #Setting up list of titles for plots
Receptors = [A1_Adenosine, A1_Adenosine, A2A_Adenosine, A2A_Adenosine, Beta1_Adrenoceptor, Beta2_Adrenoceptor,Beta2_Adrenoceptor, Cannabinoid_CB1, Cannabinoid_CB1, Cannabinoid_CB2, CCR2, CRF_1, H1_Histamine, M2_Muscarinic_acetylcholine, M3_Muscarinic_acetylcholine, Mu_Opioid_Receptor,Mu_Opioid_Receptor, Orexin_2] #setting up list of panda

for receptor_pdb_filename, Receptor, Receptor_DF in zip(receptor_pdb_filenames, Receptor_names, Receptors): #loops through all receptors
    savereceptorlines=[]
    CYXatomslist = []
    first_ter = True
    with open("Receptors/PM_EQ/AMBER/"+receptor_pdb_filename+"_EQ_Snapshot.pdb") as fname: #opening receptor-membrane file to read
            for line in fname: #loops through PDB lines
                #Adds CYX CONECT LINEs
                if line.startswith('END') == False: #saves all lines except the END line
                    if line[17:20] == 'CYX' and line[13:15] == 'SG': #checking if residue name equals CYX and atom name equals SG (Sulphur atom)
                        CYX = [line[6:11].strip(), line[30:39].strip(), line[39:47].strip(),line[47:54].strip()] #saves atom number, x, y and z coordinates in a list
                        CYXatomslist.append(CYX) #adds details to a list of a list
                    if line.startswith('TER') and first_ter == True: #checks for the first TER value
                        comb = combinations(CYXatomslist, 2) #does combination values of the all atoms in the list
                        for atom1, atom2 in comb:
                            if get_atom_distance(atom1, atom2) < 2.4: #if atom distance between two S atoms is under 2.4 angstroms must be disulphide bond
                                #making CONECT line record for PDB
                                if len(atom1[0]) == 3: #if atom number is in the hundreds,  need an extra space in CONECT line
                                    CONECTline = 'CONECT '+atom1[0]+' '+atom2[0] + '\n'
                                else:
                                    CONECTline = 'CONECT'+atom1[0]+' '+atom2[0] + '\n'
                                savereceptorlines.append(CONECTline)
                                first_ter = False
                    if line[17:20] == 'HIS' and line[13:16] == 'HD1': #avoiding writing HIS HD1 Atom as this breaks tleap
                        pass
                    else:
                        savereceptorlines.append(line) #saves the PDB line
    if path.exists("Receptors/PM_EQ/AMBER/"+receptor_pdb_filename+"_Charged_EQ_Snapshot.pdb"):
        #there is a charged PDB file
        with open("Receptors/PM_EQ/AMBER/"+receptor_pdb_filename+"_Charged_EQ_Snapshot.pdb") as fname:
            CYXatomslist = []
            first_ter = True
            savereceptorchargedlines = []
            for line in fname: #loops through PDB lines
                #Adds CYX CONECT LINEs
                if line.startswith('END') == False: #saves all lines except the END line
                    if line[17:20] == 'CYX' and line[13:15] == 'SG': #checking if residue name equals CYX and atom name equals SG (Sulphur atom)
                        CYX = [line[6:11].strip(), line[30:39].strip(), line[39:47].strip(),line[47:54].strip()] #saves atom number, x, y and z coordinates in a list
                        CYXatomslist.append(CYX) #adds details to a list of a list
                    if line.startswith('TER') and first_ter == True: #checks for the first TER value
                        comb = combinations(CYXatomslist, 2) #does combination values of the all atoms in the list
                        for atom1, atom2 in comb:
                            if get_atom_distance(atom1, atom2) < 2.4: #if atom distance between two S atoms is under 2.4 angstroms must be disulphide bond
                                #making CONECT line record for PDB
                                if len(atom1[0]) == 3: #if atom number is in the hundreds,  need an extra space in CONECT line
                                    CONECTline = 'CONECT '+atom1[0]+' '+atom2[0] + '\n'
                                else:
                                    CONECTline = 'CONECT'+atom1[0]+' '+atom2[0] + '\n'
                                savereceptorchargedlines.append(CONECTline)
                                first_ter = False
                    if line[17:20] == 'HIS' and line[13:16] == 'HD1': #avoiding writing HIS HD1 Atom as this breaks tleap
                        pass
                    else:
                        savereceptorchargedlines.append(line) #saves the PDB line
    for index, row in Receptor_DF.iterrows(): #looping through each ligand for in the dataframe
        ligand_name = row['Ligand'] #Takes the ligand name
        #Determining whether ligand is in the "Active" or "Inactive" or neither folder based on receptor name
        Activity = ''
        if 'Inactive' in receptor_pdb_filename:
            Activity = 'Inactive/'
        elif 'Active' in receptor_pdb_filename:
            Activity = 'Active/'
        #Checking charge of ligand
        charge = 0
        if ('+' in row['SMILES']) == True: #if there a + charge in the SMILES string
            charge = '+1' #charge of the ligand is +1
        ligand_savelines = []
        if path.exists('Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.pdb') == True:
            #Reads Ligand parameterisation pdb file and load all lines
            with open('Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'.pdb') as fname: #opening receptor-membrane file to read
                for line in fname:
                    ligand_savelines.append(line)
             #Add the ligand paramersiation pdb file to the receptor file in a new file *_com
            with open('Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'_com.pdb', 'w') as outfile:
                if charge == '+1':
                    outfile.writelines(savereceptorchargedlines) #writing all receptor lines (charged file)
                else:
                    outfile.writelines(savereceptorlines) #writing all receptor lines
                outfile.writelines(ligand_savelines) #writing all ligand lines

#%%Parameterisation of ligand protein complex
#determines whether ligand is charged, if so add cl- atom this will be done in tleap - need to fix this - make a charged version of each receptor, then select this version
Amber_location = '/Users/andrewpotterton/Documents/amber18/bin/'
Receptors = [A1_Adenosine, A2A_Adenosine, Beta1_Adrenoceptor, Beta2_Adrenoceptor, Cannabinoid_CB1, Cannabinoid_CB2, CCR2, CRF_1, H1_Histamine, M2_Muscarinic_acetylcholine, M3_Muscarinic_acetylcholine, Mu_Opioid_Receptor, Orexin_2] #setting up list of panda dataframes
Receptor_names =['A1', 'A2A', 'Beta1', 'Beta2', 'CB1','CB2','CCR2', 'CRF1', 'H1','M2','M3', 'Mu_Opioid', 'Orexin2']
for Receptor, Receptor_DF in zip(Receptor_names, Receptors): #loops through all receptors
    for index, row in Receptor_DF.iterrows():
        ligand_name = row['Ligand'] #Takes the ligand name
        if ('B-' in row['SMILES']) == False: #no boron atoms, that cannot be parameterised
            #Determines whether ligand is in the Active, Inactive folder or neither
            Activity = ''
            if path.exists('Conformers/'+Receptor+'/Parameterisation/Active/'+ligand_name+'.pdb') == True:
                Activity = 'Active/'
            elif path.exists('Conformers/'+Receptor+'/Parameterisation/Inactive/'+ligand_name+'.pdb') == True:
                Activity = 'Inactive/'
            #Modify the tleap script
            savelines = []
            with open("Scripts/default_PL.in") as fname: #opening file to read
                for line in fname:
                    #replaces $$$ with ligand location and name
                    newline = line.replace("$$$", 'Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name)
                    savelines.append(newline)
            #Writing modified tleap script
            with open("Scripts/PL.in", 'w') as outfile:
                outfile.writelines(savelines) #writing all lines that were saved
            #run tleap for the ligand
            subprocess.call([Amber_location+'tleap', '-f', 'Scripts/PL.in'])
#%% Copies over PDB and prmtop files for MD
#Receptors = [A1_Adenosine, A2A_Adenosine, Beta1_Adrenoceptor, Beta2_Adrenoceptor, Cannabinoid_CB1, Cannabinoid_CB2, CCR2, CRF_1, H1_Histamine, M2_Muscarinic_acetylcholine, M3_Muscarinic_acetylcholine, Mu_Opioid_Receptor, Orexin_2] #setting up list of panda dataframes
#Receptor_names =['A1', 'A2A', 'Beta1', 'Beta2', 'CB1','CB2','CCR2', 'CRF1', 'H1','M2','M3', 'Mu_Opioid', 'Orexin2']
Receptors = [FiveHT2A, FiveHT2B] #setting up list of panda dataframes
Receptor_names =['5HT2A', '5HT2B']
for Receptor, Receptor_DF in zip(Receptor_names, Receptors): #loops through all receptors
    for index, row in Receptor_DF.iterrows():
        ligand_name = row['Ligand'] #Takes the ligand name
        if ('B-' in row['SMILES']) == False: #no boron atoms, that cannot be parameterised
            #Determines whether ligand is in the Active, Inactive folder or neither
            Activity = ''
            if path.exists('Conformers/'+Receptor+'/Parameterisation/Active/'+ligand_name+'.pdb') == True:
                Activity = 'Active/'
            elif path.exists('Conformers/'+Receptor+'/Parameterisation/Inactive/'+ligand_name+'.pdb') == True:
                Activity = 'Inactive/'
            shutil.copy('Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'_complex.prmtop', 'MD/'+Receptor+'/Structures/'+ligand_name+'_complex.prmtop')
            shutil.copy('Conformers/'+Receptor+'/Parameterisation/'+Activity+ligand_name+'_complex.pdb', 'MD/'+Receptor+'/Structures/'+ligand_name+'_complex.pdb')

#%% Production MD scripts
Receptors = [FiveHT2A, FiveHT2B] #setting up list of panda dataframes
Receptor_names =['5HT2A', '5HT2B']
for Receptor, Receptor_DF in zip(Receptor_names, Receptors): #loops through all receptors
    for index, row in Receptor_DF.iterrows():
        ligand_name = row['Ligand'] #Takes the ligand name
        Activity = ''
        if ('B-' in row['SMILES']) == False: #no boron atoms, that cannot be parameterised
            #Getting PBC information
            #Determing what the PBC filename is called
            brackets = False #Sets default value of brackets to False
            if path.exists('Conformers/'+Receptor+'/Parameterisation/Active/'+ligand_name+'.pdb') == True:
                Activity ='Active'
            elif path.exists('Conformers/'+Receptor+'/Parameterisation/Inactive/'+ligand_name+'.pdb') == True:
                Activity ='Inactive'
            for filename in os.listdir('MD/'+Receptor+'/Structures/'):
                if filename.endswith('PBC.txt'):
                    if ('Inactive' in filename) == True:
                        if Activity == 'Inactive':
                            PBCfile = filename
                    elif ('Active' in filename) == True:
                        if Activity == 'Active':
                            PBCfile = filename
                    else:
                        if Activity == '':
                            PBCfile = filename
            with open('MD/'+Receptor+'/Structures/'+PBCfile) as fname:
                PBC = fname.readline()
                PBClist = PBC.split(' ') #splitting string by comma to get the individual values in a list
            #Generating input scripts for NAMD
            saveligand_name = ligand_name
            if ('[' in ligand_name) == True:
                saveligand_name = ligand_name.replace('[', '\[')
                saveligand_name = saveligand_name.replace(']', '\]')
                brackets = True
                saveligand_nobrackets = ligand_name.replace('[', '')
                saveligand_nobrackets = saveligand_nobrackets.replace(']', '')
            for i in range(1, 11): #reapting for 10 replicas
                for j in range(0, 6): #repeating for each input step
                    #Modifying the EQ scripts
                    savelines =[]
                    with open('Scripts/step7.'+str(j)+'_production.inp') as fname:
                        for line in fname:
                            newline = line.replace('$$$', saveligand_name)
                            newline = newline.replace('£X', PBClist[0])
                            newline = newline.replace('£Y', PBClist[1])
                            newline = newline.replace('£Z', PBClist[2])
                            savelines.append(newline)
                    with open('MD/'+Receptor+'/Rep'+str(i)+'/step7.'+str(j)+'_production_'+ligand_name+'.inp','w') as fname:
                        fname.writelines(savelines)
                #Modifying the EQ scripts

                savelines =[]
                with open('Scripts/step8.1_production.inp') as fname:
                    for line in fname:
                        newline = line.replace('$$$', saveligand_name)
                        savelines.append(newline)
                with open('MD/'+Receptor+'/Rep'+str(i)+'/step8.1_production_'+ligand_name+'.inp','w') as fname:
                    fname.writelines(savelines)
            #Generating Cartesius HPC scripts
            for i in range(1, 11): #reapting for 10 replicas
                savelines =[]
                if brackets == False:
                    with open('Scripts/default_myriadPL.sh') as fname:
                        for line in fname:
                            newline = line.replace('$$$', saveligand_name)
                            newline = newline.replace('???', Receptor)
                            newline = newline.replace('RepX', 'Rep'+str(i))
                            savelines.append(newline)
                else:
                    with open('Scripts/default_myriadPL_brackets.sh') as fname:
                        for line in fname:
                            newline = line.replace('$$$', saveligand_name)
                            newline = newline.replace('[[[', saveligand_nobrackets)
                            newline = newline.replace('???', Receptor)
                            newline = newline.replace('RepX', 'Rep'+str(i))
                            savelines.append(newline)
                with open('Scripts/HPC5/'+Receptor+'_'+ligand_name+'_Rep'+str(i)+'.sh','w') as fname:
                    fname.writelines(savelines)
#############################################
##    Copy Scripts and MD folders to HPC   ##
##    Run HPC scripts                      ##
#############################################
