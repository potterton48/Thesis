#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for Kinetic_data.py and KineticDataCorrelations
Created on Wed Aug 14 13:57:28 2019

@author: andrewpotterton
"""

import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Descriptors
from scipy import stats
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection

###Functions###
def get_noHalogens(SMILES):
    """Returns number of halogen atoms in a ligand.
    Required Rdkit Chem library to be imported
    Parameters
    ----------
    SMILES = SMILES string of ligand to be tested.
    """
    m = rdkit.Chem.MolFromSmiles(SMILES)
    atoms = m.GetAtoms()
    no_halogens = 0
    for atom in atoms:
        atom_name = atom.GetSymbol()
        if atom_name == 'Cl' or atom_name == 'Br' or atom_name == 'F' or atom_name == 'I':
            no_halogens = no_halogens + 1
    return no_halogens

def get_no_atoms(SMILES):
    """Returns number of atoms in a ligand.
    Required Rdkit Chem library to be imported
    Parameters
    ----------
    SMILES = SMILES string of ligand to be tested.
    """
    m = rdkit.Chem.MolFromSmiles(SMILES)
    atoms = m.GetAtoms()
    atom_count = 0
    for atom in atoms:
        atom_count=atom_count+1
    return atom_count

def remove_duplicates(Dataframe):
    """Removes duplicate values from a dataframe.
    Parameters
    ----------
    Dataframe = Kinetic data dataframe
    """
    bool_duplicates = Dataframe.duplicated(subset='SMILES', keep=False)
    duplicates = Dataframe[bool_duplicates]
    duplicate_names = duplicates['SMILES'].unique()
    column_names_to_be_averaged = ['Temp (K)', 'Kon (1/(nM*min)', 'Error', 'Koff (1/min)','Error.1', 'EA', 'Corrected Koff 298.15', 'RT (min)', 'Error.2', 'Corrected RT (298.15)', 'Error.3', 'Kd (nM)', 'Error.4', 'Kinetic Kd (nM)', 'Error.5', 'Ki (nM)', 'Error.6']
    for name in duplicate_names:
        ligand_dup = duplicates.loc[duplicates['SMILES'] == name]
        index_value = ligand_dup.index.tolist()[0]
        for column in column_names_to_be_averaged:
            average_value=ligand_dup[column].mean()
            Dataframe.loc[[index_value],[column]]= average_value
        #columns to be recalculated based on new average value
        new_Free_energy = 1.9872036*10**(-3)*Dataframe.iloc[index_value]['Temp (K)']*np.log(Dataframe.iloc[index_value]['Kinetic Kd (nM)']*10**(-9))
        Dataframe.loc[[index_value],['Free Energy (kcal/mol)']]= new_Free_energy
        new_Free_energy_error=(1/2.303)*Dataframe.iloc[index_value]['Free Energy (kcal/mol)']*(Dataframe.iloc[index_value]['Error.5']/Dataframe.iloc[index_value]['Kinetic Kd (nM)'])
        Dataframe.loc[[index_value],['Error (kcal/mol)']]= new_Free_energy_error
        new_pKi=-1*np.log10(Dataframe.iloc[index_value]['Ki (nM)']*10**(-9))
        Dataframe.loc[[index_value],['pKi']]= new_pKi
    Dataframe=Dataframe.drop_duplicates(subset=['SMILES'], keep='first')
    Dataframe.reset_index(inplace=True, drop=True)
    return Dataframe

def get_descriptors(Kinetic_dataframe):
    """Appends chemical properties to the kinetic dataframe using the SMILES data.
    And removes duplicates

    Parameters
    ----------
    Kinetic_dataframe = Kinetic data dataframe
    """
    descriptors_list=[]
    for m in Kinetic_dataframe['SMILES']:
        molecule = rdkit.Chem.MolFromSmiles(m)
        MW = rdkit.Chem.Descriptors.ExactMolWt(molecule)
        HA_MW = rdkit.Chem.Descriptors.HeavyAtomMolWt(molecule)
        logP = rdkit.Chem.Descriptors.MolLogP(molecule)
        L_HBA = rdkit.Chem.rdMolDescriptors.CalcNumLipinskiHBA(molecule)
        L_HBD = rdkit.Chem.rdMolDescriptors.CalcNumLipinskiHBD(molecule)
        HB_atoms = L_HBA+L_HBD
        no_rings = rdkit.Chem.rdMolDescriptors.CalcNumRings(molecule)
        no_arom_rings = rdkit.Chem.rdMolDescriptors.CalcNumAromaticRings(molecule)
        TPSA = rdkit.Chem.rdMolDescriptors.CalcTPSA(molecule)
        no_halogens = get_noHalogens(m)
        no_atoms = get_no_atoms(m)
        HB_atoms_pa = HB_atoms/no_atoms
        No_halogens_pa = no_halogens/no_atoms
        no_rot_bonds = rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds(molecule)
        descriptors = [MW,HA_MW,logP,L_HBA,L_HBD,HB_atoms, HB_atoms_pa, no_rings,no_arom_rings,TPSA, no_halogens, No_halogens_pa, no_atoms, no_rot_bonds]
        descriptors_list.append(descriptors)
    Descriptions = pd.DataFrame(descriptors_list, columns=['MW', 'HA_MW', 'logP', 'L_HBA','L_HBD', 'HB_atoms', 'HB_atoms_pa','no_rings', 'no_arom_rings', 'TPSA', 'No_halogens', 'No_halogens_pa','No_atoms', 'No_rot_bonds'] )
    Kinetic_dataframe= Kinetic_dataframe.join(Descriptions)
    Kinetic_dataframe= remove_duplicates(Kinetic_dataframe)
    return Kinetic_dataframe

def get_correlations(dfcolumn1, dfcolumn2):
    """Provide two dataframe columns for a correlation value to be worked out for.
    This Function checks and removes Nan values from the correlation
    Parameters
    ----------
    dfcolumn1 = first column that will be correlated e.g. A1_Adenosine['RT']
    dfcolumn2 = second column, can be from separate dataframes but must have the same number of values in each column
    """
    bad = ~np.logical_or(np.isnan(dfcolumn1), np.isnan(dfcolumn2)) #checks for nan values
    x = np.compress(bad, dfcolumn1) #removes nan values from both columns
    y = (np.compress(bad, dfcolumn2))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return (r_value, slope, intercept)

def get_properties_correlation_matrix(Kinetic_dataframe, column_name):
    "column is the name of the column in the kinetic dataframe that correlations will be provided against"
    Properties = pd.concat([Kinetic_dataframe[column_name],Kinetic_dataframe.iloc[:,28:37]], axis=1)
    Properties['HB_Atoms'] = Properties['L_HBA'] + Properties['L_HBD']
    Properties = Properties.drop(['HA_MW', 'L_HBA','L_HBD'], axis=1)
    r_values = Properties.corr()[column_name]
    r_values = r_values.drop(column_name)
    r_values= r_values.pow(2)
    return r_values
def get_properties_correlation_matrix_log(Kinetic_dataframe, column_name):
    "column is the name of the column in the kinetic dataframe that correlations will be provided against"
    Properties = pd.concat([np.log(Kinetic_dataframe[column_name]),Kinetic_dataframe.iloc[:,28:37]], axis=1)
    Properties['HB_Atoms'] = Properties['L_HBA'] + Properties['L_HBD']
    Properties = Properties.drop(['HA_MW', 'L_HBA','L_HBD'], axis=1)
    r_values = Properties.corr()[column_name]
    r_values = r_values.drop(column_name)
    r_values= r_values.pow(2)
    return r_values

def radar_factory(num_vars, frame='circle'):
    """Create a radar chart with `num_vars` axes.

    This function creates a RadarAxes projection and registers it.

    Parameters
    ----------
    num_vars : int
        Number of variables for radar chart.
    frame : {'circle' | 'polygon'}
        Shape of frame surrounding axes.

    """
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    def draw_poly_patch(self):
        # rotate theta such that the first axis is at the top
        verts = unit_poly_verts(theta + np.pi / 2)
        return plt.Polygon(verts, closed=True, edgecolor='k')

    def draw_circle_patch(self):
        # unit circle centered on (0.5, 0.5)
        return plt.Circle((0.5, 0.5), 0.5)

    patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
    if frame not in patch_dict:
        raise ValueError('unknown value for `frame`: %s' % frame)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        # define draw_frame method
        draw_patch = patch_dict[frame]

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.concatenate((x, [x[0]]))
                y = np.concatenate((y, [y[0]]))
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            return self.draw_patch()

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            # The following is a hack to get the spines (i.e. the axes frame)
            # to draw correctly for a polygon frame.

            # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
            spine_type = 'circle'
            verts = unit_poly_verts(theta + np.pi / 2)
            # close off polygon by repeating first vertex
            verts.append(verts[0])
            path = Path(verts)

            spine = Spine(self, spine_type, path)
            spine.set_transform(self.transAxes)
            return {'polar': spine}

    register_projection(RadarAxes)
    return theta


def unit_poly_verts(theta):
    """Return vertices of polygon for subplot axes.

    This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
    """
    x0, y0, r = [0.5] * 3
    verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
    return verts

def make_radar_plot(df_data, spoke_labels, title):
    "Makes radar plot of the correlation vaules provided in the dataframe (df). Title will be the title of the plot"
    theta = radar_factory(len(spoke_labels), frame='polygon')
    fig, ax = plt.subplots(subplot_kw=dict(projection='radar'))
    ax.plot(theta, df_data)
    ax.fill(theta, df_data, alpha=0.25)
    ax.set_frame_on(False)
    plt.xticks((theta), spoke_labels, fontsize='large')
    ax.set_rmax(1)
    ax.set_rgrids([0.25, 0.5, 0.75])
    ax.set_rlabel_position(-22.5)
    plt.tight_layout()
    plt.savefig('Figures/'+title+'.png', dpi=1000)

def make_radar_plot_overlay(df1, df2, df3, spoke_labels, title):
    "Makes overlay radar plot of the correlation vaules provided in the dataframes (df1-3). Title will be the title of the plot"
    theta = radar_factory(len(spoke_labels), frame='polygon')
    fig, ax = plt.subplots(subplot_kw=dict(projection='radar'))
    ax.plot(theta, df1, label='CCR2')
    ax.plot(theta, df2, label='M3')
    ax.plot(theta, df3, label='mGlu')
    plt.legend(bbox_to_anchor=(1.5,1.0), title='Receptor')
    ax.fill(theta, df1, alpha=0.15, label=None)
    ax.fill(theta, df2, alpha=0.15, label=None)
    ax.fill(theta, df3, alpha=0.15, label=None)
    ax.set_frame_on(False)
    plt.xticks((theta), spoke_labels, fontsize='large')
    ax.set_rmax(1)
    ax.set_rgrids([0.25, 0.5, 0.75])
    ax.set_rlabel_position(-80)
    plt.tight_layout()
    plt.savefig('Figures/'+title+'.png', dpi=1000)


def generateconformations(m, name, receptor_folder_name):
    """
    Generates conformers for a ligand and saves the lowest confomer as a single sdf.
    Also save all ligand conformers in a single sdf.
    The number of conformers generated depends on the number of rotatable bonds in the ligand.
    Requires rdkit
    Parameters
    ----------
    m = molecule - needs to be a rdkit molecule (e.g. from a smiles string you would need to run: m = rdkit.Chem.MolFromSmiles(row['SMILES']) )
    name = the name of the ligand as a string (which will be used to name the sdf file)
    receptor_folder_name = name of receptor as a string (A folder needs to be created with this name in the Confomer folder)
    """
    #adding Hydrogens to molecule
    m = Chem.AddHs(m)
    #determining the number of rotable bonds, then setting the no of conformers to generate based on that value
    no_rot_bonds = rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds(m)
    if no_rot_bonds <= 7:
        n = 50
    elif no_rot_bonds <= 12:
        n = 200
    else:
        n = 300
    #Generating conformations
    ids=AllChem.EmbedMultipleConfs(m, numConfs=n, pruneRmsThresh=0.5)
    #Energy Minimising each conformers and scoring each confomer
    energy_conf_list = []
    for idv in ids:
        ff = AllChem.UFFGetMoleculeForceField(m, idv) #applies forcefield to the conformer
        ff.Initialize()
        ff.Minimize() #energy minimises
        energy = float(ff.CalcEnergy()) #calculates energy of the minimised conformer
        energy_conf = (energy, idv) #saves this energy score, and the conformer id
        energy_conf_list.append(energy_conf)
    energy_conf_list.sort() #sorts energy from the lowest energy first (the best score)
    #Writing lowest energy conformers
    lowest_energy=energy_conf_list[0]
    m.SetProp('_Name', name)
    m.SetProp('ConfEnergies', str(lowest_energy[0])+ " Kcal/mol")
    writer = Chem.SDWriter('Conformers/'+receptor_folder_name+'/Lowest_energy/'+name+'.sdf')
    writer.write(m, confId=lowest_energy[1])
    writer.close()
    #Removing conformers that are similar (RMSD lowest than 0.35 angstroms after energy minimisation)
    kept_ligands = []
    kept_ligands.append(energy_conf_list[0])
    for conf in energy_conf_list:
        rms_lowest = 100.0
        for lig in kept_ligands:
            rms = AllChem.GetConformerRMS(m, lig[1], conf[1])
            if rms < rms_lowest:
                rms_lowest = rms
        if rms_lowest > 0.5:
            kept_ligands.append(conf)
    #Writing all conformers
    writer = Chem.SDWriter('Conformers/'+receptor_folder_name+'/'+name+'.sdf')
    m.SetProp('_Name', name)
    for val in kept_ligands:
        m.SetProp('ConfId', str(val[1]))
        m.SetProp('ConfEnergies', str(val[0])+ " Kcal/mol")
        writer.write(m, confId=val[1])
    writer.close()
    return lowest_energy

def get_atom_distance(atom1, atom2):
    #returns distance between atom1 and atom2 (atom1 and atom2 are lists in the format: [atomno, x, y, z])
    from math import sqrt
    distance = sqrt( (float(atom1[1]) - float(atom2[1]))**2 + (float(atom1[2]) - float(atom2[2]))**2 + (float(atom1[3]) - float(atom2[3]))**2)
    return distance
