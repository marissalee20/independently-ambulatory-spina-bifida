"""Copyright (c) 2022, Stanford Neuromuscular Biomechanics Laboratory
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import numpy as np
import os
import pandas as pd

from read_data import get_header_length, get_walkIDs

def generate_summary_waveforms(simulation_directory, print_subjects=True):
    """extract and standardize waveforms (ground reaction forces, kinematics,
    kinetics, joint reaction forces, muscle forces, and joint powers) for all
    subjects with simulation results. All standardized waveforms are resampled to
    100 timepoints.
    
    Args:
        simulation_directory (str): path to directory containing simulation results.
    
    Returns:
        summary_waveforms (pd DataFrame): dataframe containing waveforms for every
                                          subject trial.
    
    """
    subject_directories = get_subject_directories(simulation_directory)
    
    all_grfs = pd.DataFrame()
    all_kinematics = pd.DataFrame()
    all_kinetics = pd.DataFrame()
    all_jrfs = pd.DataFrame()
    all_muscle_forces = pd.DataFrame()
    all_powers = pd.DataFrame()
    for subject_directory in subject_directories:
        subject = subject_directory.split('_')[0]
        if print_subjects:
            print(subject)
            
        # get subject-specific waveform data
        subject_grfs = get_waveform_data(subject, 'grf')
        subject_kinematics = get_waveform_data(subject, 'kinematic')
        subject_kinetics = get_waveform_data(subject, 'kinetic')
        subject_jrfs = get_waveform_data(subject, 'jrf') # joint rxn
        subject_muscle_forces = get_waveform_data(subject, 'muscle force')
        subject_powers = get_waveform_data(subject, 'power')

        # append subject-specific waveform data
        all_grfs = all_grfs.append(subject_grfs)
        all_kinematics = all_kinematics.append(subject_kinematics)
        all_kinetics = all_kinetics.append(subject_kinetics)
        all_jrfs = all_jrfs.append(subject_jrfs)
        all_muscle_forces = all_muscle_forces.append(subject_muscle_forces)
        all_powers = all_powers.append(subject_powers)
        
    # combine into one dataframe
    keys = ['walk', 'side', 'subject', 'bodymass', 'height']
    all_kinetics = rename_kinetics(all_kinetics, omit_columns=keys)
    summary_waveforms = pd.merge(all_kinematics, all_kinetics, on=keys)
    summary_waveforms = pd.merge(summary_waveforms, all_jrfs, on=keys)
    summary_waveforms = pd.merge(summary_waveforms, all_grfs, on=keys)
    summary_waveforms = pd.merge(summary_waveforms, all_muscle_forces, on=keys)
    summary_waveforms = pd.merge(summary_waveforms, all_powers, on=keys)
    
    return summary_waveforms


def get_subject_directories(simulation_directory):
    """Get primary subject directories (not including EMG subject directories) held
    within the simulation directory.
    
    Args:
        simulation_directory (str): path to directory containing simulation results.
    
    Returns:
        subject_directories (list): list of paths of individual subject directories.
    
    """
    subject_directories = os.listdir(simulation_directory)
    subject_directories = [subject_directories[i]
                           for i in range(len(subject_directories))
                           if subject_directories[i].endswith('_SetupFiles')
                           and 'Template' not in subject_directories[i]
                           and not subject_directories[i].startswith('E')]
    return subject_directories


def get_waveform_data(subject, waveform_type):
    """collect 100-sample gait waveform data for a given subject.
    
    Args:
        subject (str): subject name.
        waveform_type (str): waveform type to extract. 'grf', 'kinematic', 'kinetic',
                             'jrf', 'muscle force', TODO.
    
    Returns:
        waveform_data (pd DataFrame): a dataframe containing waveform data per trial.
    """
    notes_filename = '../data/' + subject + '/' + subject + ' notes.xlsx'
    walkIDs = get_walkIDs(notes_filename)
    bodymass = get_bodymass(subject)
    height = get_height(subject)
    
    waveform_data = pd.DataFrame()
    for walkID in walkIDs:
        waveforms = get_raw_waveforms(subject, walkID, waveform_type)

        # trim and resample
        if waveform_type == 'grf':
            waveforms = trim_to_stride(waveforms, notes_filename, walkID)
            waveforms = resample(waveforms, 100)
        else:
            waveforms = resample(waveforms, 100)
        
        # stance limb waveforms only
        stance_side = get_stance_side(notes_filename, walkID)
        waveforms = keep_stance_waveform_only(waveforms, waveform_type, stance_side)
        
        waveforms = squeeze_df(waveforms)
        
        waveforms['walk'] = walkID
        waveforms['side'] = stance_side
        waveforms['bodymass'] = bodymass
        waveforms['height'] = height
        
        waveform_data = pd.concat((waveform_data, waveforms))
    waveform_data['subject'] = subject
    
    return waveform_data


def get_bodymass(subject):
    """extract bodymass from the demographics file.
    
    Args:
        subject (str): subject name.
    
    Returns:
        bodymass (float): subject bodymass in kg.
    """
    demographics = pd.read_excel('../data/demographics.xlsx', header=1)
    subject_row = demographics[demographics['subject']==subject]
    bodymass = subject_row['weight (kg)'].values[0]
    
    return bodymass


def get_height(subject):
    """extract height from the demographics file.
    
    Args:
        subject (str): subject name.
    
    Returns:
        height (float): subject height in cm.
    """
    demographics = pd.read_excel('../data/demographics.xlsx', header=1)
    subject_row = demographics[demographics['subject']==subject]
    bodymass = subject_row['height (cm)'].values[0]
    
    return bodymass


def get_raw_waveforms(subject, walkID, waveform_type):
    """extract raw gait waveform data for a given subject trial.
    
    Args:
        subject (str): subject name.
        walkID (str): walk identifier.
        waveform_type (str): waveform type to extract. 'grf', 'kinematic', 'kinetic',
                             or 'jrf'.
    
    Returns:
        waveforms (pd DataFrame): a dataframe containing raw waveform data for the trial.
    """
    subject_directory = '../simulation/' + subject + '_SetupFiles/'
    if waveform_type=='grf':
        filename = (subject_directory + '../../data/' + subject + '/' + subject
                    + '_' + walkID + '_grf_filtered.mot')
        grf_column_starts = ['1_ground_force_v', '2_ground_force_v', 
                             '3_ground_force_v']
        keywords = ([grf + 'x' for grf in grf_column_starts]
                    + [grf + 'y' for grf in grf_column_starts]
                    + [grf + 'z' for grf in grf_column_starts])
    elif waveform_type=='kinematic':
        filename = (subject_directory + '/RRA/results_rra_' + walkID
                    + '/rra_' + walkID + '_Kinematics_q.sto')
        keywords = ['ankle', 'knee', 'hip']
    elif waveform_type=='kinetic':
        filename = (subject_directory + '/RRA/results_rra_' + walkID
                    + '/rra_' + walkID + '_Actuation_force.sto')
        keywords=['ankle', 'knee', 'hip']
    elif waveform_type=='jrf':
        filename = (subject_directory + '/SO_Custom/' + walkID
                    + '_results_JointRxn_ReactionLoads.sto')
        keywords=['fx', 'fy', 'fz']
    elif waveform_type=='muscle force':
        filename = (subject_directory + '/SO_Custom/' + walkID
                    + '_results_forces.sto')
        keywords = ['gasmed_r', 'gasmed_l',
                    'gaslat_r', 'gaslat_l',
                    'soleus_r', 'soleus_l',
                    'fdl_r', 'fdl_l',
                    'fhl_r', 'fhl_l',
                    'perbrev_r', 'perbrev_l',
                    'perlong_r', 'perlong_l',
                    'tibpost_r', 'tibpost_l',

                    'recfem_r', 'recfem_l',
                    'vaslat_r', 'vaslat_l',
                    'vasint_r', 'vasint_l',
                    'vasmed_r', 'vasmed_l',
                   ]
    elif waveform_type=='power':
        pass
    else:
        raise Exception('Unknown waveform type. Expected "grf", "kinematic", \
                        "kinetic", "jrf", "muscle force", or "power" but received:',
                        waveform_type)
    
    if waveform_type=='power':
        filename = (subject_directory + '/SO_Custom/' + walkID + '_results_power.csv')
        waveforms = pd.read_csv(filename, index_col=False)
    else:
        header_length = get_header_length(filename)
        waveforms = pd.read_csv(filename, sep='\t', header=header_length, index_col=False)
        time = waveforms['time'] # save for later
    
        # keep only kw waveforms
        waveform_names = [waveform for waveform in waveforms
                          if any(kw in waveform for kw in keywords)]
        waveforms = waveforms[waveform_names]        
        waveforms['time'] = time
    return waveforms


def trim_to_stride(waveforms, notes_filename, walkID):
    """trim waveforms to only duration of stride of interest.
    
    Args:
        waveforms (pd DataFrame): dataframe containing waveforms.
        notes_filename (str): filepath to notes file containing walk identifiers.
        walkID (str): walk identifier.
    
    Returns:
        waveforms_trimmed (pd DataFrame): dataframe containing trimmed waveforms.
    """
    start_time, end_time = get_stride_endpoints(walkID, notes_filename)
    start_idx = np.where(waveforms['time'].to_numpy() > start_time)[0][0]
    stop_idx = np.where(waveforms['time'].to_numpy() > end_time)[0][0] - 1
    waveforms_trimmed = waveforms.iloc[start_idx:stop_idx]
    
    return waveforms_trimmed


def resample(waveforms, n_samples=100):
    """resample waveforms in df to have a set number of points.
    
    Args:
        waveforms (pd DataFrame): dataframe containing waveforms.
        n_samples (int): number of samples resampled waveforms should have.
        
    Returns:
        waveforms_resampled (pd DataFrame): dataframe containing resampled waveforms.
    """
    waveforms_resampled = pd.DataFrame()
    for waveform in list(waveforms):
        start = 0
        stop = len(waveforms[waveform])-1
        waveforms_resampled[waveform] = np.interp(np.linspace(start, stop,
                                                              num=n_samples),
                                           np.arange(len(waveforms)),
                                                  waveforms[waveform])
    return waveforms_resampled


def get_stance_side(notes_filename, walkID):
    """extract laterality of stance limb (limb that strikes force plate 2) from the
    subject's notes file.
    
    Args:
        notes_filename (str): filepath to notes file containing walk identifiers.
        walkID (str): walk identifier.
        
    Returns:
        stance_side (str): laterality of stance limb, 'l' or 'r'.
    """
    notes = pd.read_excel(notes_filename)
    walk_row = notes[notes['trial']==walkID.capitalize()]
    stance_side = walk_row['force plate order'].values[0][3].lower()
    return stance_side


def keep_stance_waveform_only(waveforms, waveform_type, stance_side):
    """filter and rename columns to only include waveforms of the stance limb.
    
    Args:
        waveforms (pd DataFrame): dataframe containing waveforms.
        waveform_type (str): waveform type to extract. 'grf', 'kinematic', 'kinetic',
                             or 'jrf'.
        stance_side (str): laterality of stance limb, 'l' or 'r'.
    
    Returns:
        stance_waveforms (pd DataFrame): dataframe containing stance waveforms only.
    """
    if waveform_type=='grf':
        stance_waveforms = pd.DataFrame()
    
        stance_waveforms['grfx'] = waveforms['2_ground_force_vx'].copy()
        stance_waveforms['grfy'] = waveforms['2_ground_force_vy'].copy()
        stance_waveforms['grfz'] = waveforms['2_ground_force_vz'].copy()
        if stance_side=='r':
            stance_waveforms['grfz'] *= -1 # point medially

    else:
        key1 = '_' + stance_side + '_'
        key2 = '_' + stance_side
        key3 = '_' + stance_side + '/'
        columns_to_keep = [col for col in waveforms 
                           if key1 in col or col[-2:]==key2 or key3 in col]

        # filter
        stance_waveforms = waveforms[columns_to_keep]

        # rename to omit 'r'/'l'
        name_map = {}
        for col in stance_waveforms:
            name_map[col] = col.replace(key1, '_')
            if name_map[col][-2:]==key2:
                name_map[col] = name_map[col][:-2]
            if '/' in col:
                name_map[col] = col.replace(key3, '/')

        stance_waveforms = stance_waveforms.rename(columns=name_map)
    
    return stance_waveforms


def squeeze_df(df):
    """squeeze arrays into single row of dataframe.
    
    Args:
        df (pd DataFrame): dataframe containing waveforms.
    
    Returns:
        df_squeezed (pd DataFrame): reformatted dataframe.
    """
    headers = list(df)
    data = df.to_numpy()
    df_squeezed = pd.DataFrame()
    for col_idx in range(data.shape[1]):
        array = data[:,col_idx]
        df_squeezed[headers[col_idx]] = [list(array)]
    return df_squeezed


def get_stride_endpoints(walkID, notes_filename):
    """extract stride start and end times from the subject's note file.
    
    Args:
        walkID (str): walk identifier.
        notes_filename (str): filepath to notes file containing walk identifiers.
        
    Returns:
        start_time (float): stride start time.
        end_time (float): stride end time.
    """
    notes = pd.read_excel(notes_filename)
    walk_row = notes[notes['trial']==walkID.capitalize()]
    start_time = walk_row['start'].values[0]
    end_time = walk_row['end'].values[0]
    return start_time, end_time


def rename_kinetics(df, omit_columns=[]):
    """rename data columns of an input dataframe to have a '_moment' suffix. This is
    important when merging information from a kinematics and kinetics output.
    
    Args:
        df (pd DataFrame): dataframe with columns to rename.
        omit_columns (list): list of strings of column names to omit from name
                             adjustments (typically keys on which you will later
                             merge dataframes).
    
    Returns:
        df_renamed (pd DataFrame): original dataframe with all columns that do not
                                   match 'walk', 'side', 'subject', or 'bodymass'
                                   renamed to be followed by '_moment'
    """
    name_map = {}
    for column in df:
        if column not in omit_columns:
            name_map[column] = column + '_moment'
    
    df_renamed = df.rename(columns=name_map)
    return df_renamed


def average_subject_trials(waveforms):
    """Combine trials from a single subject into an average trial.
    
    Args:
        waveforms (pd DataFrame): dataframe containing waveforms.
        
    Returns:
        avg_waveforms (pd DataFrame): dataframe containing averaged waveforms.
    """
    avg_waveforms = pd.DataFrame()
    for subject in set(waveforms['subject']):
        subject_data = waveforms[waveforms['subject']==subject]
        bodymass = list(set(list(subject_data['bodymass'])))
        height = list(set(list(subject_data['height'])))
        side = list(set(list(subject_data['side'])))

        subject_data = subject_data.drop(columns=['side', 'subject', 'bodymass', 'height'])
        headers = list(subject_data)
            
        for col in subject_data:
            subject_data[col] = subject_data[col].apply(np.array)
        subject_mean = np.expand_dims(subject_data.to_numpy(), axis=1).mean(axis=0)
        subject_mean = pd.DataFrame(subject_mean, columns=headers)
        subject_mean['subject'] = subject
        subject_mean['bodymass'] = bodymass[0]
        subject_mean['height'] = height[0]
        subject_mean['side'] = side[0]
        avg_waveforms = avg_waveforms.append(subject_mean, ignore_index=True)
    
    return avg_waveforms
