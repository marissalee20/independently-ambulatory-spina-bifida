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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal

from matplotlib.lines import Line2D

from read_data import get_analyzed_side_pf_strength

def plot_emg_comparisons(emg_demographics):
    """plot comparisons between static optimization muscle activation estimates and
    EMG signals.
    
    Args:
        emg_demographics (pd DataFrame): demographic information for subjects with
                                         EMG.
    
    Returns:
        none.
        
    """
    ## definitions ##
    PURPLE = '#882255'
    BLUE = '#0077BB'
    emg_dict = {'gasmed': 'gastrocnemius\nmedialis',
                'tibant': 'tibialis\nanterior',
                'semiten': 'semi-\ntendinosus',
                'recfem': 'rectus\nfemoris',
                'vaslat': 'vastus\nlateralis'
               }
    muscles = list(emg_dict.keys())
    n_muscles = len(muscles)
    # conversion from numbers to scores
    strength_dict = {5: '5',
                     4: '4',
                     3.67: '4-',
                     3.33: '3+',
                     3: '3',
                     2.67: '3-',
                     2.33: '2+',
                     2: '2',
                     1.67: '2-',
                     1.33: '1+',
                     1: '1',
                     0.67: '1-',
                     0: '0'
                    }
    
    ylims = {'gasmed': 0.6,
             'tibant': 0.8,
             'semiten': 0.08,
             'recfem': 0.6,
             'vaslat': 0.6
            }
    ################
    
    subjectIDs = emg_demographics['subject']
    pf_strengths, analyzed_sides = get_analyzed_side_pf_strength(emg_demographics)
    subjectIDs, pf_strengths = order_lists(subjectIDs, pf_strengths, reverse=True)

    n_subjects = len(subjectIDs)
    
    f, axes = plt.subplots(n_subjects, n_muscles, dpi=300)
    for i in range(n_subjects):
        subjectID = subjectIDs[i]
        walk_side = analyzed_sides[i]
        strength = pf_strengths[i]
        
        for j in range(n_muscles):
            muscle = muscles[j]

            plt.subplot(n_subjects, n_muscles, 5*i+j+1)
            ax = plt.gca()
            if i == 0:
                plt.title(emg_dict[muscle], fontsize=18)

            # average results over walks
            f_notes = '../data/' + subjectID + '/' + subjectID + ' notes.xlsx'
            walks, sides = get_walks_and_sides(f_notes)
            n_walks = len(walks)
            
            activations = np.empty((n_walks, 100))
            emgs = np.zeros((n_walks, 100))
            for k in range(n_walks):
                walk_id = walks[k]
                walk_side = sides[k]

                activations[k,:] = get_activation(subjectID, walk_id, walk_side,
                                                  muscle) - 0.01 # subtract offset
                emg = get_emg(subjectID, walk_id, walk_side, muscle)
                if np.any(emg):
                    emgs[k,:] = emg

            if np.any(emgs): # if there's emg data, plot
                mean_activation = np.mean(activations, axis=0)
                std_activation = np.std(activations, axis=0)
                mean_emg = np.mean(emgs, axis=0)

                # scale EMG to match SO peak
                emgs = emgs/np.max(mean_emg)*np.max(mean_activation)
                mean_emg = np.mean(emgs, axis=0)
                std_emg = np.std(emgs, axis=0)

                # plot EMG
                plt.plot(mean_emg, color=BLUE, linewidth=2)
                plt.fill_between(np.linspace(0, 100, 100), mean_emg - std_emg,
                                 mean_emg + std_emg, color=BLUE, alpha=0.25,
                                 linewidth=0)

                # plot static optimization
                plt.plot(mean_activation, color=PURPLE, linewidth=2)
                plt.fill_between(np.linspace(0, 100, 100),
                                 mean_activation - std_activation,
                                 mean_activation + std_activation,
                                 color=PURPLE, alpha=0.25, linewidth=0)

                if 'gas'in muscle or 'rec' in muscle or 'vas' in muscle:
                    yticks = [0, 0.2, 0.4, 0.6]
                elif 'tib' in muscle:
                    yticks = [0, 0.2, 0.4, 0.6, 0.8]
                elif 'semi' in muscle:
                    yticks = [0, 0.02, 0.04, 0.06, 0.08]
                plt.yticks(yticks, yticks)

            else: # remove plot if no EMG
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                plt.tick_params(bottom=False, left=False)
                plt.xticks([])
                plt.yticks([])
                plt.text(50, 0.3, 'data\nnot\ncollected', fontsize=12, ha='center', va='center')

            if j == 0:
                plt.ylabel('strength\n' + strength_dict[strength] + '/5', rotation=0, verticalalignment='center', fontsize=18)
                ax.yaxis.set_label_coords(-0.6, 0.5)
            if i == n_subjects-1 and j == 2:
                plt.xlabel('% gait cycle', fontsize=18)
                ax.xaxis.set_label_coords(0.5, -0.5)
            if i != n_subjects-1 and (i != 5 or j != 0):
                ax.xaxis.set_ticklabels([])
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.tick_params(top=False, right=False)
            plt.ylim([0, ylims[muscle]])
            plt.xlim([0, 100])

            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

    f.set_figheight(15)
    f.set_figwidth(15)
    purple_line = Line2D([0,1],[0,1], linestyle='-', linewidth=2, color=PURPLE)
    blue_line = Line2D([0,1],[0,1], linestyle='-', linewidth=2, color=BLUE)
    plt.legend([purple_line, blue_line], ['static optimization', 'EMG'], loc='center right', fontsize=18,
               bbox_to_anchor=(-1.9, -1.7, 1, 1), frameon=False, ncol=2)
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.savefig('../postprocessing/figures/emg.jpeg', bbox_inches='tight')
    plt.show()
    return


def order_lists(subjectIDs, by, reverse=False):
    """re-order two lists together.
    
    Args:
        subjectIDs (list): secondary list to order.
        by (list): list by which to order two lists.
        reverse (bool): True for reverse ordering.
        
    Returns:
        subjectIDs_ordered (list): subjectIDs, ordered by 'by' list.
        by_ordered (list): 'by' list, ordered.
    
    """
    zipped_lists = zip(by, subjectIDs)
    sorted_pairs = sorted(zipped_lists, reverse=reverse)
    tuples = zip(*sorted_pairs)
    by_ordered, subjectIDs_ordered = [list(tuple) for tuple in tuples]
    
    return subjectIDs_ordered, by_ordered


def get_walks_and_sides(notes_filename):
    """from a subject's notes, extract walk IDs and stance laterality.
    
    Args:
        notes_filename (str): path to subject notes.xlsx file.
    
    Returns:
        walks (list): list of walk identifiers.
        sides (list): list of stance laterality.
        
    """
    notes = pd.read_excel(notes_filename)
    walks = list(notes['trial'])
    
    # find second force plate sides
    sides = [patt[patt.index('2')-1] for patt in list(notes['force plate order'])]
    
    return walks, sides


def get_activation(subject_id, walk_id, walk_side, muscle):
    """get activation waveform of a particular muscle from a particular subject's particular walk.
    
    Args:
        subject_id (str): subject ID.
        walk_id (str): walk ID.
        walk_side (str): 'L' or 'R', stance side associated with walk.
        muscle (str): name of muscle.
    
    Returns:
        activation (np arr): array containing muscle activation waveform.
    
    """
    f_so = ('../simulation/' + subject_id + '_SetupFiles/SO_Custom/' + walk_id + 
                    '_results_states.sto') # static optimization states filename
    states = pd.read_csv(f_so, sep='\t', header=6)
    activation = np.array(states['/forceset/' + muscle + '_' + walk_side.lower()
                                 + '/activation'])
    return activation


def get_emg(subject_id, walk_id, walk_side, muscle):
    """get filtered emg waveform of a particular muscle from a particular subject's
    particular walk.
    
    Args:
        subject_id (str): subject ID.
        walk_id (str): walk ID.
        walk_side (str): 'L' or 'R', stance side associated with walk.
        muscle (str): name of muscle.
    
    Returns:
        emg (np arr): array containing filtered emg waveform.
    
    """
    f_emg = ('../data/' + subject_id + '/' + subject_id + '_' + walk_id +
             '_emg.xlsx')
    emgs = pd.read_excel(f_emg)
    
    muscle_name = muscle + '_' + walk_side.lower()
    if muscle_name in list(emgs):
        emg = np.array(emgs[muscle_name])
        emg = filter_emg(emg, fs=2400)
        emg = np.interp(np.linspace(0,100,100), np.linspace(0,100,len(emg)), emg)
    else:
        emg = np.zeros((1,100))
    
    return emg


def filter_emg(signal_raw, fs):
    """filter and rectify emg signal
    
    Args:
        signal_raw (np arr): array containing raw emg waveform.
        fs (int): raw signal sampling frequency (Hz).
    
    Returns:
        signal_filtered (np arr): arr containing filtered and rectified emg waveform.
        
    """
    # 4th-order bandpass filter. Note bandpass and filtfilt each double order.
    fc_low = float(50)
    fc_high = float(500)
    b, a = signal.butter(1, [fc_low/(fs/2), fc_high/(fs/2)], btype='bandpass')
    signal_filtered = signal.filtfilt(b, a, signal_raw)
    
    # rectify
    signal_filtered = np.abs(signal_filtered)
    
    # 4th-order lowpass filter
    fc = 7.5
    b, a = signal.butter(2, fc/(fs/2))
    signal_filtered = signal.filtfilt(b, a, signal_filtered)
    
    return signal_filtered
