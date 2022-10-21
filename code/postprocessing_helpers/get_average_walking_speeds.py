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
import pandas as pd

from read_data import get_header_length, get_walkIDs

def get_average_walking_speeds(subjects):
    """calculate and extract average walking speeds for given subjects.
    
    Args:
        subjects (list): list of subject identifiers.
        
    Returns:
        average_walking_speeds (np array): average walking speed associated with each
                                           subject.
    """
    n_subjects = len(subjects)
    average_walking_speeds = np.zeros(n_subjects)
    for i in range(n_subjects):
        subject = subjects[i]
        average_walking_speeds[i] = get_average_walking_speed(subject)
        
    return average_walking_speeds


def get_average_walking_speed(subject):
    """calculate average walking speed across subject trials.
    
    Args:
        subject (str): subject identifier.
    
    Returns:
        average_walking_speed (float): in m/s.
    
    """
    notes_filename = '../data/' + subject + '/' + subject + ' notes.xlsx'
    walkIDs = get_walkIDs(notes_filename)
    n_walkIDs = len(walkIDs)
    
    walking_speeds = np.zeros(n_walkIDs)
    for i in range(n_walkIDs):
        walkID = walkIDs[i]
#         position_file = ('../simulation/' + subject
#                          + '_SetupFiles/Analysis/results_analysis_' + walkID
#                          + '/analysis_' + walkID + '_BodyKinematics_pos_global.sto') # TODO if no RRA
        position_file = ('../simulation/' + subject
                         + '_SetupFiles/Analysis/PostRRA/results_analysis_' + walkID
                         + '/analysis_' + walkID + '_BodyKinematics_pos_global.sto') # with RRA

        walking_speeds[i] = calculate_walking_speed(position_file)
    
    average_walking_speed = np.mean(walking_speeds)
    return average_walking_speed


def calculate_walking_speed(position_file):
    """calculate walking speed (change in COM x-position over time of stride) from a
    given COM position file.
    
    Args:
        position_file (str): filepath to COM position file.
        
    Returns:
        walking_speed (float): walking speed in m/s.
        
    """
    header_length = get_header_length(position_file)
    df = pd.read_csv(position_file, sep='\t', header=header_length, index_col=False)
    
    com_x = df['center_of_mass_X'].to_numpy()
    delta_com_x = com_x[-1] - com_x[0]
    
    time = df['time'].to_numpy()
    delta_time = time[-1] - time[0]
    
    walking_speed = delta_com_x/delta_time
    return walking_speed