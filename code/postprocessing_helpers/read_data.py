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

def get_header_length(storage_filename):
    """Find header length of storage file.
    
    Args:
        storage_filename (str): filename.
        
    Returns:
        header_length (int): lines of header.
        
    """
    with open(storage_filename) as f:
        lines = f.readlines()
    
    end_idx = -1
    current_line = 0
    while end_idx==-1:
        if 'endheader' in lines[current_line]:
            end_idx = current_line
        current_line += 1
        if current_line > len(lines):
            raise Exception('Unable to find end of header.')

    header = lines[:end_idx+1]
    while '\n' in header:
        header.remove('\n')
        
    header_length = len(header)
    return header_length


def get_walkIDs(notes_filename):
    """extract walk identifiers from the subject's note file.
    
    Args:
        notes_filename (str): filepath to notes file containing walk identifiers.
        
    Returns:
        walkIDs (list): list of walk identifiers (e.g., 'Walk01')
    """
    notes = pd.read_excel(notes_filename)
    walkIDs = list(notes['trial'])
    return walkIDs


def get_bone_outcomes(df, bone_outcome_names):
    """Extract analyzed-side-specific bone outcomes from the demographics table.
    
    Args:
        df (pd DataFrame): dataframe containing analyzed sides and bone outcomes per
                           side, per subject.
        bone-outcome_names (list): list of strings of names of bone outcomes.
    
    Returns:
        analyzed_side_outcomes (pd DataFrame): dataframe containing analyzed-side-
                                               specific bone outcomes (cols) for each
                                               subject (rows).
    """
    n_subjects = len(df)
    n_outcomes = len(bone_outcome_names)
    analyzed_side_outcomes = np.zeros((n_subjects, n_outcomes))
    for i in range(n_subjects):
        side = df['analyzed side'].iloc[i]
        for j in range(n_outcomes):
            outcome_name = side + ' ' + bone_outcome_names[j]
            analyzed_side_outcomes[i, j] = df[outcome_name].iloc[i]
    
    analyzed_side_outcomes = pd.DataFrame(analyzed_side_outcomes, columns=bone_outcome_names)
    analyzed_side_outcomes['subject'] = list(df['subject'])
            
    return analyzed_side_outcomes


def get_analyzed_side_pf_strength(demographics):
    """return a list of plantar flexor strengths associated with the analyzed limbs.
    
    Args:
        demographics (pd DataFrame): subject demographic information. Includes
                                     'analyzed side', 'L plantar flexor strength',
                                     and 'R plantar flexor strength' column names.
    
    Returns:
        pf_strengths (list): list of plantar flexor strengths associated with the
                             analyzed limbs.
        analyzed_sides (list): list of analyzed limbs.
    
    """
    analyzed_sides = demographics['analyzed side']
    left_pf_strength = demographics['L plantar flexor strength']
    right_pf_strength = demographics['R plantar flexor strength']
    
    pf_strengths = []
    for i in range(len(analyzed_sides)):
        analyzed_side = analyzed_sides.iloc[i]
        if analyzed_side == 'L':
            pf_strengths.append(left_pf_strength.iloc[i])
        else:
            pf_strengths.append(right_pf_strength.iloc[i])
    
    return pf_strengths, analyzed_sides
