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
import spm1d
from scipy.stats import ttest_ind

def print_ttest(td_df, sb_df, label, p_digits=3):
    """print p-value of a Welch's t-test and TD and SB means and standard deviations.
    
    Args:
        td_df (pd DataFrame): dataframe containing data, including those associated
                              with the given label, from subjects with typical
                              development.
        sb_df (pd DataFrame): dataframe containing data, including those associated
                              with the given label, from subjects with spina bifida.
        label (str): label of value to compare.
        p_digits (int): number of decimal places of p-value to print.
    
    Returns:
        none.
    """
    td_data = td_df[label]
    sb_data = sb_df[label]
    stat, p = ttest_ind(td_data, sb_data, equal_var=False)
    print(('%s: %.' + str(p_digits) + 'f') %(label, p))
    print('\tTD: %.1f (%.1f)' %(np.mean(td_data), np.std(td_data)))
    print('\tSB: %.1f (%.1f)' %(np.mean(sb_data), np.std(sb_data)))
    

def run_spm(data1, data2, label, plot=True, plot_stats=False, color1='k', color2='r'):
    """Conduct SPM two-sample t-test and plot results.
    
    Args:
        data1 (ndarray): n_waveforms x n_timepoints.
        data2 (ndarray): n_waveforms x n_timepoints.
        label (str): waveform name.
        plot (bool): True to plot data and standard deviations.
        plot_stats (bool): True to plot SPM results.
        color1 (color): color assigned to data1.
        color2 (color): color assigned to data2.
    
    Returns:
        ti (spm inference).
    
    """
    alpha = 0.05
    t     = spm1d.stats.ttest2(data1, data2, equal_var=False)
    ti    = t.inference(alpha, two_tailed=True, interp=True)
    
    if plot==True:
        # plot mean and SD:
        mean_data1 = np.mean(data1, axis=0)
        mean_data2 = np.mean(data2, axis=0)

        sd_data1 = np.std(data1, axis=0)
        sd_data2 = np.std(data2, axis=0)

        x_data1 = np.linspace(0, 100, data1.shape[1])
        x_data2 = np.linspace(0, 100, data2.shape[1])

        plt.fill_between(x_data1, mean_data1-sd_data1, mean_data1+sd_data1,
                         color=color1, alpha=0.5)
        plt.fill_between(x_data2, mean_data2-sd_data2, mean_data2+sd_data2,
                         color=color2, alpha=0.5)

        plt.plot(x_data1, mean_data1, color1, linewidth=2)
        plt.plot(x_data2, mean_data2, color2, linewidth=2)
    
    if plot_stats==True:
        plt.show()
        ti.plot()
        ti.plot_threshold_label(fontsize=8)
        ti.plot_p_values(size=10, offset_all_clusters=(0,0.9))
        plt.xlabel('Time (%)')
        plt.show()
    
    return ti


def create_bar(ti):
    """generate boolean array indicating spm significance regions.
    
    Args:
        ti (spm inference).
    
    Returns:
        bar (np array): boolean array.
    
    """
    values = ti.z
    thresh = ti.zstar
    bar = np.zeros([1, len(values)])
    
    for j in range(len(values)):
        if np.abs(values[j]) >= thresh:
            bar[0, j] = 1
    return bar


def get_bar_midpoints(bar):
    """identify midpoints of regions of significance.
    
    Args:
        bar (np array): boolean array.
    
    Returns:
        midpoints (np array): midpoints of regions of significance.
    
    """
    starts = [i+1 for i in range(len(bar)-1) if bar[i+1]>bar[i]]
    if bar[0]==1:
        starts.insert(0,0)
        
    ends = [i for i in range(len(bar)-1) if bar[i+1]<bar[i]]
    if bar[-1]==1:
        ends.append(len(bar)-1)
    
    midpoints = np.mean([np.array(starts), np.array(ends)], axis=0)
    return midpoints


def plot_bar(bar):
    """plot boolean 1D array.
    
    Args:
        bar (np array): boolean array.
    
    Returns:
        none.
    
    """
    fig = plt.imshow(bar, aspect='auto', cmap='binary', vmin=0, vmax=1)
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    return


def plot_sig(ti, color='r'):
    """plot significance bar for an spm t-test.
    
    Args:
        ti (spm inference).
        color (color): color in which to draw asterisks.
    
    Returns:
        none.
    
    """
    bar = create_bar(ti)
    plot_bar(bar)
    
    midpoints = get_bar_midpoints(bar[0])
    p_vals = ti.p

    for i in range(len(midpoints)):
        p_val = p_vals[i]
        print(p_val)
        if p_val < 0.001:
            plt.scatter([midpoints[i], midpoints[i], midpoints[i]], [-0.3, 0, 0.3], marker=(6,2,0), color=color)
        elif p_val < 0.01:
            plt.scatter([midpoints[i], midpoints[i]], [-0.15, 0.15], marker=(6,2,0), color=color)
        elif p_val < 0.05:
            plt.scatter(midpoints[i], 0, marker=(6,2,0), color=color)
    return
