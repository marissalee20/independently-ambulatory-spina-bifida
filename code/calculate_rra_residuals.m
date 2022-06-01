%% calculate_ik_errors.m
% calculate residuals relative to external force and external force * 
% center of mass height magnitudes across all trials.

% Copyright (c) 2022, Stanford Neuromuscular Biomechanics Laboratory
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% input subjects to run
subjects = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', ...
    'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'T1', 'T2', 'T3', 'T4', ...
    'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', 'T13', 'T14', ...
    'T15', 'T16'};
base_filepath = '..\simulation\';

% for each walk, record peak and RMS of proportion of res force to external
% force magnitude
peak_res_force_props = [];
rms_res_force_props = [];
% for each walk, record peak and RMS of proportion of res moment to
% external force * center of mass height magnitude
peak_res_mom_props = [];
rms_res_mom_props = [];

for i = 1:length(subjects)
    subjectID = char(subjects(i));
    subject_datapath = [base_filepath '..\Data\' subjectID];
    walkIDs = get_walkIDs(subject_datapath);
    
    for j = 1:length(walkIDs)
        walkID = char(walkIDs{j});
        
        [res_forces, res_moms] = get_residuals(base_filepath, subjectID, ...
            walkID);
        [max_ext_force, max_ext_mom] = get_external_kinetics(base_filepath, ...
            subjectID, walkID);
        
        res_force_prop = reshape(res_forces, 1, 101) / max_ext_force;
        peak_res_force_prop = max(res_force_prop);
        rms_res_force_prop = compute_rmse(res_force_prop, 0);
        
        res_mom_prop = reshape(res_moms, 1, 101) / max_ext_mom;
        peak_res_mom_prop = max(res_mom_prop);
        rms_res_mom_prop = compute_rmse(res_mom_prop, 0);
        
        peak_res_force_props = [peak_res_force_props, peak_res_force_prop];
        rms_res_force_props = [rms_res_force_props, rms_res_force_prop];
        peak_res_mom_props = [peak_res_mom_props, peak_res_mom_prop];
        rms_res_mom_props = [rms_res_mom_props, rms_res_mom_prop];
    end
end

fprintf('Max Force Res: %.1f%% ± %.1f%% External Force\n', ...
    100*mean(peak_res_force_props), 100*std(peak_res_force_props));
fprintf('Res Force RMSE: %.1f%% ± %.1f%% External Force\n', ...
    100*mean(rms_res_force_props), 100*std(rms_res_force_props));
fprintf('Max Moment Res: %.1f%% ± %.1f%% COM Height * External Force\n', ...
    100*mean(peak_res_mom_props), 100*std(peak_res_mom_props));
fprintf('Res Moment RMSE: %.1f%% ± %.1f%% COM Height * External Force\n', ...
    100*mean(rms_res_mom_props), 100*std(rms_res_mom_props));

%% Helpers
function walkIDs = get_walkIDs(subject_datapath)
    % read notes file for in subject_datapath and return info.
    %
    % Args:
    %   subject_datapath (str): filepath of subject data directory.
    %
    % Returns:
    %   walkIDs (cell array): list of walk identifiers.
    
    subject_datapath_split = strsplit(subject_datapath,'\');
    subjectID = subject_datapath_split{end};
    filename = [subject_datapath '\' subjectID ' notes.xlsx'];
    T = readtable(filename);
    walkIDs = T.trial;
end

function [res_forces, res_moms] = get_residuals(base_filepath, subjectID, walkID)
    % extract residuals for the given trial.
    %
    % Args:
    %   base_filepath (str): filepath of simulation directory.
    %   subjectID (str): subject identifier.
    %   walkID (str): walk identifier.
    %
    % Returns:
    %   res_forces (array): vector of residual force magnitude over the
    %                       trial.
    %   res_moms (array): vector of residual moment magnitude over the
    %                     trial.

    filename = [base_filepath subjectID '_SetupFiles\RRA\results_rra_' ...
        walkID '\rra_' walkID '_Actuation_force.sto'];
    residuals = readtable(filename, 'filetype', 'text');
    residuals_forces = table2array(residuals(:, {'FX','FY','FZ'}));
    residuals_moments = table2array(residuals(:, {'MX','MY','MZ'}));
    
    res_forces = vecnorm(residuals_forces, 2, 2);
    res_moms = vecnorm(residuals_moments, 2, 2);
    
    % resample to 101 points
    res_forces = interp1(1:length(res_forces), res_forces, ...
        linspace(1, length(res_forces), 101));
    res_moms = interp1(1:length(res_moms), res_moms, ...
        linspace(1, length(res_moms), 101));
end

function [max_ext_force, max_ext_mom] = get_external_kinetics(base_filepath, ...
    subjectID, walkID)
    % get maximum external force and maximum external force * center of
    % mass height values across the trial.
    %
    % Args:
    %   base_filepath (str): filepath of simulation directory.
    %   subjectID (str): subject identifier.
    %   walkID (str): walk identifier.
    %
    % Returns:
    %   max_ext_force (float): maximum external force magnitude.
    %   max_ext_mom (float): maximum value of external force magnitude *
    %                        center of mass height.
    %
    % Notes:
    %   ext_moms aren't exactly external moments. They are COM height &
    %   external force magnitude.
    
    filename = [base_filepath subjectID ...
        '_SetupFiles\Analysis\PostRRA\results_analysis_' ...
        walkID '\analysis_' walkID '_BodyKinematics_pos_global.sto'];
    positions = readtable(filename, 'filetype', 'text');
    times = table2array(positions(:, {'time'}));
    start = times(1);
    stop = times(end);
    com_heights = table2array(positions(:, {'center_of_mass_Y'}));
    % resample to 101 points
    mean_com_height = mean(com_heights);
    
    % external forces
    filename = [base_filepath '..\Data\' subjectID '\' subjectID '_' ...
        walkID '_grf_filtered.mot'];
    forces = readtable(filename, 'filetype', 'text');
    times = table2array(forces(:, {'time'}));
    forces_1 = table2array(forces(:, {'x1_ground_force_vx', ...
        'x1_ground_force_vy', 'x1_ground_force_vz'}));
    forces_2 = table2array(forces(:, {'x2_ground_force_vx', ...
        'x2_ground_force_vy', 'x2_ground_force_vz'}));
    forces_3 = table2array(forces(:, {'x3_ground_force_vx', ...
        'x3_ground_force_vy', 'x3_ground_force_vz'}));
    
    forces = forces_1 + forces_2 + forces_3;
        
    % keep only relevant times and resample to 101 points
    forces = interp1(times, forces, linspace(start, stop, 101));
    
    ext_forces = vecnorm(forces, 2, 2);
    
    max_ext_force = max(ext_forces);
    max_ext_mom = mean_com_height * max_ext_force;
    % note not exactly external moment; magnitude of force * COM height
end

function rmse = compute_rmse(yhat, y)
    % compute the root-mean-square error.
    %
    % Args:
    %   yhat (array): estimates.
    %   y (array): observations.
    %
    % Returns:
    %   rmse (float): root-mean-square error.
    
    rmse = sqrt(mean((y - yhat).^2));
end