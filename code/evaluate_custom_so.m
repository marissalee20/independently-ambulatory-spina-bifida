%% evaluate_custom_so.m
% evaluate residuals of custom static optimization results.
%
% Details:
%   loop through static optimization results for a given walk and notify if
%   residuals are too large (residual forces >5% max ground reaction force,
%   residual moments >5% center of mass height * max ground reaction force.
%
% Args:
%   subjectID (str): subject identifier, e.g., 'C1'.
%   walkID (str): trial identifier, e.g., 'walk01'.
%
% Returns:
%   none.

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

%%
function evaluate_custom_so(subjectID, walkID)
    base_filepath = '..\';
    [max_ext_force, max_ext_mom] = get_ext_kinetics(base_filepath, ...
        subjectID, walkID);

    filename = [walkID '_results_forces.sto'];
    check_residuals(filename, max_ext_force, max_ext_mom);
end

%% helpers
function check_residuals(filename, max_ext_force, max_ext_moment)
    % check for large residuals (forces >5% max external force, moments >5%
    % center of mass height * max external force.
    %
    % Args:
    %   filename (str): filepath to results_forces output from static
    %       optimization.
    %   max_ext_force (double): max external force.
    %   max_ext_moment (double): center of mass height * max external
    %       force.
    %
    % Returns:
    %   none.
    
    scripts_dir = '..\..\..\code';
    residuals = readtable(filename, 'filetype', 'text');
    residual_forces = table2array(residuals(:,{'pelvis_tx_reserve', ...
        'pelvis_ty_reserve', 'pelvis_tz_reserve'}));
    residual_moments = table2array(residuals(:,{'pelvis_tilt_reserve', ...
        'pelvis_list_reserve', 'pelvis_rotation_reserve'}));
    
    total_residual_force = vecnorm(residual_forces, 2, 2);
    total_residual_moment = vecnorm(residual_moments, 2, 2);
    
    % get max/rms residual forces/moments
    max_res_force = max(total_residual_force);
    rms_res_force = rms(total_residual_force);
    max_res_moment = max(total_residual_moment);
    rms_res_moment = rms(total_residual_moment);
    
    % percents
    max_res_force_perc = max_res_force/max_ext_force*100;
    rms_res_force_perc = rms_res_force/max_ext_force*100;
    max_res_moment_perc = max_res_moment/max_ext_moment*100;
    rms_res_moment_perc = rms_res_moment/max_ext_moment*100;
    
    % check max residual force
    if max_res_force_perc > 5
        cnt = questdlg(sprintf...
                (['Max residual force too high at %.1f%% of max external' ...
                  ' force (should be <5%%). Continue?'], ...
                 max_res_force_perc), ...
                 'Max residual force is high', 'Yes', 'No', 'No');
        if ~strcmp(cnt,'Yes')
            cd(scripts_dir)
            error('Static optimization not accepted.');
        end
    end
    
    % check RMS residual force
    if rms_res_force_perc > 5
        cnt = questdlg(sprintf...
                (['RMS residual force too high at %.1f%% of max external' ...
                  ' force (should be <5%%). Continue?'], ...
                 rms_res_force_perc), ...
                 'RMS residual force is high', 'Yes', 'No', 'No');
        if ~strcmp(cnt,'Yes')
            cd(scripts_dir)
            error('Static optimization not accepted.');
        end
    end
    
    % check max residual moment
    if max_res_moment_perc > 5
        cnt = questdlg(sprintf...
                (['Max residual moment too high at %.1f%% of max ' ...
                  'external force * COM height (should be <5%%). ' ...
                  'Continue?'], max_res_moment_perc), ...
                 'Max residual moment is high', 'Yes', 'No', 'No');
        if ~strcmp(cnt,'Yes')
            cd(scripts_dir)
            error('Static optimization not accepted.');
        end
    end
    
    % check RMS residual moment
    if rms_res_moment_perc > 5
        cnt = questdlg(sprintf...
                (['RMS residual moment too high at %.1f%% of max ' ...
                  'external force * COM height (should be <5%%). ' ...
                  'Continue?'], rms_res_moment_perc), ...
                 'RMS residual moment is high', 'Yes', 'No', 'No');
        if ~strcmp(cnt,'Yes')
            cd(scripts_dir)
            error('Static optimization not accepted.');
        end
    end
    
end

function [max_ext_force, max_ext_mom] = get_ext_kinetics(sim_home_dir, ...
    subjectID, walkID)
    % extract max external force and center of mass * max external force.
    %
    % Details:
    %   note ext_moms aren't exactly external moments. They are COM height
    %   & external force magnitude.
    % 
    % Args:
    %   sim_home_dir (str): directory containing simulation setup files and
    %       RRA results.
    %   subjectID (str): subject identifier, e.g., 'C1'.
    %   walkID (str): trial identifier, e.g., 'walk01'.
    %
    % Returns:
    %   max_ext_force (double): maximum external force.
    %   max_ext_mom (double): maximum value of center of mass * maximum
    %       external force.
    
    filename = [sim_home_dir '\Analysis\PreRRA\results_analysis_' ...
        walkID '\analysis_' walkID '_BodyKinematics_pos_global.sto'];
    positions = readtable(filename, 'filetype', 'text');
    times = table2array(positions(:,{'time'}));
    start = times(1);
    stop = times(end);
    com_heights = table2array(positions(:,{'center_of_mass_Y'}));
    mean_com_height = mean(com_heights);
    
    % external forces
    filename = [sim_home_dir '\..\..\Data\' subjectID '\' subjectID '_' ...
        walkID '_grf_filtered.mot'];
    forces = readtable(filename, 'filetype', 'text');
    times = table2array(forces(:,{'time'}));
    
    forces_1 = table2array(forces(:,{'x1_ground_force_vx', ...
        'x1_ground_force_vy', 'x1_ground_force_vz'}));
    forces_2 = table2array(forces(:,{'x2_ground_force_vx', ...
        'x2_ground_force_vy', 'x2_ground_force_vz'}));
    forces_3 = table2array(forces(:,{'x3_ground_force_vx', ...
        'x3_ground_force_vy', 'x3_ground_force_vz'}));
    
    forces = forces_1 + forces_2 + forces_3;
        
    % keep only relevant times and resample to 101 points
    forces = interp1(times, forces, linspace(start, stop, 101));
    
    ext_forces = vecnorm(forces, 2, 2);
    
    max_ext_force = max(ext_forces);
    max_ext_mom = mean_com_height*max_ext_force;
    % ^ note not exactly external moment; magnitude of force * COM height
end
