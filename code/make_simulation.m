%% make_simulation.m
% make simulation for a given subject and trial.
% 
% Details:
%   scale model
%   run inverse kinematics
%   run inverse dynamics
%   run analysis to get center of mass position
%   run static optimization (developed by Uhlrich 2022)
%
% Args:
%   sim_home_dir (str): directory containing setup files for simulation.
%       This will also become the results directory.
%   walkID (str): trial identifier, e.g., 'walk01'.
%   walk_side (str): 'l' or 'r'.
%   start_time (double): simulation start time.
%   end_time (double): simulation end time.
%
% Returns:
%   none.
% 
% Notes:
%   if some steps have already been completed, they will be skipped. To
%   rerun:
%       - Scaling: delete scaled .osim file
%       - SO: delete results_forces.sto file
%       - any other process: delete out.log file
%   note that if any file is deleted, all steps proceeding that particular
%   step will be rerun. For example, if the whole pipeline has been run and
%   scaling is rerun, all other steps will be rerun. If the whole pipeline
%   has been run and IK is rerun, all steps following IK will be rerun.
%   This function is largely based on Apoorva Rajagopal's custom scripts.
%   See https://simtk.org/projects/full_body for Apoorva's original code.

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
function make_simulation(sim_home_dir, walkID, walk_side, start_time, end_time)
    import org.opensim.modeling.*;
    ModelVisualizer.addDirToGeometrySearchPaths('C:\OpenSim 4.1\Geometry');
    
    %%
    fprintf('Making simulation for %s...\n',walkID);
    overwrite_simulation = false;
    scripts_dir = cd;

    %% Get Subject ID
    sim_home_dir_split = strsplit(sim_home_dir, '\');
    subjectID = sim_home_dir_split{end};
    subjectID = subjectID(1:end-11); % remove '_SetupFiles'
    
    %% Scale model size
    scale_filename = [sim_home_dir '\Scale\' subjectID '_scaled.osim'];
    % only scale if model has not before been scaled
    if ~exist(scale_filename, 'file')
        disp('  Scaling model size...')
        cd([sim_home_dir, '/Scale'])
        [~,~] = system(['opensim-cmd run-tool ' subjectID ...
            '_Setup_ScaleTool.xml']);
        
        % check error log
        err = fileread('err.log');
        cd(scripts_dir)
        if ~isempty(err) || ~exist(scale_filename,'file')
            cd([sim_home_dir, '/Scale'])
            delete(scale_filename)
            cd(scripts_dir)
            error('Error in Scale. See error log.');
        end
        
        % warn user about iliacus and psoas muscle wrapping
        cnt = questdlg(['Open the scaled model in the OpenSim GUI. ' ...
                        'Maximally extend the hips and adjust the ' ...
                        'wrapping points of the iliacus and psoas ' ...
                        'muscles anteriorly, as needed.'], ...
                        'Check muscle wrapping', 'Complete', 'No', 'No');
        if strcmp(cnt, 'No')
            cd([sim_home_dir, '/Scale'])
            delete(scale_filename)
            cd(scripts_dir)
            error(['Ensure the iliacus and psoas muscles of the scaled' ...
                   ' model wrap correctly.'])
        end
        
        overwrite_simulation = true; % run future steps in simulation
    end
    
    %% Scale muscles
    scale_adjusted_filename = [sim_home_dir '\Scale\' subjectID '_adjusted.osim'];
    % only scale if model has not before been scaled
    if ~exist(scale_adjusted_filename, 'file')
        disp('  Scaling muscles...')
        cd([sim_home_dir, '/Scale'])
        
        generic_model = get_generic_model(sim_home_dir, subjectID);
        scaled_model = get_scaled_model(sim_home_dir, subjectID);
        subject_height = get_subject_height(subjectID);
        
        model_scaled_Fo = scale_optimal_force(generic_model, ...
            scaled_model, 1.70, subject_height);
        model_scaled_Fo = set_all_max_contraction_velocity(model_scaled_Fo, 15);
        % ^ since models that represent muscle paths as a single line tend to
        % overestimate length changes. See Thelen 2005, Arnold 2013, Ong 2019
        
        model_scaled_Fo_name = [subjectID '_adjusted'];
        model_scaled_Fo.setName(model_scaled_Fo_name);
        fclose('all');
        model_scaled_Fo.print([model_scaled_Fo_name '.osim']);
        cd(scripts_dir)
        
        overwrite_simulation = true;
    end
    
    %% IK
    ik_filename = [sim_home_dir '\IK\' walkID '_out.log'];

    if ~exist(ik_filename, 'file') || overwrite_simulation == true
        cd([sim_home_dir, '/IK'])
        disp('  Performing IK...');
        [~,~] = system(['opensim-cmd run-tool ' subjectID '_' walkID ...
            '_setupIK.xml']);
        
        % check error log
        err = fileread('err.log');
        cd(scripts_dir)
        if ~isempty(err)
            error('Error in IK. See error log.');
        end

        evaluate_ik(subjectID);
        
        cd([sim_home_dir, '/IK'])
        copyfile('out.log', [walkID '_out.log']); % save out log for each walk
        delete('out.log') % delete generic out log, which is now a copy
        cd(scripts_dir)
                
        overwrite_simulation = true;
    end

    %% ID
    id_filename = [sim_home_dir '\ID\' walkID '_out.log'];
    if ~exist(id_filename, 'file') || overwrite_simulation == true
        cd([sim_home_dir, '/ID']);
        disp('  Performing ID...');
        [~,~] = system(['opensim-cmd run-tool ' subjectID '_' walkID ...
            '_setupID.xml']);
        
        % check error log
        err = fileread('err.log');
        if ~isempty(err)
            cd(scripts_dir)
            error('Error in ID. See error log.');
        end
        
        cd([sim_home_dir, '/ID']);
        copyfile('out.log', [walkID '_out.log']); % save out log for each walk
        delete('out.log') % delete generic out log, which is now a copy
        cd(scripts_dir)
    end
    
    %% Analysis
    analysis_filename = [sim_home_dir '\Analysis\' walkID '_out.log'];
    if ~exist(analysis_filename, 'file') || overwrite_simulation == true
        cd([sim_home_dir, '/Analysis']);
        disp('  Performing BodyKinematics Analysis...');
        [~,~] = system(['opensim-cmd run-tool ' subjectID '_' ...
            walkID '_setupAnalysis.xml']);
        
        % check error log
        err = fileread('err.log');
        if ~isempty(err)
            cd(scripts_dir)
            error('Error in Analysis. See error log.');
        end
        
        cd([sim_home_dir, '/Analysis']);
        copyfile('out.log', [walkID '_out.log']); % save out log for each walk
        delete('out.log') % delete generic out log, which is now a copy
        cd(scripts_dir)
    end
        
    %% Custom SO
    so_filename = [sim_home_dir '\SO_Custom\' walkID '_out.log'];
    if ~exist(so_filename, 'file') || overwrite_simulation == true
        cd([sim_home_dir, '/SO_Custom']);
        disp('  Performing Custom SO...');
        
        run_static_optimization(sim_home_dir, subjectID, ...
            walkID, walk_side, start_time, end_time); % TODO bad out.log statements here
        
        % rename files to save for each walk
        copyfile('results_forces.sto', [walkID '_results_forces.sto']);
        copyfile('results_states.sto', [walkID '_results_states.sto']);
        copyfile('results_JointRxn_ReactionLoads.sto', ...
            [walkID '_results_JointRxn_ReactionLoads.sto']);
        
        % delete generic files, which are now copies
        delete('results_forces.sto')
        delete('results_states.sto')
        delete('results_JointRxn_ReactionLoads.sto')

        copyfile('out.log', [walkID '_out.log']); % save out log for each walk
        delete('out.log') % delete generic out log, which is now a copy
        cd(scripts_dir)
    end
    
    %% go back to scripts directory so you can start again
    cd(scripts_dir)
    fprintf('  Simulation made.\n')
end

%% helpers
function generic_model = get_generic_model(sim_home_dir, subjectID)
    % get filepath of generic model to be scaled for a given subject.
    %
    % Args:
    %   sim_home_dir (str): directory containing setup files for simulation.
    %   subjectID (str): subject identifier, e.g., 'C1'.
    %
    % Returns:
    %   generic_model (osim Model): generic model.
    
    import org.opensim.modeling.Model;
    fid = fopen([sim_home_dir '\Scale\' subjectID '_Setup_ScaleTool.xml'], 'r');
    f = fread(fid);
    s = char(f');
    start_identifier = strfind(s, '<model_file>');
    stop_identifier = strfind(s, '</model_file>');
    in_between_text = s(start_identifier+12 : stop_identifier-1);
    base_model_filepath = in_between_text;
    generic_model = Model(base_model_filepath);
end

function scaled_model = get_scaled_model(sim_home_dir, subjectID)
    % get filepath of scaled model for a given subject.
    %
    % Args:
    %   sim_home_dir (str): directory containing setup files for simulation.
    %   subjectID (str): subject identifier, e.g., 'C1'.
    %
    % Returns:
    %   scaled_model (osim Model): scaled model.
    
    import org.opensim.modeling.Model;
    scaled_model_filepath = [sim_home_dir '\Scale\' subjectID '_scaled.osim'];
    scaled_model = Model(scaled_model_filepath);
end

function subject_height = get_subject_height(subjectID)
    % read in table of subject IDs and heights and return height associated
    % with subjectID.
    %
    % Args:
    %   subjectID (str): subject identifier (e.g., 'C1').
    %
    % Returns:
    %   subject_height (double): subject height in meters.
    
    T = readtable('..\..\..\data\demographics.xlsx');
    subject_height = str2double(T.Var3(strcmp(T.Var1, subjectID)))/100;
end

function model_scaled_forces = scale_optimal_force(model_generic, ...
    model_scaled, height_generic, height_scaled)
    % return a model with muscle optimal forces scaled to subject.
    %
    % Args:
    %   model_generic (str): filepath of generic (unscaled) model.
    %   model_scaled (str): filepath of scaled model to adjust.
    %   height_generic (double): generic model height in meters.
    %   height_scaled (double): subject height in meters.
    %
    % Returns:
    %   model_scaled_forces (osim model): model with scaled peak isometric
    %       (optimal) muscle forces.
    
    mass_generic = get_model_mass(model_generic);
    mass_scaled = get_model_mass(model_scaled);
    
    % Handsfield et al. 2014
    v_total_generic = 47.05 * mass_generic * height_generic + 1289.6;
    v_total_scaled = 47.05 * mass_scaled * height_scaled + 1289.6;
    
    all_muscles_generic = model_generic.getMuscles();
    all_muscles_scaled = model_scaled.getMuscles();
    
    for i=0:all_muscles_generic.getSize()-1
        current_muscle_generic = all_muscles_generic.get(i);
        current_muscle_scaled = all_muscles_scaled.get(i);
        
        lmo_generic = current_muscle_generic.getOptimalFiberLength();
        lmo_scaled = current_muscle_scaled.getOptimalFiberLength();
        
        force_scale_factor = (v_total_scaled / v_total_generic) ...
            / (lmo_scaled / lmo_generic);
        
        current_muscle_scaled.setMaxIsometricForce(force_scale_factor ...
            * current_muscle_generic.getMaxIsometricForce());
    end
    
    model_scaled_forces = model_scaled;    
end

function mass_total = get_model_mass(model)
    % return total mass of a model.
    %
    % Args:
    %   model (osim model): model of interest.
    %
    % Returns:
    %   mass_total (double): total mass of the model, calculated from
    %       individual body masses.
    
    mass_total = 0;
    all_bodies = model.getBodySet();
    for i=0:all_bodies.getSize()-1
        current_body = all_bodies.get(i);
        mass_total = mass_total + current_body.getMass();
    end
end

function model_vmax = set_all_max_contraction_velocity(model, ...
    max_contraction_velocity)
    % update model max contraction velocities.
    %
    % Args:
    %   model (osim model): model to update.
    %   max_contraction_velocity (double): max contraction velocity value
    %       to update muscles.
    %
    % Returns:
    %   model_vmax (osim model): model with updated max contraction
    %       velocities.
    
    muscles = model.getMuscles();
    n_muscles = muscles.getSize();
    
    for i = 0:n_muscles-1
        current_muscle = muscles.get(i);
        current_muscle.setMaxContractionVelocity(max_contraction_velocity);
    end
    
    model_vmax = model;
end
