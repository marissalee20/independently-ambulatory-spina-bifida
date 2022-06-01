%% run_iterative_rra.m
% Iteratively run RRA, starting with scaled model and IK results, to
% achieve residual reduction model & adjusted kinematics.
%
% Details:
%   Aims for residual forces <5% max ground reaction force and residual
%   moments <2% center of mass height * max ground reaction force. Requires
%   some user decisions intermittently (user will be prompted).
% 
% Args:
%   sim_home_dir (str): directory containing simulation setup files and RRA
%       results.
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
function run_iterative_rra(sim_home_dir, subjectID, walkID)
    import org.opensim.modeling.*;
    scripts_dir = cd;
    cd([sim_home_dir, '/RRA'])
    disp('  Performing RRA...')

    incomplete = true;
    iteration = 1;
    
    % reset task file to default
    task_filename = 'rra_tasks_walk.xml';
    template_path = ['..\..\Template_SetupFiles\RRA\' task_filename];
    
    status = copyfile(template_path, task_filename);
    if status == 0
        cd(scripts_dir)
        error('Files not copied')
    end
    
    % replace pelvis mass center position in actuator file
    scale_model_name = ['../Scale/' subjectID '_scaled.osim'];
    update_actuators_file_pelvis_com(scale_model_name);
    
    last_mass_change = 1;
        
    while incomplete
        fprintf('    Iteration %i\n', iteration)
        
        % read in most recent model (scaled model if 1st iteration)
        model_name = get_model_name(iteration, subjectID, walkID);
        
        % adjust setup file to include proper model_file
        setup_filename = [subjectID '_' walkID '_setupRRA.xml'];
        update_setup_model(setup_filename, model_name);
        
        % replace pelvis mass center position in actuator file
        update_actuators_file_pelvis_com(model_name);
        
        % save a copy of the model so we can revert if updated model
        % kinematics don't track coordinates well
        copyfile(model_name, 'revert_model.osim');
        
        % run RRA
        [~,~] = system(['opensim-cmd run-tool ' setup_filename]);
        err = fileread('err.log');
        if ~isempty(err)
            cd(scripts_dir);
            error('Error in RRA. See error log.');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% check residuals and errors %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % check positional error. If not tracking well, double tracking 
        % weight on that coordinate and rerun with reverted model.
        err_filename = ['results_rra_' walkID '/rra_' walkID '_pErr.sto'];
        task_filename = 'rra_tasks_walk.xml';
        [pErr_flag, coord_scale] = check_tracking(err_filename);
        
        if pErr_flag == 1 % if not tracking
            fprintf('      Not tracking coordinates well. Backtracking...\n')
            
            % replace values in rra_tasks_walk based on poorly-tracked coords
            update_tasks(task_filename, coord_scale);
            
            % revert to old model for next loop by overwriting rra output
            % model with revert model
            copyfile('revert_model.osim', model_name);
            revert_model_mass(model_name, last_mass_change);
            last_mass_change = 0;
            
            % don't update iteration number (so that correct model of
            % scaled vs rra output is selected)
            
        else % if tracking
            res_err_flag = check_residuals(sim_home_dir, subjectID, walkID);
            
            if res_err_flag == 1 % if residuals too high
                fprintf('      High residuals. Iterating...\n')

                % update tracking weights
                update_tasks(task_filename, coord_scale);

                % update model's torso mass
                model_name = [subjectID '_' walkID '_rra.osim'];
                last_mass_change = update_model_mass(model_name, 'out.log');

                % update iteration number
                iteration = iteration + 1;
                if iteration > 3
                    cnt = questdlg(['Attempted 3+ iterations and could' ...
                        ' not get ideal RRA results. Continue with ' ...
                        'most recent results and move onto static ' ...
                        'optimization?'], 'Hit 3+ iterations', ...
                        'Continue to SO', 'Iterate RRA', 'Iterate RRA');
                    if strcmp(cnt, 'Continue to SO')
                        fprintf('      RRA success\n')
                        cd(scripts_dir);
                        incomplete = false;
                    end
                end
            else
                cd(scripts_dir);
                % if residuals are good, success.
                fprintf('      RRA success\n')
                incomplete = false; % stop while loop because success
            end
        end
    end
end

%% helpers

function update_actuators_file_pelvis_com(model_name)
    % update pelvis center of mass (COM) in actuator file
    %
    % Args:
    %   model_name (str): filepath of model from which to identify COM.
    %
    % Returns:
    %   none.
    
    model = org.opensim.modeling.Model(model_name);

    body_set = model.getBodySet();
    com = body_set.get('pelvis').getMassCenter(); 
    com_x = num2str(com.get(0));
    com_y = num2str(com.get(1));
    com_z = num2str(com.get(2));

    % define replace function inputs
    actuator_filename = 'rra_actuators.xml';
    string_start = '<point>';
    string_end = '</point>';
    replacement = [com_x ' ' com_y ' ' com_z];
    
    replace_in_text(actuator_filename, string_start, string_end, replacement)
end

function replace_in_text(filename, string_start, string_end, replacement)
    % replace text between and including string_start and string_end with
    % replacement in filename.
    %
    % Args:
    %   filename (str): text filepath.
    %   string_start (str): start of string to replace.
    %   string_end (str): end of string to replace.
    %   replacement (str): replacement string.
    %
    % Returns:
    %   none.
    
    % Read original file, search and replace, overwrite template
    fid  = fopen(filename,'r');
    text = fread(fid,'*char')';
    fclose(fid);

    text = replaceBetween(text, string_start, string_end, replacement);
    
    fid = fopen(filename,'w');
    fprintf(fid,'%s',text);
    fclose(fid);
end

function model_name = get_model_name(iteration, subjectID, walkID)
    % get model name.
    %
    % Args:
    %   iteration (int): iteration number. Gets scaled model if RRA has not
    %       been run yet, most recent RRA model if RRA has been run.
    %   subjectID (str): subject identifier, e.g., 'C1'.
    %   walkID (str): trial identifier, e.g., 'walk01'.
    %
    % Returns:
    %   model_name (str): path of model of interest.
    
    if iteration == 1
        model_name = ['../Scale/' subjectID '_scaled.osim'];
    else
        model_name = [subjectID '_' walkID '_rra.osim'];
    end
end

function update_setup_model(setup_filename, model_name)
    % update RRA setup file with correct model_name
    %
    % Args:
    %   setup_filename (int): setup file name containing a model_file
    %   model_name (str): replacement model name.
    %
    % Returns:
    %   none.
    
    string_start = '<model_file>';
    string_end = '</model_file';
    replace_in_text(setup_filename, string_start, string_end, model_name)
end

function [pErr_flag, coord_scale] = check_tracking(err_filename)
    % checks for tracking errors from RRA results. If tracking errors
    % exist, plan to revert non-tracking weights to former (double) values.
    % If no tracking errors exist, plan to halve all weights for rotational
    % coordinates that tracked well.
    %
    % Args:
    %   err_filename (str): filepath to err file
    %
    % Returns:
    %   pErr_flag (0/1): 1 if tracking error(s) exist.
    %   coord_scale (table): recommended scale factors for each coordinate
    %       task.
    
    pErr_flag = 0; % default flag for errors
        
    % get translational and rotational errors
    perrs = readtable(err_filename, 'filetype', 'text');
    perrs = removevars(perrs, {'time'}); % time doesn't matter
    coords = perrs.Properties.VariableNames;
    translation_coords = {'pelvis_tx','pelvis_ty','pelvis_tz'};
    perrs = table2array(perrs);
    max_errs = max(abs(perrs));
    rms_errs = rms(perrs);
    
    % read initial coordinate weights
    n_coords = length(coords);
    coord_scale = ones(1, n_coords);
    coord_scale = array2table(coord_scale, 'VariableNames', coords);
    
    for i = 1:n_coords
        coord_name = char(coords{i});
        max_err = max_errs(i);
        rms_err = rms_errs(i);
        
        % translational: max & rms should be < 5 cm and < 4 cm
        if any(strcmp(translation_coords, coord_name))
            % if errors are too large, double tracking weight, raise flag
            if max_err > 0.05 || rms_err > 0.04
                coord_scale.(coord_name) = 2;
                pErr_flag = 1;
            end
        % rotational: max and rms should be < 5 deg each
        else
            % if errors are too large, double tracking weights, raise flag
            if max_err > deg2rad(5) || rms_err > deg2rad(5)
                coord_scale.(coord_name) = 2;
                pErr_flag = 1;
            end
        end
    end
    
    % if no tracking issues, decrease weight of rotational coordinates
    % tracked w/in 1 deg and translational coordinates tracked w/in 1 cm by
    % factor of 2
    if pErr_flag == 0
        for i = 1:n_coords
            coord_name = char(coords{i});
            max_err = max_errs(i);
            % if rotational and trackign w/in 1 deg
            if ~any(strcmp(translation_coords,coord_name))
                if max_err < deg2rad(1)
                    coord_scale.(coord_name) = 0.5;
                end
            % if translational and tracking w/in 1 cm
            else
                if max_err < 0.01
                    coord_scale.(coord_name) = 0.5;
                end
            end
        end
    end
end

function update_tasks(task_filename, coord_scale)
    % update task file with scaled values
    %
    % Args:
    %   task_filename (str): filepath to task file (for RRA setup)
    %   coord_scale (table): recommended scale factors for each coordinate
    %       task.
    %
    % Returns:
    %   none.
    
    % read old weights
    old_weights = get_old_task_weights(task_filename);
    coord_names = old_weights.Properties.VariableNames;
    
    % calculate and write new weights by using scaling factors
    for i = 1:length(coord_names)
        var_name = char(coord_names(i));
        new_weight = old_weights.(var_name)*coord_scale.(var_name);
        
        % write
        string_start = [var_name '"> <weight>'];
        string_end = '</weight>';
        replace_in_text(task_filename, string_start, string_end, ...
            num2str(new_weight))
    end
end

function weights = get_old_task_weights(task_filename)
    % extract task weights from task file.
    %
    % Args:
    %   task_filename (str): filepath to task file (for RRA setup)
    %
    % Returns:
    %   weights (table): task weights.
    
    % read in task weights to table
    fid  = fopen(task_filename, 'r');
    text = fread(fid, '*char')';
    fclose(fid);
    
    % ignore text above "</defaults>"
    text = extractAfter(text, '</defaults>');
    
    % for remainder of document, extract joint names and weights
    subtexts = extractBetween(text, 'name="', '</weight>');
    n_coords = length(subtexts);
    var_names = cell(1,n_coords);
    weights = zeros(1,n_coords);
    for i = 1:length(subtexts)
        subtext = subtexts{i};
        var_names{1,i} = extractBefore(subtext, '">');
        weights(1,i) = str2double(extractAfter(subtext, '<weight>'));
    end
    
    weights = array2table(weights, 'VariableNames', var_names);
end

function res_err_flag = check_residuals(sim_home_dir, subjectID, walkID)
    % check for residuals that are too large (forces > 5% max external
    % force, moments > 5% center of mass * max external force)
    %
    % Args:
    %   sim_home_dir (str): directory containing simulation setup files and
    %       RRA results.
    %   subjectID (str): subject identifier, e.g., 'C1'.
    %   walkID (str): trial identifier, e.g., 'walk01'.
    %
    % Returns:
    %   res_err_flag (0/1): 1 if residual errors are too large.
    
    force_filename = ['results_rra_' walkID '/rra_' walkID ...
        '_Actuation_force.sto'];
    
    res_err_flag = 0; % default error flag
    
    residuals = readtable(force_filename, 'filetype', 'text');
    residual_forces = table2array(residuals(:,{'FX','FY','FZ'}));
    residual_moments = table2array(residuals(:,{'MX','MY','MZ'}));
    
    total_residual_force = vecnorm(residual_forces, 2, 2);
    total_residual_moment = vecnorm(residual_moments, 2, 2);
    
    % get max/rms residual forces/moments
    max_res_force = max(total_residual_force);
    rms_res_force = rms(total_residual_force);
    max_res_moment = max(total_residual_moment);
    rms_res_moment = rms(total_residual_moment);
    
    % get external kinetics
    [max_ext_force, max_ext_moment] = ...
        get_external_kinetics(sim_home_dir, subjectID, walkID);
    
    % percents
    max_res_force_perc = max_res_force/max_ext_force*100;
    rms_res_force_perc = rms_res_force/max_ext_force*100;
    max_res_moment_perc = max_res_moment/max_ext_moment*100;
    rms_res_moment_perc = rms_res_moment/max_ext_moment*100;
    
    fprintf('\t\t\tMax residual force:  %.1f%% max ext force\n', ...
        max_res_force_perc)
    fprintf('\t\t\tRMS residual force:  %.1f%% max ext force\n', ...
        rms_res_force_perc)    
    fprintf('\t\t\tMax residual moment: %.1f%% max ext mom\n', ...
        max_res_moment_perc)
    fprintf('\t\t\tRMS residual moment: %.1f%% max ext mom\n', ...
        rms_res_moment_perc)

    % check values for issues
    if max_res_force_perc > 5 || rms_res_force_perc > 5 || ...
            max_res_moment_perc > 2 || rms_res_moment_perc > 2
        res_err_flag = 1;
    end
end

function [max_ext_force, max_ext_mom] = get_external_kinetics(sim_home_dir, ...
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

function last_mass_change = update_model_mass(model_name, log_filename)
    % update model body masses based on mass change in out.log.
    % 
    % Args:
    %   model_name (str): filepath of OpenSim model for which to update the
    %       model mass using RRA log file results.
    %   log_filename (str): filepath of RRA log file containing model mass
    %       change.
    %
    % Returns:
    %   none.
    
    out_log = fileread(log_filename);
    mass_change = str2double(out_log(regexpi(out_log, 'total mass change: ', ...
        'end'):regexpi(out_log, 'total mass change: .?[0-9]+[.][0-9]+', 'end')));
    
    % adjust body masses based on recommended total mass change
    model = org.opensim.modeling.Model(model_name);
%     if abs(mass_change) > 2 % reduce mass adjustment if > 2 kg change
%         mass_change = 2*sign(mass_change);
%     end
    [model, last_mass_change] = update_body_masses(model, ...
        mass_change);
    model.print(model_name);
end

function revert_model_mass(model_name, last_mass_change)
    % revert model mass to previous model mass.
    %
    % Args:
    %   model_name (str): path to model file.
    %   last_mass_change (float): inverse of mass change to apply.
    %
    % Returns:
    %   none.
    
    model = org.opensim.modeling.Model(model_name);
    [model, ~] = update_body_masses(model, ...
        -last_mass_change);
    model.print(model_name);
end

function [model, mass_change] = update_body_masses(model, mass_change)
    % update model body masses using total mass_change.
    % 
    % Args:
    %   model (osim model): model for which to update body masses using RRA
    %       suggested mass change.
    %   mass_change (double): RRA suggested mass change.
    %
    % Returns:
    %   model (osim model): updated model with body masses scaled to
    %       reflect the RRA mass change.
    
    orig_mass_total = get_model_mass(model);
    new_mass_total = orig_mass_total + mass_change;
    fprintf('\t\t\tModel mass: %.2f kg\n', new_mass_total)
    mass_scale_factor = new_mass_total/orig_mass_total;
    
    all_bodies = model.getBodySet();
    for i = 0:all_bodies.getSize()-1
        orig_mass_of_body = all_bodies.get(i).getMass();
        new_mass_of_body = orig_mass_of_body*mass_scale_factor;
        all_bodies.get(i).setMass(new_mass_of_body);
    end
end

function mass_total = get_model_mass(model)
    % get total mass of OpenSim model.
    % 
    % Args:
    %   model (osim model): model for which to extract total mass.
    %
    % Returns:
    %   mass_total (double): model total mass, based on individual body
    %       masses.
    
    mass_total = 0;
    all_bodies = model.getBodySet();
    for i=0:all_bodies.getSize()-1
        curr_body = all_bodies.get(i);
        mass_total = mass_total + curr_body.getMass();
    end
end
