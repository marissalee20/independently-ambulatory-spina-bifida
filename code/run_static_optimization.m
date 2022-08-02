%% run_static_optimization
% run static optimization for a given trial, accounting for tendon
% compliance and passive muscle forces.
%
% Details:
%   run static optimization, accounting for tendon compliance and passive
%   forces.
%
% Args:
%   sim_home_dir (str): directory containing setup files for simulation.
%       This will also become the results directory.
%   subjectID (str): subject identifier, e.g., 'C1'.
%   walkID (str): trial identifier, e.g., 'walk01'.
%   walk_side (str): 'l' or 'r'.
%   start_time (double): simulation start time.
%   end_time (double): simulation end time.
%
% Returns:
%   none.
% 
% Notes:
%   this function is largely a pared down version of Scott Uhlrich's custom
%   static optimization. See https://simtk.org/projects/coordretraining/
%   for Scott's original code.

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
function run_static_optimization(sim_home_dir, subjectID, walkID, ...
    walk_side, start_time, end_time)

    import org.opensim.modeling.*
    
    %% setup
    % paths
    grf_filepath = [sim_home_dir '\..\..\data\' subjectID '\' ...
        subjectID '_' walkID '_grf_filtered.mot'];
    coordinate_filepath = [sim_home_dir '\IK\' subjectID '_' walkID '_ik.mot'];
    actuation_force_filepath = [sim_home_dir '\ID\results_' walkID ...
        '\inverse_dynamics_' walkID '.sto'];
    output_filepath = [sim_home_dir '\SO_Custom\'];
    model_name = [sim_home_dir '\Scale\' subjectID '_adjusted.osim'];

    % degrees of freedom to ignore during moment matching constraint. This
    % is specific to the model.
    fixed_dofs = {'knee_angle_r_beta', 'knee_angle_l_beta'};

    % External Forces Definitions
    grf_names = {'GRF_FP1', 'GRF_FP2', 'GRF_FP3'};
    if strcmp(walk_side, 'l')
        grf_application_bodies = {'calcn_r', 'calcn_l', 'calcn_r'};
    elseif strcmp(walk_side, 'r')
        grf_application_bodies = {'calcn_l', 'calcn_r', 'calcn_l'};
    else
        error("Unrecognized input for walkSide. Must be 'l' or 'r'");
    end

	grf_expression_frames =  {'ground', 'ground', 'ground'};
	grf_identifiers = {'1_ground_force_v', '2_ground_force_v', ...
        '3_ground_force_v'};
	gr_point_expression_frames = {'ground', 'ground', 'ground'};
    gr_point_identifiers = {'1_ground_force_p', '2_ground_force_p', ...
        '3_ground_force_p'};

    % joint reaction Fields
    jrxn_joints = {'walker_knee_r', 'ankle_r', 'walker_knee_l', 'ankle_l'}; % these joints
    jrxn_bodies = {'child', 'parent', 'child', 'parent'}; % on these bodies
    jrxn_expression_frames = {'child', 'parent', 'child', 'parent'}; % in these frames

    % create output directory and log output
    mkdir(output_filepath);
    diary([output_filepath 'out.log']) ;
    fprintf(['\n\n\tBeginning Optimization for %s from %.2f to %.2f ' ...
             'seconds. \n\n'], walkID, start_time, end_time);
    
    %% organize model & get free coordinates
    model = Model(model_name);
    coords = model.getCoordinateSet();
    n_coords = coords.getSize() ;
    forces = model.getForceSet() ;
    n_forces = forces.getSize ;
    
    % get free coordinates (not knee_angle_beta). Note specific to
    % Rajagopal2015 model.
    free_coord_names = {} ; 
    for i = 0:n_coords-1
        if coords.get(i).get_locked == 0 ... 
           && coords.get(i).get_prescribed == 0 ... 
           && sum(strcmp(fixed_dofs, char(coords.get(i).getName))) == 0
            free_coord_names{end+1} = char(coords.get(i).getName) ;
        end
    end
    n_free_coords = length(free_coord_names) ;
    
    %% organize actuators
    % remove coordinate actuators. NOTE: you must remove all coordinate
    % actuators.
    force_names_to_delete = {}; 
    for i = 1:n_forces
        actuator_name = char(forces.get(i-1).getName()); 
        if contains(convertCharsToStrings(actuator_name), ...
                ["lumbar", "shoulder", "elbow", "pro_sup", "wrist"])
             force_names_to_delete{end+1} = char(forces.get(i-1).getName);
        end        
    end
    for i = 1:length(force_names_to_delete)
        forces.remove(forces.get(force_names_to_delete{i}));
    end
    
    % append residual and reserve actuators to free coords
    for i = 0:n_coords-1
        if sum(strcmp(char(coords.get(i).getName), free_coord_names))
            new_actuator = CoordinateActuator(char(coords.get(i).getName));
            new_actuator.setName([char(coords.get(i).getName()) '_reserve']);
            new_actuator.set_min_control(-Inf);
            new_actuator.set_max_control(Inf);
            coord_name = char(coords.get(i).getName);
            if contains(coord_name,'pelvis') || contains(coord_name,'lumbar') 
                new_actuator.set_optimal_force(100); % cheap for dynamic consistency
            else
                new_actuator.set_optimal_force(1);
            end    
            model.addForce(new_actuator);

            % add prescribed controllers for coord actuators (necessary for
            % joint reaction force to work properly)

            % construct piecewise constant function
            constFxn = Constant(0);
            constFxn.setName([char(coords.get(i).getName) '_constFxn']);

            % construct prescribed controller
            pController = PrescribedController();
            pController.setName([char(coords.get(i).getName) '_controller']);
            pController.addActuator(new_actuator);
            pController.prescribeControlForActuator(0,constFxn);
            model.addController(pController);
        end % freeCoords
    end
    
    controllers = model.getControllerSet;
    n_controllers = controllers.getSize;
    actuators = forces.updActuators;
    
    %% organize muscles
    muscles = model.updMuscles();
    n_muscles = muscles.getSize();
    muscleNames = cell(1, n_muscles);
    for i = 1:n_muscles
        muscleNames{i} = char(muscles.get(i-1).getName);
    end
    
    %% append ground reaction forces
    data_source = Storage(grf_filepath);
    for i = 1:length(grf_application_bodies) % nExternalForces
        new_force = ExternalForce();
        new_force.setName(grf_names{i});
        new_force.set_applied_to_body(grf_application_bodies{i});
        new_force.set_force_expressed_in_body(grf_expression_frames{i});
        new_force.set_force_identifier(grf_identifiers{i});
        new_force.set_point_expressed_in_body(gr_point_expression_frames{i});
        new_force.set_point_identifier(gr_point_identifiers{i});
        new_force.setDataSource(data_source);
        model.addForce(new_force);
    end
    
    % % % % NO MORE EDITING THE MODEL AFTER THIS POINT % % % % %
    
    %% get model state
    state = model.initSystem;
    n_forces = forces.getSize;
    n_actuators = actuators.getSize;
    n_muscles = muscles.getSize;
    
    % state names
    states = model.getStateVariableNames;
    n_states = getSize(states);
    state_names = cell(n_states,1);
    for i = 1:n_states
        full_state_name = char(states.get(i-1));
        inds_slash = strfind(full_state_name,'/');
        state_names{i} = full_state_name(inds_slash(end-1)+1:end);
    end
    
    % force names
    force_names = cell(n_forces,1);
    for i = 1:n_forces
       force_names{i} = char(forces.get(i-1).getName());
    end

    % actuator names
    actuator_names = cell(n_actuators,1);
    for i = 1:n_actuators
       actuator_names{i} = char(actuators.get(i-1).getName());
    end
    
    % indices referring to coord velocities in state vector (every coord
    % that isn't fixed or constrained)
    coord_vel_indices = zeros(n_free_coords,1);
    coord_vel_names = cell(n_free_coords,1);
    for i = 1:n_free_coords
        coord_vel_indices(i) = find(strcmp([free_coord_names{i} ...
            '/speed'], state_names));
        coord_vel_names{i} = state_names{coord_vel_indices(i)};
    end
    
    %% get coordinates
    coordinate_sto = Storage(coordinate_filepath);
    
    time_coords = ArrayDouble();
    coordinate_sto.getTimeColumn(time_coords);
    time_coords = str2num(time_coords);
    sample_freq = 1 / (time_coords(2) - time_coords(1)); % coordinate sample rate

    coord_names_q = cell(1, n_coords);
    q = zeros(length(time_coords), n_coords);
    free_coord_inds_in_coord_matrix = zeros(1, n_free_coords);
    counter_free_coords = 1;
    for i = 1:n_coords
        data_col = ArrayDouble();
        coordinate_sto.getDataColumn(coords.get(i-1).getName(), data_col);
        q(:,i) = str2num(data_col);
        coord_names_q{i} = char(coords.get(i-1).getName()) ;
        
        % deg to rad for rotational dofs
        if coordinate_sto.isInDegrees()
            if strcmp(coords.get(i-1).getMotionType(), 'Rotational')
                q(:,i) = deg2rad(q(:,i)); 
            end
        end
        
        % indicies (1-indexing) of coord matrix that correspond to free coords
        if sum(strcmp(char(coords.get(i-1).getName()), free_coord_names))
            free_coord_inds_in_coord_matrix(counter_free_coords) = i;
            counter_free_coords = counter_free_coords + 1;
        end
    end

    % filter IK results with 4th order, zero-lag, butterworth lowpass (6Hz)
    cutoff_freq = 6; % Hz
    [b,a] = butter(2, cutoff_freq/(sample_freq/2)); % (filtfilt doubles order)
    q = filtfilt(b,a,q) ;

    qd = zeros(size(q));
    qdd = qd;
    for i = 1:n_coords
        qd(:,i) = gradient(q(:,i), 1/sample_freq);
        qdd(:,i) = gradient(qd(:,i), 1/sample_freq);
    end
    
    %% get actuation forces (inverse dynamics moments and residuals)
    actuation_sto = Storage(actuation_force_filepath);
    orig_actuation_labels = actuation_sto.getColumnLabels();
    
    time_actuation = ArrayDouble();
    actuation_sto.getTimeColumn(time_actuation);
    time_actuation = str2num(time_actuation) ;
    
    % rename pelvis residual actuators
    n_actuation_labels = orig_actuation_labels.getSize - 1; % don't want time
    actuation_labels = cell(1, n_actuation_labels);
    for i = 1:n_actuation_labels
        actuation_label = char(orig_actuation_labels.get(i));
        if strcmp(actuation_label,'FX')
            actuation_label = 'pelvis_tx_force';
        elseif strcmp(actuation_label,'FY')
            actuation_label = 'pelvis_ty_force';
        elseif strcmp(actuation_label,'FZ')
            actuation_label = 'pelvis_tz_force';
        elseif strcmp(actuation_label,'MX')
            actuation_label = 'pelvis_list_moment';
        elseif strcmp(actuation_label,'MY')
            actuation_label = 'pelvis_rotation_moment';
        elseif strcmp(actuation_label,'MZ')
            actuation_label = 'pelvis_tilt_moment';
        elseif strcmp(actuation_label, 'knee_angle_r_beta_force')
            actuation_label = '';
        elseif strcmp(actuation_label, 'knee_angle_l_beta_force')
            actuation_label = '';
        end
        actuation_labels{i} = actuation_label;
    end
    
    % store values
    coord_names_id = cell(1, n_free_coords);
    moments_id = zeros(length(time_actuation), n_free_coords);
    for i = 1:length(free_coord_names)
        data_col = ArrayDouble();
        id_col_idx = find(contains(actuation_labels, ...
            free_coord_names{i})) - 1; % to OpenSim indexing
        actuation_sto.getDataColumn(id_col_idx, data_col);
        moments_id(:,i) = str2num(data_col);
        coord_names_id{i} = free_coord_names{i};
    end
    
    %% interpolate so kinematics & kinetics align
    start_time = max([time_coords(1) time_actuation(1)]);
    end_time = min([time_coords(end) time_actuation(end)]);
    time_vec = linspace(start_time, end_time, 101);

    q = interp1(time_coords, q, time_vec);
    qd = interp1(time_coords, qd, time_vec);
    moments_id = interp1(time_actuation, moments_id, time_vec);
        
    %% set up analyses
    % ForceReporter
    force_report = ForceReporter(model);
    force_report.includeConstraintForces(1);
    model.addAnalysis(force_report);

    % StateReporter
    state_report = StatesReporter(model);
    state_report.setInDegrees(true);
    model.addAnalysis(state_report);

    % JointReaction
    in_frames = ArrayStr;
    on_bodies = ArrayStr;
    joint_names = ArrayStr;
    for i = 1:length(jrxn_expression_frames)
        in_frames.set(i-1, jrxn_expression_frames{i});
        on_bodies.set(i-1, jrxn_bodies{i});
        joint_names.set(i-1, jrxn_joints{i});
    end
    jointRxn = JointReaction();
    jointRxn.setName('JointRxn');
    jointRxn.setInFrame(in_frames);
    jointRxn.setOnBody(on_bodies);
    jointRxn.setJointNames(joint_names);
    jointRxn.setStartTime(start_time);
    model.addAnalysis(jointRxn);
    jointRxn.setModel(model);
    
    %% optimization
    % Find start and end rows
    [~,start_row] = min(abs(time_vec-start_time));
    [~,end_row] = min(abs(end_time-time_vec));
    n_time_steps = end_row - start_row + 1 ; % # of iterations through time loop
    
    % parameters that don't change with each timestep
    
    coeffs_initial = 0.01*ones(n_actuators, 1);
    coeffs_last_step = coeffs_initial ; % init for tendon force estimation

    % linear equality and inequality constraints
    A = [];
    Aeq = [];
    b = [];
    beq = [];
    
    % boundary constraints
    lb = zeros(1, n_free_coords);
    ub = lb; % just initialization
    for i = 1:n_muscles+n_free_coords
        if strcmp(forces.get(i-1).getConcreteClassName(), 'CoordinateActuator')
            lb(i) = -inf;
            ub(i) = inf; 
        elseif contains(char(forces.get(i-1).getConcreteClassName()), 'Muscle')
            lb(i) = 0;
            ub(i) = 1;            
        else
            warning([char(forces.get(i-1).getConcreteClassName()) ...
                'is not a muscle or coordinate actuator'])
            lb(i) = -inf;
            ub(i) = inf; 
        end
    end
    
    % set optimizer options
    options_sqp = optimoptions('fmincon', 'Display', 'off', ...
         'TolCon', 1e-4, 'TolFun', 1e-12, 'TolX', 1e-8, ...
         'MaxFunEvals', 20000, 'MaxIter', 5000, 'Algorithm', 'sqp');
    options_ip = optimoptions('fmincon', 'Display', 'off', ...
         'TolCon', 1e-4, 'TolFun', 1e-12, 'TolX', 1e-8, ...
         'MaxFunEvals', 1000000, 'MaxIter', 10000, 'Algorithm', ...
         'interior-point');
    
    n_sample_steps = 100;
    state_vector = zeros(n_states, 1);
    
    % unchanging vars to pass to optimizer
    params.coords = coords;
    params.actuators = actuators;
    params.muscles = muscles;
    params.n_muscles = n_muscles;
    params.free_coord_names = free_coord_names;
    params.n_free_coords = n_free_coords;
    
    %% time-stepping Optimization Loop
    timesteps = linspace(1, n_time_steps, n_sample_steps);
    timesteps = round(timesteps);
    for j = 1:n_sample_steps
        t_ind = timesteps(j); % Matlab indexing
        row_ind = start_row - 1 + t_ind; % Matlab indexing
        t = time_vec(row_ind);
        fprintf('\tOptimizing t=%.3fs\n',t)

        state.setTime(t) ; % clears all position & velocity info in state

        % set coordinate qs (even locked and constrained coords)
        for i = 1:n_coords
            coords.get(coord_names_q{i}).setValue(state, q(row_ind, i));
            coords.get(coord_names_q{i}).setSpeedValue(state, qd(row_ind, i));
        end
        model.realizeVelocity(state);

        % set state vector values
        for i = 1:n_states
            state_vector(i) = model.getStateVariableValue(state, ...
                states.get(i-1));
        end

        % changing variables to pass to optimizer
        params.model = model;
        params.state = state;
        params.moments_id = moments_id(row_ind,:);

        % compute fixed muscle params (moment arms, active force
        % multipliers, etc.)
        params.muscle_params = get_muscle_params(params, coeffs_last_step);

        % nonlinear constraint function
        nonlcon = @(coeffs0) dynamics_constraint(coeffs0,params) ;

        % call constrained optimization
        try % try SQP optimizer first
            [coeffs_final, fval, exit_flag, ~] = fmincon(@(coeffs0) ...
                cost_function(coeffs0), coeffs_initial, A, b, Aeq, beq, ...
                lb, ub, nonlcon, options_sqp);
            % look for spikes in cost function - use IP if identified
            if t_ind > 1 && fval > 2*fval_last
                disp(['Cost Value was ' num2str(ceil(fval/fval_last)) ...
                      'x last time, will try to use IP'])
                error('Cost Too High')
            elseif exit_flag < 1
                disp(['SQP optimizer had an issue - output state was ' ...
                      num2str(exit_flag)])
                error('Optimizer Failed')
            end
        catch % use IP if having issues with SQP
            disp('Trying IP instead') ;
            [coeffs_final, fval, exit_flag, ~] = fmincon(@(coeffs0) ...
                cost_function(coeffs0), coeffs_initial, A, b, Aeq, beq, ...
                lb, ub, nonlcon, options_ip);
        end
        fval_last = fval ;

        %% Post-hoc Analyses
        % set prescribed controller constant value to control value.
        % Controls don't live through joint reaction. A new state is
        % declared, so the cache is deleted, and this is a workaround.
        for i = 1:n_controllers
             this_controller = PrescribedController...
                 .safeDownCast(controllers.get(i-1));
             this_constant_function = Constant.safeDownCast(this_controller...
                 .get_ControlFunctions(0).get(0));
             this_constant_function.setValue(coeffs_final(n_muscles+i));
        end 

        % set controls
        controls = Vector(n_actuators,0);
        for i = n_muscles:n_actuators-1
            controls.set(i, coeffs_final(i+1));
        end
        model.setControls(state,controls);

        % set muscle activations
        for i = 1:n_muscles
            muscles.get(i-1).setActivation(state,coeffs_final(i))
        end

        model.equilibrateMuscles(state);
        model.realizeAcceleration(state);

        % run analyses
        if t_ind == 1
            force_report.begin(state) ;
            state_report.begin(state) ;
            jointRxn.begin(state) ;
        else
            force_report.step(state, t_ind-1) ;
            state_report.step(state, t_ind-1) ;
            jointRxn.step(state, t_ind-1) ;
        end

        % store activations
        coeffs_last_step = coeffs_final ; % used to compute tendon force
        if num2str(exit_flag) <= 0
            error('Optimization failed')
        end
    end % t_ind - this is the time loop
    
    %% save analyses
    force_report.getForceStorage.print([output_filepath 'results_forces.sto']);
    state_report.getStatesStorage.print([output_filepath 'results_states.sto']);
    jointRxn.end(state);
    jointRxn.printResults('results', output_filepath, -1, '.sto');

    diary off
end

%% helpers
function muscle_params = get_muscle_params(params, coeffs_initial)
    % extract muscle parameters from the state.
    %
    % Args:
    %   params (struct): contains model and state info (muscles, n_muscles,
    %       coords, free_coord_names, n_free_coords, actuators, state,
    %       model).
    %   coeffs_initial (arr): initial coefficients for optimization.
    %
    % Returns:
    %   muscle_params (struct): muscle parameters (passive forces, active
    %       force multipliers, cos(pennation angle), moment arms, optimal
    %       force for coordinate actuators)
    
    import org.opensim.modeling.*
    
    muscles = params.muscles;
    n_muscles = params.n_muscles;
    coords = params.coords;
    free_coord_names = params.free_coord_names;
    n_free_coords = params.n_free_coords;
    actuators = params.actuators;
    state = params.state;
    model = params.model;

    % initialize
    passive_force = zeros(1,n_muscles);
    active_force_multiplier = passive_force;
    cos_alpha = passive_force;
    moment_arms = zeros(n_muscles, n_free_coords);
    coord_opt_force = zeros(1, n_free_coords);
    
    for i = 1:n_muscles
        muscles.get(i-1).setActivation(state, coeffs_initial(i));
    end
    model.equilibrateMuscles(state);

    for i = 1:n_muscles
        passive_force(i) = max([muscles.get(i-1).getPassiveFiberForce(state),...
            0]);
        active_force_multiplier(i) = muscles.get(i-1)...
            .getActiveForceLengthMultiplier(state)...
            * muscles.get(i-1).getForceVelocityMultiplier(state)...
            * muscles.get(i-1).getMaxIsometricForce;
        cos_alpha(i) = muscles.get(i-1).getCosPennationAngle(state);
        for j = 1:n_free_coords
            moment_arms(i,j) = muscles.get(i-1).computeMomentArm(state, ...
                coords.get(free_coord_names{j})) ;
        end
    end

    for i = 1:n_free_coords
        this_coord_actuator = CoordinateActuator.safeDownCast(...
            actuators.get(n_muscles-1+i));
        coord_opt_force(i) = this_coord_actuator.get_optimal_force;
    end
    
    muscle_params.passive_force = passive_force;
    muscle_params.active_force_multiplier = active_force_multiplier;
    muscle_params.cos_alpha = cos_alpha;
    muscle_params.moment_arms = moment_arms;
    muscle_params.coord_opt_force = coord_opt_force;
end

function [c, ceq] = dynamics_constraint(coeffs, params)
    % establish the dynamics constraint for the optimization.
    %
    % Args:
    %   coeffs (struct): coefficients to optimize.
    %   params (arr): contains model and constraint info (n_muscles,
    %   n_free_coords, ID moments, muscle parameters).
    
    n_muscles = params.n_muscles;
    n_free_coords = params.n_free_coords;
    moments_id = params.moments_id;
    MP = params.muscle_params;
    
    % compute model moments and match with ID moments
    moments_sim_muscles = ((coeffs(1:n_muscles)' ...
        .* MP.active_force_multiplier + MP.passive_force) ...
        .* MP.cos_alpha) * MP.moment_arms;
    moments_sim_actuators = coeffs(n_muscles+1:n_muscles+n_free_coords)' ...
        .* MP.coord_opt_force;
    moments_sim = moments_sim_muscles + moments_sim_actuators;

    % compute differences in accelerations
    ceq = moments_sim - moments_id;
    c = [];
end

function c = cost_function(coeffs)
    % establish the cost function for the optimization.
    %
    % Args:
    %   coeffs (struct): coefficients to optimize.
    %
    % Returns:
    %   c (double): cost value to minimize.
    
    c = sum(coeffs.^2);
end    