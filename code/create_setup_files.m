% create_setup_files.m
% copy template setup files and make them subject-specific.
%
% Details:
%   copy the template setup files
%   replace subject-specific info
%   This function does not update COM location in RRA and SO actuator
%       files. Updates are done in make_simulation.m after scaling step is
%       complete.
% 
% Args:
%   subject_datapath (str): directory containing raw datafiles.
%
% Returns:
%   setup_path (str): directory containing simulation setup files and
%       RRA results.
%   walkIDs (cell array): trial identifiers, e.g., {'walk01', 'walk02'}.
%   fp2_sides (cell array): stance limbs per trial, e.g., {'l', 'r'}.
%   start_times (cell array): stride start times.
%   end_times (cell array): stride end times.

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
function [setup_path, walkIDs, fp2_sides, start_times, end_times] = create_setup_files(subject_datapath)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    datapath = '..\data';
    template_setup_path = '..\simulation\Template_SetupFiles';
    
    %% get filepath if not supplied
    if nargin == 0
        subject_datapath = uigetdir(datapath);
    end
    
    % get subject ID
    subject_datapath_split = strsplit(subject_datapath, '\');
    subjectID = subject_datapath_split{end};
    
    %% copy template setup files
    setup_path = [cd '\..\simulation\' subjectID '_SetupFiles'];
    
    if exist(setup_path, 'dir') == 0
        fprintf('Creating setup files for subject %s...\n', subjectID)
        status = copyfile(template_setup_path, setup_path);
        if status == 0
            error('Files not copied')
        end
    end
    
    %% get subject-specific info
    [walkIDs, start_times, end_times, fp1_sides, fp2_sides, ...
        fp3_sides] = read_notes(subject_datapath);
    mass = get_mass(subject_datapath);

    %% replace scale template setup file
    template_name = fullfile(setup_path,...
        'Scale/SubjectID_Setup_ScaleTool.xml');
    new_filename = fullfile(setup_path, 'Scale',...
        [subjectID '_Setup_ScaleTool.xml']);
    if exist(template_name, 'file')
        % replace 'SubjectID' and mass defaults with proper values
        template_to_specific(template_name, new_filename,...
            {'SubjectID', '<mass>-1</mass>'},...
            {subjectID, ['<mass>' mass '</mass>']})
        delete(template_name)
        disp('  Template scale setup file overwritten with subject-specific.')
    end

    %% IK files
    template_name = fullfile(setup_path, ...
        'IK/SubjectID_walkID_setupIK.xml');
    if exist(template_name, 'file')
        for i=1:size(walkIDs, 1)
            walkID = lower(walkIDs{i});
            new_filename = fullfile(setup_path, 'IK', ...
                [subjectID '_' walkID '_setupIK.xml']);
            template_to_specific(template_name, new_filename,...
                {'SubjectID', 'walkID', ...
                    '<time_range> 4.0 5.0 </time_range>'}, ...
                {subjectID, walkID, ['<time_range>' start_times{i} ...
                    ' ' end_times{i} '</time_range>']});
        end
        delete(template_name)
        disp('  Template IK setup file overwritten with subject-specific.')
    end

    %% ID files
    template_name = fullfile(setup_path, ...
        'ID/SubjectID_walkID_setupID.xml');
    if exist(template_name, 'file')
        for i=1:size(walkIDs, 1)
            walkID = lower(walkIDs{i});
            % setup file
            new_filename = fullfile(setup_path, 'ID',...
                [subjectID '_' walkID '_setupID.xml']);
            template_to_specific(template_name, new_filename,...
                {'SubjectID', 'walkID', ...
                    '<time_range> 4.0 5.0 </time_range>'},...
                {subjectID, walkID, ['<time_range>' start_times{i} ...
                    ' ' end_times{i} '</time_range>']});
            % grf file
            new_filename = fullfile(setup_path, 'ID', ...
                ['grf_' walkID '.xml']);
            grfTemplateName = fullfile(setup_path, 'ID/grf_walkID.xml');
            template_to_specific(grfTemplateName, new_filename,...
                {'SubjectID', 'walkID', 'calcn_lr1', 'calcn_lr2', ...
                    'calcn_lr3'}, ...
                {subjectID, walkIDs{i}, ...
                    ['calcn_' fp1_sides{i}], ['calcn_' fp2_sides{i}], ...
                    ['calcn_' fp3_sides{i}]});
        end
        delete(template_name)
        delete(grfTemplateName)
        disp('  Template ID setup file overwritten with subject-specific.')
    end
    
    %% RRA files
    template_name = fullfile(setup_path, ...
        'RRA/SubjectID_walkID_setupRRA.xml');
    if exist(template_name, 'file')
        for i=1:size(walkIDs, 1)
            walkID = lower(walkIDs{i});
            % setup file
            new_filename = fullfile(setup_path, 'RRA',...
                [subjectID '_' walkID '_setupRRA.xml']);
            template_to_specific(template_name, new_filename,...
                {'SubjectID', 'walkID', ...
                    '<initial_time>4.0</initial_time>',...
                    '<final_time>5.0</final_time>'},...
                {subjectID, walkID,...
                    ['<initial_time>' start_times{i} '</initial_time>'], ...
                    ['<final_time>' end_times{i} '</final_time>']});
            % grf file
            new_filename = fullfile(setup_path, 'RRA',...
                ['grf_' walkID '.xml']);
            grfTemplateName = fullfile(setup_path, 'RRA/grf_walkID.xml');
            template_to_specific(grfTemplateName, new_filename,...
                {'SubjectID', 'walkID', 'calcn_lr1', 'calcn_lr2', ...
                    'calcn_lr3'}, ...
                {subjectID, walkIDs{i}, ...
                    ['calcn_' fp1_sides{i}], ['calcn_' fp2_sides{i}], ...
                    ['calcn_' fp3_sides{i}]});
            % actuators and tasks already copied in
        end
        delete(template_name)
        delete(grfTemplateName)
        disp('  Template RRA setup file overwritten with subject-specific.')
    end
    
    %% Analysis files
    template_name = fullfile(setup_path, ...
        'Analysis/PreRRA/SubjectID_walkID_setupAnalysis.xml');
    if exist(template_name, 'file')
        for i=1:size(walkIDs, 1)
            walkID = lower(walkIDs{i});
            % setup file
            new_filename = fullfile(setup_path, 'Analysis/PreRRA',...
                [subjectID '_' walkID '_setupAnalysis.xml']);
            template_to_specific(template_name, new_filename,...
                {'SubjectID', 'walkID', ...
                    '<initial_time>4.0</initial_time>', ...
                    '<final_time>5.0</final_time>', ...
                    '<start_time>4.0</start_time>', ...
                    '<end_time>5.0</end_time>'}, ...
                {subjectID, walkID, ...
                    ['<initial_time>' start_times{i} '</initial_time>'], ...
                    ['<final_time>' end_times{i} '</final_time>'], ...
                    ['<start_time>' start_times{i} '</start_time>'], ...
                    ['<end_time>' end_times{i} '</end_time>']});
        end
        delete(template_name)
        disp('  Template Pre-RRA Analysis setup file overwritten with subject-specific.')
    end
    
    template_name = fullfile(setup_path, ...
        'Analysis/PostRRA/SubjectID_walkID_setupAnalysis.xml');
    if exist(template_name, 'file')
        for i=1:size(walkIDs, 1)
            walkID = lower(walkIDs{i});
            % setup file
            new_filename = fullfile(setup_path, 'Analysis/PostRRA',...
                [subjectID '_' walkID '_setupAnalysis.xml']);
            template_to_specific(template_name, new_filename,...
                {'SubjectID', 'walkID', ...
                    '<initial_time>4.0</initial_time>', ...
                    '<final_time>5.0</final_time>', ...
                    '<start_time>4.0</start_time>', ...
                    '<end_time>5.0</end_time>'}, ...
                {subjectID, walkID, ...
                    ['<initial_time>' start_times{i} '</initial_time>'], ...
                    ['<final_time>' end_times{i} '</final_time>'], ...
                    ['<start_time>' start_times{i} '</start_time>'], ...
                    ['<end_time>' end_times{i} '</end_time>']});
        end
        delete(template_name)
        disp('  Template Post-RRA Analysis setup file overwritten with subject-specific.')
    end
    
    
    %% No files needed for SO (all generated from scratch from script)

end

function [walkIDs, start_times, end_times, fp1_sides, fp2_sides, ...
    fp3_sides] = read_notes(subject_datapath)
    % read notes file for in subject_datapath and return info.
    %
    % Args:
    %   subject_datapath (str): directory containing raw datafiles.
    %
    % Returns:
    %   walkIDs (cell array): trial identifiers, e.g., {'walk01', 'walk02'}.
    %   start_times (cell array): stride start times.
    %   end_times (cell array): stride end times.
    %   fp1_sides (cell array): laterality of force plate 1 limbs per
    %       trial, e.g., {'r', 'l'}.
    %   fp2_sides (cell array): laterality of force plate 2 limbs per
    %       trial, e.g., {'l', 'r'}.
    %   fp3_sides (cell array): laterality of force plate 3 limbs per
    %       trial, e.g., {'r', 'l'}.
    
    subject_datapath_split = strsplit(subject_datapath, '\');
    subjectID = subject_datapath_split{end};
    filename = [subject_datapath '\' subjectID ' notes.xlsx'];
    T = readtable(filename);
    
    % get walk IDs 
    walkIDs = T{:,1};
    
    % get start/end times and force plate laterality for each walk
    n_walks = height(T);
    start_times = cell(1, n_walks);
    end_times = cell(1, n_walks);
    fp1_sides = cell(1, n_walks);
    fp2_sides = cell(1, n_walks);
    fp3_sides = cell(1, n_walks);
    for i=1:n_walks
        start_times{1,i} = char(string(T{i,3}));
        end_times{1,i} = char(string(T{i,4}));
        step_order = char(T{i,2}); % L3 R2 L1, L1 R2 L3, etc.
        
        fp1_side = lower(step_order(strfind(step_order, '1')-1)); % 'l'/'r'
        fp2_side = lower(step_order(strfind(step_order, '2')-1));
        fp3_side = lower(step_order(strfind(step_order, '3')-1));
        
        fp1_sides{1,i} = fp1_side;
        fp2_sides{1,i} = fp2_side;
        fp3_sides{1,i} = fp3_side;
            
    end
end

function mass = get_mass(subject_datapath)
    % read demographics file to get mass.
    %
    % Args:
    %   subject_datapath (str): directory containing raw datafiles.
    %
    % Returns:
    %   mass (double): mass of subject in kg.

    subject_datapath_split = strsplit(subject_datapath, '\');
    subjectID = subject_datapath_split{end};
    filename = [subject_datapath '\..\demographics.xlsx'];
    T = readtable(filename, 'ReadVariableNames', false);
    
    weight_col = 4;
    mass = char(T{strcmp(T.Var1, subjectID), weight_col}); % kg
end

function template_to_specific(template_filename, new_filename, ...
    old_strings, new_strings)
    % read template file, search and replace, and write new file.
    %
    % Args:
    %   template_filename (str): path of template file.
    %   new_filename (str): path of file to write.
    %   old_strings (cell array): strings to replace.
    %   new_strings (cell array): replacement strings.
    %
    % Returns:
    %   none.
    
    fid = fopen(template_filename, 'r');
    f = fread(fid, '*char')';
    fclose(fid);

    for i=1:size(old_strings, 2)
        f = regexprep(f, old_strings{i}, new_strings{i});
    end
    
    fid = fopen(new_filename, 'w');
    fprintf(fid, '%s', f);
    fclose(fid);
end
