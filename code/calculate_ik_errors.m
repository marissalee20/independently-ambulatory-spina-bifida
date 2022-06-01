%% calculate_ik_errors.m
% calculate average marker RMS and maximum marker errors across all trials.

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

mean_rmses = [];
max_errs = [];
for i = 1:length(subjects)
    subjectID = char(subjects(i));
    file_list = dir(['..\simulation\' subjectID '_SetupFiles\IK\' ...
        '*_out.log']);
    for j = 1:length(file_list)
        logname = [file_list(j).folder '\' file_list(j).name];
        [mean_rmse, max_err] = get_errors(logname);
        mean_rmses = [mean_rmses, mean_rmse];
        max_errs = [max_errs, max_err];
    end
end

fprintf('Average RMSE: %.1f ± %.1f\n', 100*mean(mean_rmses), 100*std(mean_rmses));
fprintf('Max Error: %.1f ± %.1f\n', 100*mean(max_errs), 100*std(max_errs));

%% Helpers
function [mean_rmse, max_err] = get_errors(logname)
    % extract mean marker RMSE and maximum marker error from an inverse
    % kinematics out.log file.
    %
    % Args:
    %   logname (str): filepath of out.log file.
    %
    % Returns:
    %   mean_rmse (float): mean marker RMSE across the trial.
    %   max_err (float): maximum individual marker error across the trial.
    
    fileID = fopen(logname,'r');
    log = fscanf(fileID,'%c');
    fclose(fileID);
    [~, mean_rmse] = read_value(log, 'RMS');
    [max_err, ~] = read_value(log, 'max');
end

function [max_val, mean_val] = read_value(log, string)
    % extract maximum and mean values of a specified string from the
    % out.log file.
    %
    % Args:
    %   log (str): out.log file contents.
    %   string (str): metric of interest.
    %
    % Returns:
    %   max_val (float): maximum metric value.
    %   mean_val (float): mean metric value.
    
    string_offset = length(string)+1;
    string_starts = strfind(log,[string '=']) + string_offset;
    string_ends = string_starts + 6; % 6 gets enough decimal places to evaluate.

    n_timesteps = length(string_starts);
    vals = zeros(1, n_timesteps);
    for i = 1:n_timesteps
        string_start = string_starts(i);
        string_end = string_ends(i);
        vals(1,i) = str2double(log(string_start:string_end));
    end
    
    max_val = max(vals);
    mean_val = mean(vals);
end
