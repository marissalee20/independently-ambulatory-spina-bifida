%% evaluate_ik.m
% evaluate IK results for a given subject, checking for max and RMS marker
% errors.
%
% Details:
%   loop through all IK results for a subject and notify if max errors are
%   larger than 0.04m or if RMS is larger than 0.02m
% 
% Args:
%   subjectID (str): subject identifier, e.g., 'C1'.
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
function evaluate_ik(subjectID)
    log_name = [cd '\..\simulation\' subjectID '_SetupFiles\IK\out.log'];
    check_errors(log_name);
end

function check_errors(log_name)
    % check max and RMS errors for IK evaluation.
    %
    % Args:
    %   log_name (str): filepath of IK log file.
    %
    % Returns:
    %   none.
    
    fileID = fopen(log_name, 'r');
    log = fscanf(fileID, '%c');
    fclose(fileID);
    read_value(log, 'max', 0.04);
    read_value(log, 'RMS', 0.02);
end

function read_value(log, metric_name, thresh)
    % read values of interest.
    %
    % Args:
    %   log (str): text from which to read values.
    %   metric_name (str): name of metric of interest in log.
    %   thresh (double): error threshold above which to raise an error.
    %
    % Returns:
    %   none.
    
    string_offset = length(metric_name) + 1;
    string_starts = strfind(log, [metric_name '=']) + string_offset;
    string_ends = string_starts + 6; % 6 = enough decimal places to evaluate

    n_timesteps = length(string_starts);
    vals = zeros(1, n_timesteps);
    for i = 1:n_timesteps
        string_start = string_starts(i);
        string_end = string_ends(i);
        vals(1,i) = str2double(log(string_start:string_end));
    end
    
    n_error = sum(vals>=thresh);
    percent_error = n_error/n_timesteps*100;
    max_val = max(vals);
    if n_error > 0
        cnt = questdlg(sprintf...
                (['IK %s error is greater than %.0f cm for %.0f%% of ' ...
                  'the walk. %s val of %.1f cm. Continue?'], metric_name, ...
                 thresh*100, percent_error, metric_name, max_val*100), ...
                'IK Errors', 'Yes', 'No', 'No');
        if strcmp(cnt, 'No')
            error('IK canceled.')
        end
    end
end