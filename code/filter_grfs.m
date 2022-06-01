%% filter_grfs.m
% filter ground reaction forces for all trials associated with subject,
% using a 4th-order low-pass Butterworth filter, and save results.
%
% Args:
%   subject_datapath (str): directory containing raw datafiles.
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
function filter_grfs(subject_datapath)
    raw_grf_files = dir(fullfile(subject_datapath, '*_grf.mot'));
    filtered_grf_files = dir(fullfile(subject_datapath, ...
        '*_grf_filtered.mot'));
    
    n_raw_grf_files = length(raw_grf_files);
    n_filtered_grf_files = length(filtered_grf_files);
    
    if n_raw_grf_files ~= n_filtered_grf_files
        disp('Filtering raw GRFs...')
        for i = 1:length(raw_grf_files)
            grf_filename = raw_grf_files(i).name;
            full_filename = [raw_grf_files(i).folder '\' grf_filename];
            % if filtered file doesn't exist, make it
            if ~isfile([full_filename(1:end-4) '_filtered.mot'])
                grf_filename = fullfile(subject_datapath, grf_filename);
                fprintf('  %s\n',grf_filename)
                filter_grf(grf_filename);
            end
        end
    end
end


function filter_grf(filename)
    % filter ground reaction forces for a single trial, using a 4th-order
    % low-pass Butterworth filter, and save results.
    %
    % Args:
    %   filename (str): path to .mot file containing grf data.
    %
    % Returns:
    %   none.

    %% Read in data
    [data,~,header] = read_sto(filename);

    %% for all but time column, filter data
    time = data(:,1);
    data_to_filter = data(:,2:37);

    time_step = time(2) - time(1);
    sample_rate = 1 / time_step;

    lowpass_freq = 6;

    filtered_data = lowpass(data_to_filter, sample_rate, lowpass_freq);

    output_data = [time filtered_data];
    n_cols = size(output_data, 2);
    format_spec = ['%.19f\t' repmat('%d\t', [1 n_cols-1])];
    
    %% Write out data
    write_filename = [filename(1:end-4) '_filtered.mot'];
    write_sto(header, output_data', format_spec, write_filename);
    
end

function data = lowpass(data, sample_rate, lowpass_freq)
    % filter data using a 4th-order low-pass Butterworth filter.
    %
    % Args:
    %   data (array): data to filter, with time moving along the rows.
    %   sample_rate (double): sample rate of data array.
    %   lowpass_freq (double): cutoff frequency.
    %
    % Returns:
    %   data (array): filtered data.
    %
    % Notes:
    %   this code is based heavily on code written by Scott Uhlrich.
    
    half_sample_rate = sample_rate / 2;
    cutoff_freq = lowpass_freq / half_sample_rate;
    [b,a] = butter(4, cutoff_freq, 'low');
    % filter with 'filtfilt', which filters in both the forward and reverse
    % directions to eliminate phase distortion
    for i = 1:size(data,2)
        data(:,i) = filtfilt(b, a, data(:,i));
    end
end

function [data, labels, header] = read_sto(filename)
    % read data from a .sto (or .mot) file.
    %
    % Args:
    %   filename (str): path of the file of interest.
    %
    % Returns:
    %   data (array): data from the .sto file.
    %   labels (cell array): labels associated with data.
    %   header (str): full file header, including labels.
    
    % Initialize variables.
    delimiter = '\t';

    % Open the text file.
    fileID = fopen(filename, 'r');

    header_lines = 0;
    end_of_header = false;
    header = '';
    while ~end_of_header && ~feof(fileID)
        thisLine = fgetl(fileID);
        header = [header sprintf(thisLine) '\n'];
        header_lines = header_lines + 1;
        end_of_header = contains(thisLine, 'endheader');
    end
    if feof(fileID)
        fprintf('EOF reached. No ''endheader'' line found.\n')
        data = [];
        labels = [];
        return
    end
    
    % Read line containing variable labels
    labels_line = fgetl(fileID);
    header = [header sprintf(labels_line) '\n'];
    labels = strsplit(labels_line, delimiter);
    nCols = length(labels);

    % Read columns of data according to format string
    format_spec = repmat('%f', [1 nCols]);
    data_array = textscan(fileID, format_spec, Inf, 'Delimiter', delimiter, 'ReturnOnError', false);

    % Close the text file.
    fclose(fileID);

    % Create output data
    data = [data_array{1:end}];
end

function write_sto(header, output_data, format_spec, filename)
    % write data to a .sto (or .mot) file.
    %
    % Args:
    %   header (str): full file header, including labels.
    %   data (array): data to write to the .sto file.
    %   format_spec (str): format specifications for data columns.
    %   filename (str): path of output file.
    %
    % Returns:
    %   none.
    
    fileID = fopen(filename, 'w');
    fprintf(fileID, header);
    fprintf(fileID,[format_spec '\n'], output_data);
    fclose(fileID);
end
