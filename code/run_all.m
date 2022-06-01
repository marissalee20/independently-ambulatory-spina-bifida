%% run_all.m
% using raw datafiles, runs entire simulation pipeline for the selected
% subject.
%
% Details:
%   filter GRF files
%   create setup files
%   run simulations for each walk
% 
% Args:
%   subject_datapath (str): directory containing raw datafiles. If not
%       supplied, user will be prompted to select a datapath.
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
function run_all(subject_datapath)
    % add code directory to path
    addpath(cd);
    
    % prompt user to choose datapath if missing
    if nargin == 0
        subject_datapath = uigetdir('..\data');
    end
    
    filter_grfs(subject_datapath);
    [sim_home_dir, walkIDs, walk_sides, start_times, end_times] = ...
        create_setup_files(subject_datapath);

    %% make simulation for each walk
    for i=1:size(walkIDs, 1)
        walkID = lower(char(walkIDs{i}));
        walk_side = char(walk_sides{i});
        start_time = str2double(start_times{i});
        end_time = str2double(end_times{i});
        make_simulation(sim_home_dir, walkID, walk_side, start_time, ...
            end_time);
    end    
end