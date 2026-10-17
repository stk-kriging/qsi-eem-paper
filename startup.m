% STARTUP path setup and dependency bootstrap of qsi-eem-paper

% Copyright Notice
%
%    Copyright (C) 2024, 2026 CentraleSupelec
%
%    Author(s):  Julien Bect  <julien.bect@centralesupelec.fr>
%                Romain Ait Abdelmalek-Lomenech <romain.ait@centralesupelec.fr>

% Copying Permission Statement
%
%    This file is part of qsi-eem-paper
%       (https://github.com/stk-kriging/qsi-eem-paper)
%
%    qsi-eem-paper is free software: you can redistribute it and/or modify it under
%    the terms of the GNU General Public License as published by the Free
%    Software Foundation,  either version 3  of the License, or  (at your
%    option) any later version.
%
%    qsi-eem-paper is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See  the GNU  General Public
%    License for more details.
%
%    You should  have received a copy  of the GNU  General Public License
%    along with qsi-eem-paper.  If not, see <http://www.gnu.org/licenses/>.

function startup ()

% Locate project root
root = fileparts (mfilename ('fullpath'));
here = pwd ();

% Download BSS (if needed)
contrib_bss = fullfile (root, 'lib', 'algorithms', 'bss', 'contrib-bss');
if ~ exist (contrib_bss, 'dir')
    % commit 86aaee3 / date: Thu Jul 16 08:29:31 2026 +0200
    git_clone_dependency ('contrib-bss', contrib_bss, ...
        'https://github.com/stk-kriging/contrib-bss.git', '86aaee3');
end

% Download contrib-qsi (if needed)
contrib_qsi = fullfile (root, 'lib', 'algorithms', 'qsi');
if ~ exist (contrib_qsi, 'dir')
    % commit a217469 / date: Fri Jul 24 20:05:35 2026 +0200
    git_clone_dependency ('contrib-qsi', contrib_qsi, ...
        'https://github.com/stk-kriging/contrib-qsi.git', 'a217469');
end


%% Set path for this repo

lib_dir = fullfile (root, 'lib');
addpath (lib_dir);
addpath (fullfile (lib_dir, 'utils'));

algo_dir = fullfile (lib_dir, 'algorithms');
addpath (fullfile (algo_dir, 'bss'));
addpath (genpath (fullfile (algo_dir, 'qsi_eem')));

addpath (fullfile (root, 'demo'));
addpath (genpath (fullfile (root, 'test_cases')));

experiments = fullfile (root, 'experiments');
addpath (fullfile (experiments, 'eem_sur_doe_crit_comparison'));
addpath (fullfile (experiments, 'quadrature_size_comparison'));
addpath (fullfile (experiments, 'qsi_comparison'));
addpath (fullfile (experiments, 'small_rotor37'));
addpath (fullfile (experiments, 'small_synthetic'));
addpath (fullfile (experiments, 'eem_smc_parameters_comparison'));
addpath (fullfile (experiments, 'time_eem_sur'));


%% Set path for dependencies

% Set contrib-qsi path (& initialize STK)
run (fullfile (contrib_qsi, 'startup.m'));

% Set BSS path
addpath (genpath (fullfile (contrib_bss, 'functions')));

% Avoid concurrent compilation when startup runs on pool workers.
bss_mex_dir = fullfile (contrib_bss, 'functions', 'bss');
bss_mex_src = fullfile (bss_mex_dir, 'replicate_particles.c');
bss_mex_ext = mexext ();
bss_mex_bin = fullfile (bss_mex_dir, ['replicate_particles.' bss_mex_ext]);
src_info = dir (bss_mex_src);
bin_info = dir (bss_mex_bin);

if isempty (bin_info) || (src_info.datenum > bin_info.datenum)
    try
        cd (bss_mex_dir);
        mex replicate_particles.c
        cd (here);
    catch e
        cd (here);
        rethrow (e);
    end % try-catch
end


%% Load statistics package if using Octave

if exist ('OCTAVE_VERSION', 'builtin') == 5
    pkg load statistics
end

end % function


function git_clone_dependency (name, dst, url, sha1_or_tag)

if exist (dst, 'dir')
    error (sprintf ('Directory already exists: %s\n', dst)); %#ok<SPERR>
end

fprintf ('Cloning %s... ', name);

here = pwd ();

try

    gitcmd = create_git_call (sprintf ('clone %s %s', url, dst));
    evalc (sprintf ('[status, output] = system (''%s'')', gitcmd));
    if status ~= 0
        error ([ ...
            'git-clone failed with the following ' ...
            'error message:\n\n%s\n\n'], output);
    end

    cd (dst);

    gitcmd = create_git_call (sprintf ('checkout %s', sha1_or_tag));
    evalc (sprintf ('[status, output] = system (''%s'')', gitcmd));
    if status ~= 0
        error ([ ...
            'git-checkout failed with the following ' ...
            'error message:\n\n%s\n\n'], output);
    end

    cd (here);

catch e

    cd (here);

    % Remove partial/failed install
    if exist (dst, 'dir')
        rmdir (dst, 's');
    end

    rethrow (e);

end % try-catch

fprintf ('OK\n');

end % function


function gitcmd = create_git_call (cmd)

if isunix ()
    prefix = 'env -u LD_LIBRARY_PATH ';
else
    prefix = '';
end

gitcmd = sprintf ('%s git %s', prefix, cmd);

end % function
