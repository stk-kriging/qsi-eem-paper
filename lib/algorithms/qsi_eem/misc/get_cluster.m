% GET_CLUSTER local cluster with a private job storage location

% Copyright Notice
%
%    Copyright (C) 2025, 2026 CentraleSupelec
%
%    Author(s): Romain Ait Abdelmalek-Lomenech <romain.ait@centralesupelec.fr>
%               Julien Bect <julien.bect@centralesupelec.fr>
%               Emmanuel Vazquez <emmanuel.vazquez@centralesupelec.fr>

% Copying Permission Statement
%
%    This file is part of qsi-eem-paper (https://github.com/stk-kriging/qsi-eem-paper)
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
function [local_cluster, cleanUp] = get_cluster (profile)

    if nargin < 1
        profile = 'local';
    end

    iodir_base = fullfile (tempdir (), 'cluster_tmp');
    if ~ exist (iodir_base, 'dir')
        mkdir (iodir_base);
    end
    if ~ exist (iodir_base, 'dir')
        error ('This should never happen...');
    end

    iodir = tempname (iodir_base);
    cleanUp = onCleanup(@()rmdir(iodir, 's')); %delete tempdir at the end.

    if exist (iodir, 'dir')
        error ('This should never happen...');
    end
    mkdir (iodir);
    if ~ exist (iodir, 'dir')
        error ('This should never happen...');
    end

    local_cluster = parcluster(profile);

    %solve bug multiple parallel jobs
    local_cluster.JobStorageLocation = iodir;

end % function

