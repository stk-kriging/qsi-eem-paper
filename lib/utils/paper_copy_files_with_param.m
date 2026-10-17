% PAPER_COPY_FILES_WITH_PARAM Recursively scans a folder and copies the files
% suffixed with a given parameter/value pair to another pair.
%
%   paper_copy_files_with_param(rootFolder, fromName, fromValue, toName, toValue)
%
% Inputs:
%   rootFolder             - path of the starting folder
%   fromName, fromValue    - suffix of the files to copy
%   toName, toValue        - suffix of the copies
%
% Suffixes are built as: ['__' name '_' value '.csv']

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
function paper_copy_files_with_param(rootFolder, fromName, fromValue, toName, toValue)

    fromWord = ['__' fromName '_' fromValue '.csv'];
    toWord = ['__' toName '_' toValue '.csv'];

    % Recursively find all files (folder + subfolders)
    listing = dir(fullfile(rootFolder, '**', '*'));
    listing = listing(~[listing.isdir]);   % keep files only

    for k = 1:numel(listing)
        currentName = listing(k).name;

        % If the suffix to copy is not in the name -> nothing to do
        if ~contains(currentName, fromWord)
            continue;
        end

        % New name
        newName = strrep(currentName, fromWord, toWord);

        oldPath = fullfile(listing(k).folder, currentName);
        newPath = fullfile(listing(k).folder, newName);

        % Copy
        [ok, msg] = copyfile(oldPath, newPath);
        if ok
            fprintf('Copied: %s -> %s\n', currentName, newName);
        else
            warning('Failed to copy %s: %s', oldPath, msg);
        end
    end
end % function
