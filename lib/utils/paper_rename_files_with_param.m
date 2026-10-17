% PAPER_RENAME_FILES_WITH_PARAM suffixes result files with a parameter value

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
function paper_rename_files_with_param(rootFolder, expressionToAvoid, name, value)

    wordToReplace = '.csv';
    replacementWord = ['__' name '_' value '.csv'];

    % Recursively find all files (folder + subfolders)
    listing = dir(fullfile(rootFolder, '**', '*'));
    listing = listing(~[listing.isdir]);   % keep files only

    for k = 1:numel(listing)
        currentName = listing(k).name;

        % If the NAME contains the expression to avoid -> skip
        if contains(currentName, expressionToAvoid)
            continue;
        end

        % If the word to replace is not in the name -> nothing to do
        if ~contains(currentName, wordToReplace)
            continue;
        end

        % New name
        newName = strrep(currentName, wordToReplace, replacementWord);

        oldPath = fullfile(listing(k).folder, currentName);
        newPath = fullfile(listing(k).folder, newName);

        % Rename
        [ok, msg] = movefile(oldPath, newPath);
        if ok
            fprintf('Renamed: %s -> %s\n', currentName, newName);
        else
            warning('Failed to rename %s: %s', oldPath, msg);
        end
    end
end % function

