% WRITEMATRIX_ writematrix, with a dlmwrite fallback for Octave

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
function writematrix_ (data, filename)

try

    % This is the recommended way in Matlab since R2019
    writematrix (data, filename);

catch err

    if exist ('OCTAVE_VERSION', 'builtin') == 5
        undef = strcmp (err.identifier, 'Octave:undefined-function');
    else
        undef = strcmp (err.identifier, 'MATLAB:UndefinedFunction');
    end

    if undef
        % This should work if writematrix is not available
        dlmwrite (filename, data, 'delimiter', ',', 'precision', 16); %#ok<DLMWT>
    else
        rethrow (err);
    end

end % try-catch

end % function

