% PAPER_RUN_TIMED_ALGORITHM runs an algorithm, times it, logs the outcome

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
function wall_time = paper_run_timed_algorithm (label,run_func,time_file)

run_timer = tic();
captured_output = '';
try
  captured_output = evalc('feval(run_func);');
  wall_time = toc(run_timer);
catch run_error
  wall_time = toc(run_timer);
  fprintf('  %-20s failed in %.2f s\n',label,wall_time);
  if ~isempty(captured_output)
    fprintf('%s',captured_output);
  end
  rethrow(run_error);
end

fprintf('  %-20s done in %.2f s',label,wall_time);

if nargin >= 3 && ~isempty(time_file) && exist(time_file,'file')
  step_time = readmatrix_(time_file);
  step_time = step_time(:);
  step_time = step_time(isfinite(step_time));
  if ~isempty(step_time)
    step_word = 'steps';
    if numel(step_time) == 1
      step_word = 'step';
    end
    fprintf(' | %d %s, step sum %.2f s, mean %.2f s', ...
      numel(step_time),step_word,sum(step_time),mean(step_time));
  end
end

fprintf('\n');

end % function

