% PAPER_EXTRACT_DEVIATION_VARIABLE_PARAMS grid-based deviation, suffixed runs

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
function paper_extract_deviation_variable_params ...
    (problem, config, name, parameter, value, it, data_dir, mean_plugin)

fprintf('Run %d\n', it);

if nargin < 7
    here = fileparts (mfilename ('fullpath'));
    data_dir = fullfile (here, '..', 'data');
end

if nargin < 8
    mean_plugin = false;
end

[prm, f, s_trnsf] = problem ();
config = config();

PTS_X = config.pts_eval_x;
PTS_S = config.pts_eval_s;

dim_tot = prm.dim_x + prm.dim_s;

xf = stk_sampling_sobol(PTS_X, prm.dim_x, prm.BOXx);
sf = stk_sampling_sobol(PTS_S, prm.dim_s, prm.BOXs);
sf = s_trnsf(sf);
df = adapt_set (xf, sf);

file = sprintf ('results_grid_%s.csv', prm.name);
file = fullfile (data_dir, 'grid', file);
if ~ exist (file, "file")
    zf = f(df);
    csvwrite (file, zf);
else
    zf = csvread(file);
end
zf = double(zf);

trueSet = get_true_quantile_set(zf, PTS_X, PTS_S, prm.alpha, prm.const);

file_design = sprintf('doe_%s_%s_%d__%s_%s.csv', name, prm.name, it, parameter, value);
design = readmatrix_(fullfile(data_dir, 'results/design', file_design));

para = zeros(config.T+1, dim_tot+1, prm.M);
file_cov = zeros(config.T+1,1,prm.M);

for m =1:prm.M
    filename_para = sprintf('param_%s_%d_%s_%d__%s_%s.csv', name, m, prm.name, it, parameter, value);
    para(:,:,m) = readmatrix_(fullfile(data_dir, 'results/param/', filename_para));
    filename_cov = sprintf('cov_%s_%d_%s_%d__%s_%s.csv', name, m, prm.name, it, parameter, value);
    file_cov(:,:,m) = readmatrix_(fullfile(data_dir, 'results/param/', filename_cov));
end

dev = [];

for j = 1:config.axT:config.T+1
    dt = design(1:config.pts_init+j-1,:);
    zt = f(dt);
    Model = stk_model ();

    for m = 1:prm.M
        cov = prm.list_cov{file_cov(j,:,m)};
        Model(m) = stk_model(cov, dim_tot);
        Model(m).param = para(j,:,m);
    end

    if mean_plugin
        zp = stk_predict(Model, dt, zt, df);
        zp = zp.mean;
        approxSet = get_true_quantile_set(zp, PTS_X, PTS_S, prm.alpha, prm.const);
    else
        approxSet = get_expected_quantile_set(Model,df,PTS_X, PTS_S,dt,zt,prm.const,prm.alpha);
    end

    dev = [dev, lebesgue_deviation(trueSet,approxSet)];

end

filename_dev = sprintf('dev_%s_%s_%d__%s_%s.csv', name, prm.name, it, parameter, value);
writematrix_(dev,fullfile(data_dir, 'results/deviations', filename_dev));

end % function

