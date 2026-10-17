% ROTOR37_CHECK_AVAILABLE errors out if the proprietary ROTOR37 code is missing

function rotor37_check_available ()

% functions of the proprietary code used by rotor37_struct
required = {'rotor37', 'rotor37_prepare_s', 'rotor37_eval_x'};

missing = {};
for i = 1:(numel (required))
    if exist (required{i}) ~= 2
        missing{end+1} = required{i};  %#ok<AGROW>
    end
end

if ~ isempty (missing)
    code_dir = fullfile (fileparts (mfilename ('fullpath')), 'code');
    error (['The proprietary ROTOR37 code is not available ' ...
        '(missing: %s).\nIt is not publicly available at the moment; ' ...
        'if you have it, place it in %s\nand run startup again ' ...
        '(see README.md in that directory).'], ...
        strtrim (sprintf ('%s ', missing{:})), code_dir);
end

end % function
