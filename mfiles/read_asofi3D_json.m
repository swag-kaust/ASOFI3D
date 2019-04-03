function config = read_asofi3D_json(filename)
%READ_ASOFI3D_JSON  Read configuration file into `struct`.
%   json_config = read_asofi3D_json('in_and_out/sofi3D.json') reads file
%   'in_and_out/sofi3D.json' relative to the current directory.

json_text = fileread(filename);
i = find(json_text=='{');
j = find(json_text=='}');
b = json_text(i:j);
config = jsondecode(b);





% ASOFI3D JSON parameter file contains all variables in string datatype.
% We convert numerical parameters to double datatype, this does not over
% parsing

field_list = fieldnames(config);

for field_n = 1:length(field_list)
    config.(field_list{field_n}) = ...
        str2num_if_num(config.(field_list{field_n}));
end
%write_asofi3D_json([filename, '_snap'], config);

end

function numval = str2num_if_num(strval)
%makes double from string only if it is a number

numval = str2double(strval);

if isnan(numval)
    numval = strval;
end
end
