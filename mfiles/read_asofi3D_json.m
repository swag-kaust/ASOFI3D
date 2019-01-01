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
% Convert necessary numerical parameters to double datatype.
field_list = [
    "NX", "NY", "NZ", ...
    "NPROCX", "NPROCY", "NPROCZ", ...
    "IDX", "IDY", "IDZ", ...
    "DX", "DY", "DZ" ...
	"TSNAP1", "TSNAP2", "TSNAPINC", "TIME", ...
];
for field = field_list
	config.(field) = cast_json_param(config.(field));
end

end

function numval = cast_json_param(strval)
%CAST_JSON_PARAM  Convert string `strval` to its numerical value.
if ischar(strval)
    numval = str2double(strval);
else
    numval = strval;
end

end
