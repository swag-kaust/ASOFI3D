function config = read_asofi3D_json(filename)
%READ_ASOFI3D_JSON  Read configuration file into `struct`.
%   json_config = read_asofi3D_json('in_and_out/sofi3D.json') reads file
%   'in_and_out/sofi3D.json' relative to the current directory.

json_text = fileread(filename);
i = find(json_text=='{');
j = find(json_text=='}');
b = json_text(i:j);
config = jsondecode(b);

end