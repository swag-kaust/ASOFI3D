function config = read_asofi3D_json(filename)

json_text = fileread(filename);
i = find(json_text=='{');
j = find(json_text=='}');
b = json_text(i:j);
config = jsondecode(b);

end