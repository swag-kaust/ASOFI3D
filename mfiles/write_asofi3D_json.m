function write_asofi3D_json(filename, config) 
% writes to json file

%% make all numbers strings 
field_list = fieldnames(config);

for field_n = 1:length(field_list)
	config.(field_list{field_n}) = ...
        num2str(config.(field_list{field_n}));
end

%% encode to json

bb = jsonencode(config);

% replace manually the field names altered by MATLAB
bb = replace(bb, '"REFRECX_REFRECY_REFRECZ"','"REFRECX, REFRECY, REFRECZ"');
bb = replace(bb,'"XREC1_YREC1_ZREC1"','"XREC1,YREC1, ZREC1"');
bb = replace(bb,'"XREC2_YREC2_ZREC2"','"XREC2,YREC2, ZREC2"');

bb = replace(bb, '"SOURCE_ALPHA_SOURCE_BETA"', '"SOURCE_ALPHA, SOURCE_BETA"');
bb = replace(bb, '"AMON_STR_DIP_RAKE"', '"AMON, STR, DIP, RAKE"');
bb = replace(bb, '"AMON_M11_M12_M13_M22_M23_M33"',...
    '"AMON, M11, M12, M13, M22, M23, M33"');

bb = replace(bb, '","', ['",', newline,'"']);
bb = sprintf(bb);

bb = replace(bb, '"NDT_NDTSHIFT"', '"NDT, NDTSHIFT"');
filewrite(filename,bb);

end

function filewrite(file_name, text)
    fid = fopen(file_name,'wt');
    fprintf(fid, '%s', text);
    fclose(fid);
end
