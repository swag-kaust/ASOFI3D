function [tr]=su2matlab(filename);

% input: 
%   filename    :     complete filename of su-file
%
% output:
%   tr          :     struct containing the header information and traces
%                     of the su-file; 
%                     the fields of the struct containing the header information
%                     are named after the corresponding su header words;
%                     the traces itself are saved in the struct field 'trace'



% Determine important parameters of the given file
fileparameters = dir(filename);


% Open the file
fid=fopen(filename,'r','ieee-le');


% Creating an empty struct with the needed fields
tr=create_su_struct();


% Reading first trace
tracl1=1;
tr(tracl1).tracl = fread(fid,1,'int');
tr(tracl1).tracr = fread(fid,1,'int');
tr(tracl1).fldr = fread(fid,1,'int');
tr(tracl1).tracf = fread(fid,1,'int');
tr(tracl1).ep = fread(fid,1,'int');
tr(tracl1).cdp = fread(fid,1,'int');
tr(tracl1).cdpt = fread(fid,1,'int');
tr(tracl1).trid = fread(fid,1,'short');
tr(tracl1).nvs = fread(fid,1,'short');
tr(tracl1).nhs = fread(fid,1,'short');
tr(tracl1).duse = fread(fid,1,'short');
tr(tracl1).offset = fread(fid,1,'int');
tr(tracl1).gelev = fread(fid,1,'int');
tr(tracl1).selev = fread(fid,1,'int');
tr(tracl1).sdepth = fread(fid,1,'int');
tr(tracl1).gdel = fread(fid,1,'int');
tr(tracl1).sdel = fread(fid,1,'int');
tr(tracl1).swdep = fread(fid,1,'int');
tr(tracl1).gwdep = fread(fid,1,'int');
tr(tracl1).scalel = fread(fid,1,'short');
tr(tracl1).scalco = fread(fid,1,'short');
tr(tracl1).sx = fread(fid,1,'int');
tr(tracl1).sy = fread(fid,1,'int');
tr(tracl1).gx = fread(fid,1,'int');
tr(tracl1).gy = fread(fid,1,'int');
tr(tracl1).counit = fread(fid,1,'short');
tr(tracl1).wevel = fread(fid,1,'short');
tr(tracl1).swevel = fread(fid,1,'short');
tr(tracl1).sut = fread(fid,1,'short');
tr(tracl1).gut = fread(fid,1,'short');
tr(tracl1).sstat = fread(fid,1,'short');
tr(tracl1).gstat = fread(fid,1,'short');
tr(tracl1).tstat = fread(fid,1,'short');
tr(tracl1).laga = fread(fid,1,'short');
tr(tracl1).lagb = fread(fid,1,'short');
tr(tracl1).delrt = fread(fid,1,'short');
tr(tracl1).muts = fread(fid,1,'short');
tr(tracl1).mute = fread(fid,1,'short');
tr(tracl1).ns = fread(fid,1,'ushort');
tr(tracl1).dt = fread(fid,1,'ushort');
tr(tracl1).gain = fread(fid,1,'short');	
tr(tracl1).igc = fread(fid,1,'short');
tr(tracl1).igi = fread(fid,1,'short');
tr(tracl1).corr = fread(fid,1,'short');
tr(tracl1).sfs = fread(fid,1,'short');
tr(tracl1).sfe = fread(fid,1,'short');
tr(tracl1).slen = fread(fid,1,'short');
tr(tracl1).styp = fread(fid,1,'short');
tr(tracl1).stas = fread(fid,1,'short');
tr(tracl1).stae = fread(fid,1,'short');
tr(tracl1).tatyp = fread(fid,1,'short');
tr(tracl1).afilf = fread(fid,1,'short');
tr(tracl1).afils = fread(fid,1,'short');
tr(tracl1).nofilf = fread(fid,1,'short');
tr(tracl1).nofils = fread(fid,1,'short');
tr(tracl1).lcf = fread(fid,1,'short');
tr(tracl1).hcf = fread(fid,1,'short');
tr(tracl1).lcs = fread(fid,1,'short');
tr(tracl1).hcs = fread(fid,1,'short');
tr(tracl1).year = fread(fid,1,'short');
tr(tracl1).day = fread(fid,1,'short');
tr(tracl1).hour = fread(fid,1,'short');
tr(tracl1).minute = fread(fid,1,'short');
tr(tracl1).sec = fread(fid,1,'short');
tr(tracl1).timbas = fread(fid,1,'short');
tr(tracl1).trwf = fread(fid,1,'short');
tr(tracl1).grnors = fread(fid,1,'short');
tr(tracl1).grnofr = fread(fid,1,'short');
tr(tracl1).grnlof = fread(fid,1,'short');
tr(tracl1).gaps = fread(fid,1,'short');
tr(tracl1).otrav = fread(fid,1,'short');
tr(tracl1).d1 = fread(fid,1,'float');
tr(tracl1).f1 = fread(fid,1,'float');
tr(tracl1).d2 = fread(fid,1,'float');
tr(tracl1).f2 = fread(fid,1,'float');
tr(tracl1).ungpow = fread(fid,1,'float');
tr(tracl1).unscale = fread(fid,1,'float');
tr(tracl1).ntr = fread(fid,1,'int');
tr(tracl1).mark = fread(fid,1,'short');
tr(tracl1).shortpad = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');
tr(tracl1).unass = fread(fid,1,'short');

%for ii=1:tr(tracl1).ns
%   tr(tracl1).trace(ii) = fread(fid,1,'float');
%end
tr(tracl1).trace = fread(fid,tr(tracl1).ns,'float');


% Calculation of the number of traces in the su-file
no_of_traces=fileparameters.bytes/(240+tr(1).ns*4);
if mod(no_of_traces,1)~=0     
   error('number of bytes differ from estimated number')
else
   no_of_traces=round(no_of_traces);
   disp(['File contains ' int2str(no_of_traces) ' traces.'])
   disp(' ')
end



% reading rest of traces
for tracl1=2:no_of_traces
	tr(tracl1).tracl = fread(fid,1,'int');
	tr(tracl1).tracr = fread(fid,1,'int');
	tr(tracl1).fldr = fread(fid,1,'int');
	tr(tracl1).tracf = fread(fid,1,'int');
	tr(tracl1).ep = fread(fid,1,'int');
	tr(tracl1).cdp = fread(fid,1,'int');
	tr(tracl1).cdpt = fread(fid,1,'int');
	tr(tracl1).trid = fread(fid,1,'short');
	tr(tracl1).nvs = fread(fid,1,'short');
	tr(tracl1).nhs = fread(fid,1,'short');
	tr(tracl1).duse = fread(fid,1,'short');
	tr(tracl1).offset = fread(fid,1,'int');
	tr(tracl1).gelev = fread(fid,1,'int');
	tr(tracl1).selev = fread(fid,1,'int');
	tr(tracl1).sdepth = fread(fid,1,'int');
	tr(tracl1).gdel = fread(fid,1,'int');
	tr(tracl1).sdel = fread(fid,1,'int');
	tr(tracl1).swdep = fread(fid,1,'int');
	tr(tracl1).gwdep = fread(fid,1,'int');
	tr(tracl1).scalel = fread(fid,1,'short');
	tr(tracl1).scalco = fread(fid,1,'short');
	tr(tracl1).sx = fread(fid,1,'int');
	tr(tracl1).sy = fread(fid,1,'int');
	tr(tracl1).gx = fread(fid,1,'int');
	tr(tracl1).gy = fread(fid,1,'int');
	tr(tracl1).counit = fread(fid,1,'short');
	tr(tracl1).wevel = fread(fid,1,'short');
	tr(tracl1).swevel = fread(fid,1,'short');
	tr(tracl1).sut = fread(fid,1,'short');
	tr(tracl1).gut = fread(fid,1,'short');
	tr(tracl1).sstat = fread(fid,1,'short');
	tr(tracl1).gstat = fread(fid,1,'short');
	tr(tracl1).tstat = fread(fid,1,'short');
	tr(tracl1).laga = fread(fid,1,'short');
	tr(tracl1).lagb = fread(fid,1,'short');
	tr(tracl1).delrt = fread(fid,1,'short');
	tr(tracl1).muts = fread(fid,1,'short');
	tr(tracl1).mute = fread(fid,1,'short');
	tr(tracl1).ns = fread(fid,1,'ushort');
	tr(tracl1).dt = fread(fid,1,'ushort');
	tr(tracl1).gain = fread(fid,1,'short');	
	tr(tracl1).igc = fread(fid,1,'short');
	tr(tracl1).igi = fread(fid,1,'short');
	tr(tracl1).corr = fread(fid,1,'short');
	tr(tracl1).sfs = fread(fid,1,'short');
	tr(tracl1).sfe = fread(fid,1,'short');
	tr(tracl1).slen = fread(fid,1,'short');
	tr(tracl1).styp = fread(fid,1,'short');
	tr(tracl1).stas = fread(fid,1,'short');
	tr(tracl1).stae = fread(fid,1,'short');
	tr(tracl1).tatyp = fread(fid,1,'short');
	tr(tracl1).afilf = fread(fid,1,'short');
	tr(tracl1).afils = fread(fid,1,'short');
	tr(tracl1).nofilf = fread(fid,1,'short');
	tr(tracl1).nofils = fread(fid,1,'short');
	tr(tracl1).lcf = fread(fid,1,'short');
	tr(tracl1).hcf = fread(fid,1,'short');
	tr(tracl1).lcs = fread(fid,1,'short');
	tr(tracl1).hcs = fread(fid,1,'short');
	tr(tracl1).year = fread(fid,1,'short');
	tr(tracl1).day = fread(fid,1,'short');
	tr(tracl1).hour = fread(fid,1,'short');
	tr(tracl1).minute = fread(fid,1,'short');
	tr(tracl1).sec = fread(fid,1,'short');
	tr(tracl1).timbas = fread(fid,1,'short');
	tr(tracl1).trwf = fread(fid,1,'short');
	tr(tracl1).grnors = fread(fid,1,'short');
	tr(tracl1).grnofr = fread(fid,1,'short');
	tr(tracl1).grnlof = fread(fid,1,'short');
	tr(tracl1).gaps = fread(fid,1,'short');
	tr(tracl1).otrav = fread(fid,1,'short');
	tr(tracl1).d1 = fread(fid,1,'float');
	tr(tracl1).f1 = fread(fid,1,'float');
	tr(tracl1).d2 = fread(fid,1,'float');
	tr(tracl1).f2 = fread(fid,1,'float');
	tr(tracl1).ungpow = fread(fid,1,'float');
	tr(tracl1).unscale = fread(fid,1,'float');
	tr(tracl1).ntr = fread(fid,1,'int');
	tr(tracl1).mark = fread(fid,1,'short');
	tr(tracl1).shortpad = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	tr(tracl1).unass = fread(fid,1,'short');
	
	%for ii=1:tr(tracl1).ns
 	%	tr(tracl1).trace(ii) = fread(fid,1,'float');
	%end
	tr(tracl1).trace = fread(fid,tr(tracl1).ns,'float');
	
	if tr(tracl1).ns~=tr(1).ns
		error(['number of samples differ in trace no: ' int2str(tracl1)])
	end
	
end

fclose(fid);


   
