function outseis(outfilename,data,source_x,source_y,switch_comp);

fid=fopen(outfilename,'w');

for tracl1=1:length(data)      % /*SEGY (without file-header)*/
   xr=source_x+data(tracl1).offset;
   yr=source_y;
   x=xr;
   y=yr;
   tr.tracl=tracl1;     
   tr.cdp=0;
   tr.trid=1;           %/* trace identification code: 1=seismic*/
   tr.offset=data(tracl1).offset*10^3;
   tr.gelev=yr*1000.0;
   tr.sdepth=source_y*1000.0;  % /* source depth (positive) */
   % /* angle between receiver position and reference point
   % (sperical coordinate system: swdep=theta, gwdep=phi) */
   tr.swdep=0;
   tr.scalel=-1000;
   tr.scalco=-1000;
   tr.sx=source_x*1000.0;         % /* X source coordinates */
   % /* group coordinates */
   tr.gx=(xr*1000.0);

   tr.ns=data(tracl1).npts;        % /* number of samples in this trace */
   tr.dt=data(tracl1).delta*1.0e6; % /* sample interval in micro-seconds */
   tr.d1=data(tracl1).delta;       % /* sample spacing for non-seismic data */

			
   tr.tracr=0;                     % 	/* trace sequence number within reel */

			tr.fldr=0; %	/* field record number */

			tr.tracf=0; %      ;	/* trace number within field record */

			tr.ep=0;%         ;	/* energy source point number */

			tr.cdpt=0;%       ;	/* trace number within CDP ensemble */

	
			tr.nvs=0;%       	;   /* number of vertically summed traces (see vscode
			   		%	in bhed structure) */

			tr.nhs=0; %       	;   /* number of horizontally summed traces (see vscode
			  		%	in bhed structure) */

			tr.duse=0;%     	;   /* data use:
					%	1 = production
						%2 = test */

			tr.selev=0;%     	; /* source elevation from sea level
			   		%	(above sea level is positive) */

			tr.gdel=0;%      	; /* datum elevation at receiver group */

			tr.sdel=0;%      	; /* datum elevation at source */

			tr.gwdep=0;%    	; /* water depth at receiver group */

			
			tr.sy=0;%      	;   /* Y source coordinate */

	
			tr.gy=0;%       	;   /* Y group coordinate */

			tr.counit=0;%    	;   /* coordinate units code:
					%	for previous four entries
						%1 = length (meters or feet)
						%2 = seconds of arc (in this case, the
						%X values are longitude and the Y values
						%are latitude, a positive value designates
						%the number of seconds east of Greenwich
						%or north of the equator */

			tr.wevel=0;%     ;	/* weathering velocity */

			tr.swevel=0;%    ;	/* subweathering velocity */

			tr.sut=0;%       ;	/* uphole time at source */

			tr.gut=0;%       ;	/* uphole time at receiver group */

			tr.sstat=0;%     ;	/* source static correction */

			tr.gstat=0;%     ;	/* group static correction */

			tr.tstat=0;%     ;	/* total static applied */

			tr.laga=0;%      ; /* lag time A, time in ms between end of 240-
			   		%	byte trace identification header and time
			   			%break, positive if time break occurs after
			   			%end of header, time break is defined as
			   			%the initiation pulse which maybe recorded
			   			%on an auxiliary trace or as otherwise
			   			%specified by the recording system */

			tr.lagb=0;%     	; /* lag time B, time in ms between the time break
			   		%	and the initiation time of the energy source,
			   			%may be positive or negative */

			tr.delrt=0;%     	; /* delay recording time, time in ms between
			   		%	initiation time of energy source and time
			   		%	when recording of data samples begins
			   		%	(for deep water work if recording does not
			   		%	start at zero time) */

			tr.muts=0;%      ; /* mute time--start */

			tr.mute=0;%      ; /* mute time--end */

	
			tr.gain=0;%      ; /* gain type of field instruments code:
					%	1 = fixed
					%	2 = binary
					%	3 = floating point
					%	4 ---- N = optional use */

			tr.igc=0;%       ; /* instrument gain constant */

			tr.igi=0;%       ; /* instrument early or initial gain */

			tr.corr=0;%      ; /* correlated:
					%	1 = no
					%	2 = yes */

			tr.sfs=0;%       ; /* sweep frequency at start */

			tr.sfe=0;%       ; /* sweep frequency at end */

			tr.slen=0;%      ; /* sweep length in ms */

			tr.styp=0;%      ; /* sweep type code:
				%		1 = linear
				%		2 = cos-squared
				%		3 = other */

			tr.stas=0;%;      	; /* sweep trace length at start in ms */

			tr.stae=0;%      	; /* sweep trace length at end in ms */

			tr.tatyp=0;%     	; /* taper type: 1=linear, 2=cos^2, 3=other */

			tr.afilf=0;%     	; /* alias filter frequency if used */

			tr.afils=0;%     	; /* alias filter slope */

			tr.nofilf=0;%    	; /* notch filter frequency if used */

			tr.nofils=0;%   	; /* notch filter slope */

			tr.lcf=0;%      	; /* low cut frequency if used */

			tr.hcf=0;%      	; /* high cut frequncy if used */

			tr.lcs=0;%       	; /* low cut slope */

			tr.hcs=0;%       	; /* high cut slope */

			tr.year=0;%      	; /* year data recorded */

			tr.day=0;%       	; /* day of year */

			tr.hour=0;%     	; /* hour of day (24 hour clock) */

			tr.minute=0;%    	; /* minute of hour */

			tr.sec=0;%       	; /* second of minute */

			tr.timbas=0;%    	; /* time basis code:
						% 1 = local
						% 2 = GMT
						% 3 = other */

			tr.trwf=0;%      	; /* trace weighting factor, defined as 1/2^N
			   			% volts for the least sigificant bit */

			tr.grnors=0;%   	; /* geophone group number of roll switch
			   			% position one */
% 
			tr.grnofr=0;%    	; /* geophone group number of trace one within
			   			% original field record */

			tr.grnlof=0;%    	; /* geophone group number of last trace within
			   			% original field record */

			tr.gaps=0;%      	;  /* gap size (total number of groups dropped) */

			tr.otrav=0;%     	;  /* overtravel taper code:
						% 1 = down (or behind)
						% 2 = up (or ahead) */

			%/* local assignments */

			tr.f1=0.0;	%/* first sample location for non-seismic data */

			tr.d2=0.0;	%/* sample spacing between traces */

			tr.f2=0.0;	%/* first trace location */

			tr.ungpow=0.0;	%/* negative of power used for dynamic
			   		%	range compression */

			tr.unscale=0.0;	%/* reciprocal of scaling factor to normalize
			   		%	range */
			tr.ntr=0; %  /* number of traces */

			tr.mark=0;
			
			tr.unass=0;
			
			if strcmp(switch_comp,'vr')
			   tr.data=data(tracl1).trace_vr;
			elseif strcmp(switch_comp,'ur')
			   tr.data=data(tracl1).trace_ur;
			else
			   tr.data=data(tracl1).trace_p;
			end


% Schreiben der Spur
count = fwrite(fid,tr.tracl,'int');
count = fwrite(fid,tr.tracr,'int');
count = fwrite(fid,tr.fldr,'int');
count = fwrite(fid,tr.tracf,'int');
count = fwrite(fid,tr.ep,'int');
count = fwrite(fid,tr.cdp,'int');
count = fwrite(fid,tr.cdpt,'int');
count = fwrite(fid,tr.trid,'short');
count = fwrite(fid,tr.nvs,'short');
count = fwrite(fid,tr.nhs,'short');
count = fwrite(fid,tr.duse,'short');
count = fwrite(fid,tr.offset,'int');
count = fwrite(fid,tr.gelev,'int');
count = fwrite(fid,tr.selev,'int');
count = fwrite(fid,tr.sdepth,'int');
count = fwrite(fid,tr.gdel,'int');
count = fwrite(fid,tr.sdel,'int');
count = fwrite(fid,tr.swdep,'int');
count = fwrite(fid,tr.gwdep,'int');
count = fwrite(fid,tr.scalel,'short');
count = fwrite(fid,tr.scalco,'short');
count = fwrite(fid,tr.sx,'int');
count = fwrite(fid,tr.sy,'int');
count = fwrite(fid,tr.gx,'int');
count = fwrite(fid,tr.gy,'int');
count = fwrite(fid,tr.counit,'short');
count = fwrite(fid,tr.wevel,'short');
count = fwrite(fid,tr.swevel,'short');
count = fwrite(fid,tr.sut,'short');
count = fwrite(fid,tr.gut,'short');
count = fwrite(fid,tr.sstat,'short');
count = fwrite(fid,tr.gstat,'short');
count = fwrite(fid,tr.tstat,'short');
count = fwrite(fid,tr.laga,'short');
count = fwrite(fid,tr.lagb,'short');
count = fwrite(fid,tr.delrt,'short');
count = fwrite(fid,tr.muts,'short');
count = fwrite(fid,tr.mute,'short');
count = fwrite(fid,tr.ns,'ushort');
count = fwrite(fid,tr.dt,'ushort');
count = fwrite(fid,tr.gain,'short');	
count = fwrite(fid,tr.igc,'short');
count = fwrite(fid,tr.igi,'short');
count = fwrite(fid,tr.corr,'short');
count = fwrite(fid,tr.sfs,'short');
count = fwrite(fid,tr.sfe,'short');
count = fwrite(fid,tr.slen,'short');
count = fwrite(fid,tr.styp,'short');
count = fwrite(fid,tr.stas,'short');
count = fwrite(fid,tr.stae,'short');
count = fwrite(fid,tr.tatyp,'short');
count = fwrite(fid,tr.afilf,'short');
count = fwrite(fid,tr.afils,'short');
count = fwrite(fid,tr.nofilf,'short');
count = fwrite(fid,tr.nofils,'short');
count = fwrite(fid,tr.lcf,'short');
count = fwrite(fid,tr.hcf,'short');
count = fwrite(fid,tr.lcs,'short');
count = fwrite(fid,tr.hcs,'short');
count = fwrite(fid,tr.year,'short');
count = fwrite(fid,tr.day,'short');
count = fwrite(fid,tr.hour,'short');
count = fwrite(fid,tr.minute,'short');
count = fwrite(fid,tr.sec,'short');
count = fwrite(fid,tr.timbas,'short');
count = fwrite(fid,tr.trwf,'short');
count = fwrite(fid,tr.grnors,'short');
count = fwrite(fid,tr.grnofr,'short');
count = fwrite(fid,tr.grnlof,'short');
count = fwrite(fid,tr.gaps,'short');
count = fwrite(fid,tr.otrav,'short');
count = fwrite(fid,tr.d1,'float');
count = fwrite(fid,tr.f1,'float');
count = fwrite(fid,tr.d2,'float');
count = fwrite(fid,tr.f2,'float');
count = fwrite(fid,tr.ungpow,'float');
count = fwrite(fid,tr.unscale,'float');
count = fwrite(fid,tr.ntr,'int');
count = fwrite(fid,tr.mark,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');
count = fwrite(fid,tr.unass,'short');

for ii=1:length(data(tracl1).trace_ur)
   count = fwrite(fid,tr.data(ii),'float');
end
end
fclose(fid);
