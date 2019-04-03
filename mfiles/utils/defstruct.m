function defstruct(name,fields,values)
% DEFSTRUCT(name,fields,values)
%
% Assigns default values to a named structure
%
% INPUT:
% 
% name    A string, enclosed in single quotes, with a structure name
% fields  A cell array with fields that you want the structure to have 
% values  A cell array with the values by which you populate those fields
%
% OUTPUT:
%
%      None. The variable appears as if by magic into your workspace or
%      will be available inside your function.
%
% SEE ALSO: 
%
% CELLNAN, STRUCTNAN, STRUCT
%
% Last modified by fjsimons-at-alum.mit.edu, 07/02/2014

if ~ischar(name),
  error(sprintf(['The first argument of DEFSTRUCT',...
		'has to be a string with a structure name']));
end

% Always do it is our default here
si=1;

% Maybe it's only a single field name and value pair
if isstr(fields); 
  fields={fields}; 
end  
if ~iscell(values); 
  values={values}; 
end

% Commented out by Vladimir Kazei Check input lengths
%difer(length(values)-length(fields),[],1,NaN)

% If the STRUCTURE exists at all...
if evalin('caller',[ 'exist(''' name ''',''var'')']);
  % ... and it's empty, do it; but don't do it if it's non empty
  % This as the the whole structure variable, not any piece of it
  si=evalin('caller',[ 'isempty(' name ')']);
  % Now if the STRUCTURE is not empty, SOME of it may be empty
  if ~si
    % Which of the fieldnames exist?
    fn=evalin('caller',['fieldnames(' name ')'])';
    % For now they all exist
    isv=logical(zeros(1,length(fn)));
    % But now we know which ones are the empty ones
    for index=1:length(fn)
      isv(index)=evalin('caller',['isempty(' name '.' fn{index} ')']);
    end

    % Then there are the fields that did exist and were nonempty
    oldvalues=cellnan(1,1,sum(~isv)); ondex=0;
    for index=find(~isv)
      ondex=ondex+1;
      oldvalues{ondex}=evalin('caller',[ name '.' fn{index} ]);
    end

    % Variables that did exist but were empty for which we have values 
    [ifn,ivs]=intersect(fields,fn(isv),'stable');

    % Variables that existed but were empty for which we don't have news 
    bfn=setdiff(fn(isv),ifn);

    % Variables that didn't exist for which we have values
    [newfields,ivf]=setdiff(fields,fn,'stable');

    % Combine the new fields, the old values, and the new values and
    % remember that where we have no information the default is [] anyhow
    % Remember the difference between smooth and curly parentheses!
    fields={newfields{:} fn{~isv}      ifn{:}      bfn{:}};
    values={values{ivf}  oldvalues{:}  values{ivs}       };
    
    % So in the end you always have to do something!
    si=1;
  end
end

% Do it or not, either the whole thing or just the pieces of it
if si
  % Do it exactly like OPTIMSET
  structinput=cell(2,length(fields));
  % Fields go in first row
  structinput(1,:)=fields;
  % []'s go in second row, remember there may be hollow empties at the end
  structinput(2,1:length(values))=values;
  % Turn it into correctly ordered comma separated list and call struct
  S=struct(structinput{:});
  % Order output. NOT. Usually, we expect a certain ordering!
  % S=orderfields(S);
  % Put in the caller workspace  
  assignin('caller',name,S);
end
