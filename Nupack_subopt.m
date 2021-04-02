function [mfe, structure, bplist] = Nupack_subopt(sequence,energyGap,varargin)
%Last updated: 2015-10-20
%Usage: [structure, bplist, mfe] = nupackfold(sequence, ...)
%
%IMPORTANT: please make sure the main NUPACK folder (should contain the folders bin, doc, lib, parameters, src, etc.) is in the same directory as this script!
%IMPORTANT: please make sure that you are calling this script in its directory, not while in another folder!!!
%
%Input: sequence - either a string or a cell array of strings, for multi-strand
%       folding and/or hybridization.
%
%Optional Parameters:
%   "home" - homepath
%   "temp" - temperature of folding in celsius, defaults to 37C
%   "material" - DNA or RNA? defaults to DNA if not specified.
%   "sodium" - sodium concentration of solution, in Mol/L.
%   "Keepfiles" - Should nupack's output file be kept (Yes/no)? Defaults to
%       no. If yes is selected, file will be outputted to "out.mfe". WARNING:
%       If "out.mfe" already exists, it will overwritten if keepfiles = yes!!!!
%   "verbose" - Whether or not to display elapsed time for folding. By default is 'yes'. To turn off, use 'no'


%EXAMPLE:

%Finds MFE of the sequence using default settings (37 Degrees, assumes DNA, assumes 1.0M sodium)
%  mfe = nupackMFE('AGAGTTGGAGGGAAGAGG');

%Finds MFE, and structure of the sequence where the material is DNA, temperature is 25 degrees, and sodium is 0.15M
%  [mfe, structure] = nupackMFE('AGAGTTGGAGGGAAGAGG','material','dna','temperature',25,'sodium',0.15);

    TEMPERATURE = 37;
    MATERIAL = 'dna';
    IS_KEEP_FILES = 0;
    IS_VERBOSE = 0;
    SODIUM = 1;
    structure = '';
    bplist = [];
  	homepath = '/Users/jiaming/Dropbox\ \(Nablab\)/Low_Yield_Bisulfite_Conversion/nupack';
  	homepath2 = '/Users/jiaming/Dropbox (Nablab)/Low_Yield_Bisulfite_Conversion/nupack';
    
    argTable = vararginToTable(varargin,nargin-2);
    for i=1:size(argTable,1)
        switch lower(argTable{i,1});
			case {'homepath'}
				homepath = argTable{i,2};
            case {'temp', 'temperature'}
                if isnumeric(argTable{i,2})
                    TEMPERATURE = argTable{i,2};
                else
                    disp('Invalid temperature input!');
                    return
                end 
            case {'material','mat','type'}
                if sum(strcmpi({'dna','rna'},argTable{i,2}))
                    MATERIAL = lower(argTable{i,2});
                else
                    disp('Invalid material type!');
                    return
                end
            case {'keepfiles','keepfile'}
                if sum(strcmpi({'yes','y','true','1'}, argTable{i,2}))
                    IS_KEEP_FILES = 1;
                elseif sum(strcmpi({'no','n','false','0'}, argTable{i,2}))
                    IS_KEEP_FILES = 0;
                else
                    disp('Invalid choice for keepfiles option!');
                    return
                end
            case {'verbose'}
                if sum(strcmpi({'yes','y','true','1'}, argTable{i,2}))
                    IS_VERBOSE = 1;
                else
                    IS_VERBOSE = 0;
                end
            case {'salt','sodium'}
                if isnumeric(argTable{i,2})
                    SODIUM = argTable{i,2};
                else
                    disp('Invalid sodium concentration input!');
                    return
                end
            otherwise
                disp(['Invalid option ' argTable{i,1} '!']);
        end
	end
	
	setenv('NUPACKHOME', [homepath2 '/nupack']);
    
    sodiumString = '';
    %Check to make sure the option "sodium" is only used on DNA
    if SODIUM~=1 && ~strcmpi(MATERIAL,'dna')
        error('Cannot use "sodium" option on non-DNA material! NUPACK supports different sodium concentrations for DNA Only.');
    end
    if SODIUM~=1
        sodiumString = ['-sodium ' num2str(SODIUM)];
    end
    
    %Check if files need to be kept. Otherwise, create a temporary file.
    if IS_KEEP_FILES
        filename = 'out';
    else
        filename = 'nupack_tempfile';
    end
    
    if ~iscell(sequence)
        sequence = {sequence};
    end
    
    fid = fopen([filename '.in'],'w');
    if strcmpi(MATERIAL,'rna');
        MATERIAL = 'rna';
        for i=1:length(sequence)
            sequence{i} = regexprep(sequence{i},'T','U');
        end
    else
        MATERIAL = 'dna';
    end
    
    
    if length(sequence) > 1
        multisettings = '-multi ';
        fprintf(fid,'%d\n',length(sequence));
        for i=1:length(sequence)
            fprintf(fid,'%s\n',upper(sequence{i}));
        end
        fprintf(fid,[repmat('%d ',1,length(sequence)) '\n'],1:length(sequence));
        
    else
        multisettings = '';
        fprintf(fid,'%s\n',upper(sequence{1}));
    end
    fprintf(fid,'%.f',energyGap);
    fclose(fid);
    
    
    systemcallstr = [homepath, '/subopt ' multisettings '-T ' num2str(TEMPERATURE) ' -material ' MATERIAL ' ' sodiumString ' ' filename];
	
    if IS_VERBOSE
        disp(['System call: ' systemcallstr]);
    end
    tic;
    
    %Call NUPACK
    system(systemcallstr);
    finishtime = toc;
    if IS_VERBOSE
        disp(['Folding finished. Time elapsed: ' num2str(finishtime) ' seconds']);
    end
    
    %Parse NUPACK results from temporary output to MATLAB
    fid = fopen([filename '.subopt']);
    substructCount = 0;
    mfe = zeros(100,1);
    structure = cell(100,1);
    bplist = cell(100,1);
    while 1
        line = fgets(fid);
        if ~ischar(line); break; end;
        trimmedLine = strtrim(line);
        if isSeparator(trimmedLine) %Found beginning of substructure block
            clear substructure;
            substructCount = substructCount + 1;
            bplist{substructCount} = zeros(100,2);
            lineCount = 0;
            while 1
                lineCount = lineCount + 1;
                line = fgets(fid);
                if ~ischar(line); error('Error: encountered unexpected EOF while parsing NUPACK output'); end;
                if isSeparator(line) %End of substructure block
                    bplist{substructCount} = bplist{substructCount}(1:lineCount-4,:); 
                    break;
                end
                if lineCount==2
                    mfe(substructCount) = str2double(strtrim(line));
                elseif lineCount==3
                    structure{substructCount} = strtrim(line);
                elseif lineCount>3
                    bplist{substructCount}(lineCount-3,:) = str2num(strtrim(line)); %#ok<ST2NM>
                end
            end
        end
    end
    mfe = mfe(1:substructCount);
    structure = structure(1:substructCount);
    bplist = bplist(1:substructCount);
    fclose(fid);
    %Delete temporary files
    if ~IS_KEEP_FILES
        delete([filename '.in']);
        delete([filename '.subopt']);
    end
end

function boolean = isSeparator(str)
boolean = (length(str) > 2 && sum(str(2:end)=='%') > 5);
end

function argTable = vararginToTable(vararg,narg)
%Parses varargin options in the format "option1",option1value,.... into a
%table of ["option1" option1value;"option2" option2value;....]
    if mod(narg,2)==1
        disp('Error: Incorrect optional argument format!');
        disp(vararg);
        return
    elseif narg == 0
        argTable = {};
        return
    else
        argTable = reshape(vararg,2,narg/2)';
    end
end