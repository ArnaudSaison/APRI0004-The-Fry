function Polar = XFOIL2Mat(inputDir, inputFiles, autosave)
% XFOIL2MAT Aggregate multiple XFOIL polars into a single structure
%   This function aggregates different airfoil polars obtained with XFOIL
%   into a single structure for easier handling in other scripts.
%   To facilitate comparisons, the range of AOA possible for each polar is
%   trimmed to the largest common range of AOA. Data are interpolated within
%   this range in order to have matching AOA for all polars.
%
%   The function can be run without arguments, in this case it will prompt the
%   user for the input.
%   Three inputs arguments can also be provided to automate the function.
%   In this case, the user should specifiy the directory containing the XFOIL
%   polars, the specific files to load and if they want the resulting structure
%   to be saved automatically.
%
% <a href="https://gitlab.uliege.be/thlamb/xfoil2mat">Documentation (README)</a>
% <a href="https://gitlab.uliege.be/thlamb/xfoil2mat/-/issues">Report an issue</a>
% -----
%
% Synopsis:
%   XFOIL2MAT prompts the user for the files to aggregate
%
%   Polar = XFOIL2MAT('inputDir') try to combine all files from inputDir into a
%   single polar strcture.
%
%   Polar = XFOIL2MAT('inputDir','inputFiles') combines the files matching
%   inputFiles from inputDir into a single polar structure.
%
%   Polar = XFOIL2MAT('inputDir','inputFiles', true) combines the files matching
%   inputFiles from inputDir into a single polar structure and automatically
%   save the Polar structure.
%
% Inputs:
%   inputDir   : Path of the directory with the XFOIL data
%   inputFiles : Files to select (ex: '*' (default), '*0012*', ...)
%   autosave   : Automatically save results (true / false(default))
%
% Output:
%   Polar : Structure collecting the results
%
% Example:
%   XFOIL2MAT
%   XFOIL2MAT('test_data')
%   XFOIL2MAT('test_data', '*0012*')
%   XFOIL2MAT('test_data', {'0012_re1e5', '0012_re1e6'})
%   XFOIL2MAT('test_data', '*0012*', true)

% -----
% XFOIL (GPL v2): http://web.mit.edu/drela/Public/web/xfoil/
% -----
% Author: Thomas Lambert <t.lambert@uliege.be>
% ULiege - Aeroelasticity and Experimental Aerodynamics
% MIT License
% https://gitlab.uliege.be/thlamb/xfoil2mat
% ------------------------------------------------------------------------------


%% Constants

DEFAULT_INPUTFILES = '*';  % If not specified, all files of the dir included by default
DEFAULT_AUTOSAVE = false;  % Autosave disabled by default
DEFAULT_AUTONAME = true;  % Name saved file automatically
PREVENT_OVERWRITING = true; % Prevents overwriting of the files during autosave


%% Check inputs and assign defaults

if nargin >= 1
    
    % Check if inputDir is a string and if the directory exist
    validateattributes(inputDir, {'char','string'},{'nonempty','vector'}, ...
        mfilename(), 'inputDir', 1)
    
    if ~exist(inputDir, 'dir')
        error(['XFOIL2MAT:DirectoryNotFound: '...
            'The directory specified as inputDir (%s) can not be found!\n'], ...
            inputDir);
    end
    
    % Check if inputFiles is a string (and set the default value if needed)
    if nargin >= 2
        validateattributes(inputFiles, {'char','cell'},{'nonempty','vector'}, ...
            mfilename(), 'inputFiles', 2)
        
        if iscell(inputFiles)
            if ~iscellstr(inputFiles)
                error(['If inpuFiles is given as a cell array, it must '...
                    'contain only character vectors.']);
            end
        end
        
    else
        inputFiles = DEFAULT_INPUTFILES;
    end
    
end

% Check if autosave is properly given and set its default value if not
if nargin == 3
    validateattributes(autosave, {'logical'},{'scalar'}, mfilename(), ...
        'autosave', 3)
else
    autosave = DEFAULT_AUTOSAVE;
end


%% Look for the files

if nargin == 0
    % Prompt user for the files to read
    [allFileNames, fullpath] = uigetfile('*.txt', 'Select all XFOIL files to aggregate', ...
        'MultiSelect', 'on','textFiles');
    
else
    
    % Get absolute path
    fullpath = fullfile(pwd, inputDir);
    
    % Convert inputFiles to string for simpler handling
    inputFiles = string(inputFiles);
    
    % Get all files
    allFiles = [];
    for i = 1:length(inputFiles)
        inputFiles(i) = appendextension(inputFiles(i),'.txt'); % Add extension
        dummy = dir(fullfile(fullpath,inputFiles(i)));
        if isempty(dummy)
            warning('XFOIL2Mat:FileNotFound: Could not find file %s.', fullfile(fullpath,inputFiles(i)));
        end
        allFiles = [allFiles; dummy];
    end
    
    if isempty(allFiles)
        error(['XFOIL2Mat:NoFilesFound: Impossible to find ANY file ' ...
            'matching the input arguments.']);
    end
    
    allFileNames = cell(1,length(allFiles));
    for i = 1:length(allFiles)
        allFileNames{i} = allFiles(i).name;
    end
    
end

% Convert filenames to string array
allFileNames = string(allFileNames);

%% Process XFoil results

% Number of files found
nbFiles=length(allFileNames);


% Initialize structures
Polar = struct('airfoil', [], ...
    'reynolds',zeros(1,nbFiles), ...
    'mach', zeros(1,nbFiles), ...
    'nCrit', [], ...
    'alpha', [], ...
    'cl', [], ...
    'cd', [], ...
    'cm', []);
Tmp(nbFiles).resultArray=[];

% Initialize AOA range
aoa_first = -180;
aoa_last = 180;
aoa_step = 360;

for i=1:nbFiles
    
    % Load file and extract parameters and results
    fileID = fopen(fullfile(fullpath,allFileNames{i}));
    tmpAirfoil = textscan(fileID,'%s %s', 1, 'Delimiter','\n:','HeaderLines',3);
    tmpParam = textscan(fileID, '%s %s %f %s %s %f %s %f %s %s %f %*[^\n]', 1, 'HeaderLines', 4); % ! HeaderLines is a cursor. This will read line 9
    Tmp(i).resultArray = cell2mat(textscan(fileID,'%f %f %f %f %f %*[^\n]', 'HeaderLines',4)); % ! HeaderLines is a cursor. This will start reading at line 13
    fclose(fileID); clear fileID;
    
    % Extract the airfoil name
    airfoil_name = char(tmpAirfoil{2});
    Polar.airfoil{1,i} = sanitizestring(airfoil_name);
    
    % Get Mach, Reynolds and Ncrit numbers corresponding to the processed file
    Polar.mach(1,i) = tmpParam{3};
    Polar.reynolds(1,i) = tmpParam{6}*10^tmpParam{8};
    Polar.nCrit(1,i) = tmpParam{11};
    
    % Get the larger common AOA range and the smallest step
    if Tmp(i).resultArray(1,1) > aoa_first
        aoa_first = Tmp(i).resultArray(1,1);
    end
    if Tmp(i).resultArray(end,1) < aoa_last
        aoa_last = Tmp(i).resultArray(end,1);
    end
    if min(diff(Tmp(i).resultArray(:,1))) < aoa_step
        aoa_step = min(diff(Tmp(i).resultArray(:,1)));
    end
end

% Common range for the AOA
aoaRange = aoa_first:aoa_step:aoa_last;
Polar.alpha = deg2rad(aoaRange');

% Interpolate results over the same angle of attacks
for i=1:nbFiles
    Polar.cl(:,i) = interp1(Tmp(i).resultArray(:,1), Tmp(i).resultArray(:,2), aoaRange);
    Polar.cd(:,i) = interp1(Tmp(i).resultArray(:,1), Tmp(i).resultArray(:,3), aoaRange);
    Polar.cm(:,i) = interp1(Tmp(i).resultArray(:,1), Tmp(i).resultArray(:,5), aoaRange);
end


%% Output

% Sort everything by increasing Reynolds
[Polar.reynolds, indSort] = sort(Polar.reynolds);
Polar.airfoil(:) = Polar.airfoil(indSort);
Polar.mach(:)  = Polar.mach(indSort);
Polar.nCrit(:) = Polar.nCrit(indSort);
Polar.cl(:,:)  = Polar.cl(:,indSort);
Polar.cd(:,:)  = Polar.cd(:,indSort);
Polar.cm(:,:)  = Polar.cm(:,indSort);


%% Save results automatically
if nargin == 0 % Prompt only if user input is to be expected
    prompt = 'Save the resulting polar? Y/N [Y]: ';
    str = input(prompt,'s');
else
    str = 'N';
end

if strcmpi(str,'Y') || strcmpi(str,'yes') || isempty(str) || autosave
    
    % Prompt user for filename
    if ~DEFAULT_AUTONAME
        prompt = 'Filename for saved polar?: ';
        filename = input(prompt,'s');
    end
    
    % Automatic filename if user did not input one
    if DEFAULT_AUTONAME || isempty(filename)
        % Different names if only one or multiple files loaded
        if nbFiles == 1
            filename_base = [Polar.airfoil{1},'-single'];
        else
            filename_base = [Polar.airfoil{1},'-polars'];
        end
        
        % Prevent overwriting files
        filename = filename_base;
        if ~PREVENT_OVERWRITING
            i=0;
            while isfile(appendextension(filename,'.mat'))
                i=i+1;
                filename = [filename_base,'_',num2str(i)];
            end
        end
        
    end
    filename = appendextension(filename,'.mat');
    save(filename,'Polar')
end


end


%% HELPERS

function str = appendextension(str,ext)
% APPENDEXTENSION Append an extension to a string

if isempty(regexp(str,[ext,'$'],'once'))
    str = strcat(str,ext);
end

end


function str = sanitizestring(str)
%SANITIZESTRING Santitizes a string by replacing unwanted characters by better ones

%% Mapping table
% Maps unwanted chars with their replacements
mappingTable = {' ', '_';
    '%', 'pc'};


%% String sanitization

% Remove whitespaces at the beginning and end of the string
str = strtrim(str);

% Remove multiple consecutive whitespaces
str = regexprep(str, ' +', ' ');

% Replace characters according to mapping table
for i = 1:size(mappingTable,1)
    str = replace(str, mappingTable{i,1},mappingTable{i,2});
end

end
