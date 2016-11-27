%% READENVIRASTER - Read an image of ENVI standard type to a MATLAB array.
%
%% Syntax
%       image = READENVIRASTER(fname);
%       [image,p,t] = READENVIRASTER(fname, ext);
%
%% Input
% *|fname|* : string giving the full pathname of the ENVI image to read.
%
%% Outputs
% *|image|* : array of size |(c,l,b)| containing the ENVI image values organised 
%	  in |c| cols, |l| lines and |b| bands.
%
% *|p|* : vector of size |(1,3)| that contains (1) the nb of cols, (2) the 
%     number of lines and (3) the number of bands of the opened image.
%
% *|t|* : string describing the image data type string in MATLAB conventions.
%
%% Remark
% |READENVIRASTER| needs the corresponding image header file generated 
% automatically by ENVI. The ENVI header file must have the same name as the 
% ENVI image file + the '.hdf' exention.
%
%% See also 
% Related:
% <READENVIROI.html |READENVIROI|>.

%% Function implementation
function [image,p,t] = readenviraster(fname,varargin)

%%
% parameters initialization
elements={'samples ' 'lines   ' 'bands   ' 'bands ' 'data type '};
d={'bit8' 'int16' 'int32' 'float32' 'float64' 'uint16' 'uint32' 'int64' 'uint64'};

%%
% check user input
if ~ischar(fname)
    error('readenviraster:errorinput', 'fname should be a char string');
end

if nargin>=2
    ext=varargin{1};
    if ~isempty(ext) && ~strcmp(ext(1),'.')
        ext=strcat('.',ext);
    end
    if nargin==3,  disp = varargin{2};  
    else           disp = false;        end
else % default: empty string
    ext = '';
    disp = false;
end

%%
% open ENVI header file to retreive s, l, b & d variables
rfid = fopen(strcat(fname,'.hdr'),'r');

%%
% check if the header file is correctely open
if rfid == -1
    error('readenviraster:readerror','input header file does not exist');
end;

%%
% read ENVI image header file and get p(1) : nb samples,
% p(2) : nb lines, p(3) : nb bands and t : data type
while 1
    tline = fgetl(rfid);
    if ~ischar(tline), break, end
    [first,second]=strtok(tline,'=');
    
    switch first
        case elements(1)
            [~,s] = strtok(second);
            p(1) = str2double(s);
        case elements(2)
            [~,s] = strtok(second);
            p(2) = str2double(s);
        case elements(3)
            [~,s] = strtok(second);
            p(3) = str2double(s);
        case elements(4)
            [~,s] = strtok(second);
            p(3) = str2double(s);
        case elements(5)
            [~,s] = strtok(second);
            t = str2double(s);
            switch t
                case 1
                    t = d(1);
                case 2
                    t = d(2);
                case 3
                    t = d(3);
                case 4
                    t = d(4);
                case 5
                    t = d(5);
                case 12
                    t = d(6);
                case 13
                    t = d(7);
                case 14
                    t = d(8);
                case 15
                    t = d(9);
                otherwise
                    error('readenviraster:errorinput', 'unknown image data type');
            end
            
    end
end

fclose(rfid);

t = t{1,1};

%%
% open the ENVI image and store it in the 'image' MATLAB array
if disp
    disp([('Opening '), ...
        (num2str(p(1))),('cols x '),...
        (num2str(p(2))),('lines x '),...
        (num2str(p(3))),('bands'),...
        ('of type '), (t), (' image...')]);
end

%%
% read and store into an image
fid = fopen(strcat(fname,ext));
image = fread(fid,t);
image = reshape(image,[p(1),p(2),p(3)]);

%%
% close the file
fclose(fid);
% disp([('Image data type : '),(t)])
end % end of readenviraster