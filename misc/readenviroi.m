%% READENVIROI - Read ENVI files.
%
%% Description
% This function reads in an ENVI ROI file and outputs a Struct containing
%
%   ROI Name:         ROI.'ROI Number'.'Name'     (Char)
%   ROI Color:        ROI.'ROI Number'.'Color'    (vector)
%   ROI Num Points:   ROI.'ROI Number'.numPoints  (scalar)
%   ROI Data:         ROI.'ROI Number'.Points     (matrix)
%
%% Syntax 
%     ROI = READENVIROI(roiFile);
%
%% Output
% *|ROI|* : the data in each cell of the struct is as outlined above. the matrix
%      is of size # of bands plus the ID number, and (X,Y) coordinates by the
%      number of points in the |ROI|.
%
%% Remark
% Note that this function assumes a strict formatting standard for the ROI
% file.  IT IS DESIGNED TO ONLY WORK ON THE ENVI ROI FILES SAVED IN ASCII 
% FORMAT.  The ROI points must follow immediately after the last ROI
% information line in the file top matter.  There is a commented out while
% loop that was intended to test for the start of the ROI points, but was
% removed for convience.  However, a more robust algorithm will test for
% the start of the ROI points so as to avoid strict formatting errors.
% Future improvements to this code may incorporate reading additional ROI 
% file formats and/or the test for teh start of the roi points.
%
%% Acknowledgment
% Written by Jared Herweg, Rochester Institute of Technology, May 2010
% jxh6389@rit.edu
%
%% See also 
% Related:
% <READENVIRASTER.html |READENVIRASTER|>.

% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met: 
%   - Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer. 
%   - Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in the 
%     documentation and/or other materials provided with the distribution.
%   - Neither the name of Rochester Institute of Technology nor the names 
%     of its contributors may be used to endorse or promote products 
%     derived from this software without specific prior written permission.
% 
% Disclamer 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Function implementation
function ROI = readenviroi(roiFile)

%%
% open File
fid = fopen(roiFile);

%%
% get ROI header information
ROI = ''; % Initialize Struct
test = [];
while isempty(test)
    fline = fgetl(fid);
    % find number of ROIs
     test = find(strncmp('; Number of', fline, 11));
    % test = strmatch('; Number of',fline); 
end
ROI.('NumROI') = str2double(strtrim(fline(strfind(fline,':')+1:end)));
test = [];

%%
% for each ROI, get associated info
for count = 1:ROI.NumROI
    while isempty(test)
        fline = fgetl(fid);
        test = find(strncmp('; ROI name', fline, 10));
       % test = strmatch('; ROI name',fline);
    end
    ROIName = strtrim(fline(strfind(fline,':')+1:end));
    fline = fgetl(fid);
    ROIcolor = str2double(strtrim(fline(strfind(fline,'{')+1:end-1)));
    fline = fgetl(fid);
    ROI_numPts = strtrim(fline(strfind(fline,':')+1:end));
    % Dynamically expand struct and assign ROI information to cell
    % fields
    ROIcnt = ['R' num2str(count)]; % Each independent ROI name assigned
    ROI.(ROIcnt).('Name') = ROIName;
    ROI.(ROIcnt).('Color') = ROIcolor;
    ROI.(ROIcnt).('numPts') = str2double(ROI_numPts);
    test = [];
end

%%
% get the ROI Points for each ROI

% Advance to the start of the ROI points
%     while isempty(test)
fline = fgetl(fid);
%         test = strmatch(';   ID',fline);
%     end

%%
% get number of bands in ROI
NumBands = length(find(fline == 'B'));

%%
% add 3 fields to num bands to account for ROI point ID and the
% x,y coordinates (see ROI ASCII file).
NumBands = NumBands + 3;

%%
% for each ROI, get the associative ROI points
for count = 1:ROI.NumROI
    ROIcnt = ['R' num2str(count)];
    ROI.(ROIcnt).('Points') = zeros(ROI.(ROIcnt).numPts,NumBands);
    for count2 = 1:ROI.(ROIcnt).numPts
        fline = fgetl(fid);
        ROI.(ROIcnt).Points(count2,:) = str2double(fline);
    end
    % fline = fgetl(fid); % skips the blank line
end

%%
% close ROI File
fclose(fid);

end % end of readenviroi