%% PETROUEDGE - Edge detection based on optimal ramp filter.
%
%% Description
% Implements the edge detector of [PK91] derived from optimal filters for 
% ramp edges.
%
%% Syntax
%     emap = PETROUEDGE(I);
%     emap = PETROUEDGE(I, mask);
%
%% Inputs
% *|I|* : input image, possibly multichannel.
%      
% *|mask|* : ((optional) indice of the mask used by Petrou & Kittler edge
%     detector for smoothing the input image; it is either 0, 1 or 2; 
%     default: |mask=0|.
%
%% Outputs
% *|emap|* : edge map.
%
%% References
% [Spac86]  L.A. Spacek: "Edge detection and motion detection", Image &
%      Vision Computing, 4:43-56, 1986.
%      <http://www.sciencedirect.com/science/article/pii/0262885686900077>
%      
% [PK91]  M. Petrou and J. Kittler: "Optimal edge detectors for ramp
%      edges", IEEE Trans on Pattern Analysis and Machine Intelligence, 
%      13(5):483-491, 1991.
%      <http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=134047>
%
%% See also
% Related:
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <CANNYEDGEMAP.html |CANNYEDGEMAP|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <CONGRUENCYEDGE.html |CONGRUENCYEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>.
% Called: 
% <PETROUEDGE_BASE.html |PETROUEDGE_BASE|>.

%% Function implementation
function edgemap = petrouedge(I,varargin)

%%
% parsing and checking parameters

error(nargchk(1, 10, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(I)
    error('petrouedge:inputparameter','matrix required in input'); 
end

p = createParser('PETROUEDGE');   % create an instance of the inputParser class.
p.addOptional('mask',0, @(x) isscalar(x) && any(x==0 |x==1|x==2)); 

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);                                                            


%% 
% main calculation
edgemap = petrouedge_base(I, p.mask);

%%
% display

if p.disp
    figure, imagesc(edgemap), axis image off, title('Petrou&Kittler edge map');
    if size(edgemap,3) == 1, colormap gray;  end;
end

end % end of petrouedge
