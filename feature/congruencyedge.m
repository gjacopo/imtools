%% CONGRUENCYEDGE - Wrapper for edge detection based on phase congruency.
%
%% Description
% Nothing else than a wrapper function for P.Kovesi's edge and corner detector
% based on phase congrunecy [Koves99,Koves02,Koves03]. 
% See function |PHASECONG3| in distribution [Kovesi] if you'd rather use the
% original function.
%
%% Syntax
%     [M, m]  = CONGRUENCYEDGE(I);
%     [M m or ft pc EO T]  = ...
%          CONGRUENCYEDGE(I, nscale, norient, minwave, mult, ...
%          'Property', propertyvalue, ...);
%
%% Inputs
% See |PHASECONG3| for input parameters.
%
% *|I|* : input image, possibly multichannel.
%      
% *|nscale|* : number of wavelet scales; default: |nscale=4|. 
%      
% *|norient|* : number of filter orientations; default: |norient=6|.
%      
% *|minwave|* : wavelength of smallest scale filter; default: |minwave=3|.
%      
% *|mult|* : scaling factor between successive filters; default: |mult=2.1|.
%
%% Property [propertyname  propertyvalues]
% *|'sigonof'|* : ratio of the standard deviation of the Gaussian describing
%     the log Gabor filter's transfer function in the frequency domain to
%     the filter center frequency; default: |sigonof=0.55|.
%      
% *|'k'|* : standard deviations of the noise energy beyond the mean at which
%     the noise threshold point is set; default: |k=2|.
%
% *|'cutoff'|* : fractional measure of frequency spread below which phase 
%     congruency values get penalized; default: |cutoff=0.5|.
%      
% *|'g'|* : sharpness of the transition in the sigmoid function used to 
%     weight phase congruency for frequency spread; default: |g=10|.
%      
% *|'noise'|* : parameter specifying the method used to determine noise 
%     statistics; it is either: 
%      
% * -1 to use median of smallest scale filter responses,
% * -2 to use mode of smallest scale filter responses,
% * or any >=0 value to set a fixed noise threshold;
%      
% default: |noise=-1|.
%
%% Outputs
% *|M|* : maximum moment of phase congruency covariance; used as an 
%     indicator of edge strength.
%      
% *|m|* : minimum moment of phase congruency covariance; used as an 
%     indicator of corner strength.
%      
% *|or|* : orientation image in integer degrees 0-180, positive 
%     anticlockwise; 0 corresponds to a vertical edge, 90 is horizontal.
%      
% *|ft|* : local weighted mean phase angle at every point in the image; a 
%     value of pi/2 corresponds to a bright line, 0 corresponds to a step 
%     and -pi/2 is a dark line.
%      
% *|pc|* : cell array of phase congruency images (values between 0 and 1) for  
%     each orientation
%      
% *|EO|* : a 2D cell array of complex valued convolution result
%      
% *|T|* : calculated noise threshold
%
%% Acknowledgment
% This function is nothing else than a wrapper for P.Kovesi functions
% implementing congruency analysis.
% 
%% References
% [Koves99]  P.D. Kovesi: "Image features from phase congruency", Videre:
%      Journal of Computer Vision Research, 1:1-26, 1999.
%      <http://mitpress.mit.edu/e-journals/Videre/001/articles/v1n3001.pdf>
%      
% [Koves02]  P.D. Kovesi: "Edges are not just steps", Proc. ACCV, pp. 822-
%      827, 2002.
%      <www.csse.uwa.edu.au/~pk/research/pkpapers/ACCV62.pdf>
%      
% [Koves03]  P.D. Kovesi: "Phase congruency detects corners and edges",
%      Proc. DICTA, pp. 309-318, 2003. 
%      <www.csse.uwa.edu.au/~pk/research/pkpapers/phasecorners.pdf>
%
% [Kovesi]   P.D. Kovesi: "MATLAB and Octave Functions for Computer Vision
%      and Image Processing", The University of Western Australia, available
%      at <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
%
%% See also
% <../../../toolboxes/kovesi/PhaseCongruency/html/PHASECONG3.html |PHASECONG3|>,
% <EDGECORNER.html |EDGECORNER|>,
% <CANNYEDGE.html |CANNYEDGE|>,
% <CANNYEDGEPROD.html |CANNYEDGEPROD|>,
% <SDGDEDGE.html |SDGDEDGE|>,
% <KOETHEDGE.html |KOETHEDGE|>,
% <COMPASSEDGE.html |COMPASSEDGE|>,
% <ANISOEDGE.html |ANISOEDGE|>,
% <ELDERZUCKEREDGE.html |ELDERZUCKEREDGE|>,
% <PETROUEDGE.html |PETROUEDGE|>,
% <ROTHWELLEDGE.html |ROTHWELLEDGE|>.
% Called: 
% <CONGRUENCYEDGE_BASE.html |CONGRUENCYEDGE_BASE|>.

%% Function implementation
function [M, m, varargout] = congruencyedge(I, varargin)

%%
% check if possible

if ~exist('phasecong3','file')
    error('congruencyedge:libraryerror','Kovesi''s library required');
end

%%
% parsing and checking parameters

error(nargchk(1, 23, nargin, 'struct'));
error(nargoutchk(1, 7, nargout, 'struct'));

if ~isnumeric(I)
    error('congruencyedge:inputerror','matrix required in input'); 
end

p = createParser('CONGRUENCYEDGE');   
% optional parameters
% See PHASECON3 for default parameters' values
p.addOptional('nscale', 4, @(x)isscalar(x) && round(x)==x && x>0);
p.addOptional('norient', 6, @(x)isscalar(x) && round(x)==x && x>0);
p.addOptional('minwave', 3, @(x)isscalar(x) && round(x)==x && x>0);
p.addOptional('mult', 2.1, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('sigonf', 0.55, @(x)isscalar(x) && isfloat(x) && x>=0 && x<=1);
p.addParamValue('k', 2, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('cutoff', 0.5, @(x)isscalar(x) && isfloat(x) && x>=0 && x<=1);
p.addParamValue('g', 10, @(x)isscalar(x) && isfloat(x) && x>0);
p.addParamValue('noise', -1, @(x)isscalar(x) && ...
    isfloat(x) && (x>0 || x==-1 || x==-2));

% parse and validate all input arguments
p.parse(varargin{:});
p = getvarParser(p);                                                            

%%
% main processing

[M, m, or, ft, pc, EO, T]  = ...
    congruencyedge_base(I, p.nscale, p.norient, p.minwave, ...
    p.mult, p.sigonf, p.k, p.cutoff, p.g, p.noise);

%%
% set outputs

if nargout>=3,
    varargout{1} = or;
    if nargout>=4,
        varargout{2} = ft;
        if nargout>=5,
            varargout{3} = pc;
            if nargout>=6,
                varargout{4} = EO;
                if nargout==7,
                    varargout{5} = T;
                end
            end
        end
    end
end

%%
% display 

if p.disp
    C = size(I,3);
    figure,
    subplot(1,nargout,1), imagesc(rescale(M,0,1)), axis image off
    title('phase congruency covariance - max');
    if nargout==2
        subplot(1,2,2), imagesc(rescale(m,0,1)), axis image off
        title('phase congruency covariance - min');
    end
    if C==1,  colormap gray;   end;
end

end % end of congruencyedge
