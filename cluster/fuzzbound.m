%% FUZZBOUND - Fuzzy boundaries extraction.
%
%% Description
% Compute the fuzzy boundaries on a fuzzy categorical map following the 
% approach in [ZK00].
% 
%% Syntax
%     fmap = FUZZBOUND(fmemb);
%     fmap = FUZZBOUND(fmemb, 'propertyname',propertyvalue,...);
%
%% Input
% *|fmemb|* : array |(X,Y,N)| of |N| categorical map providing the |N| fuzzy
%     memberships of 2D points.
%
%% Properties [propertyname  propertyvalues]
% *|'thmax'|* : threshold on the membership value, see part 1 of the
%     algorithm described below.
%
% *|'thcon'|* : threshold on the confusion index, see part 2.
%
% *|'thent'|* : entropy threshold, see part 3.
%
%% Output
% *|fmap|* : map of fuzzy boundaries.
%
%% Algorithm
% # an individual location |x| is selected as belonging to a fuzzy boundary 
% if the value of its maximum fuzzy membership |p_max| is less than a
% threshold value |thmax| [Zhang96], 
% # the confusion index for defining fuzzy boundaries, involves two fuzzy
% membership values for each location: the confusion index is evaluated by
% 1 minus the difference between the fuzzy membership values of location |x|
% belonging to the first most likely and the second most likely classes; a
% threshold |T| is applied so that location |x| defines a fuzzy boundary if
% the confusion index is greater than a pre-defined threshold |thcon| [Bur96],
% # the measure of entropy, for defining fuzzy boundaries, which uses the 
% complete fuzzy membership values for each location. Large values indicate
% low accuracy in classification, while small values indicate high accuracy
% in classification. Then, it is logical to assert that boundaries usually
% occur where locations have high degrees of fuzziness, that is, big values
% of entropy, larger than a threshold |thent| [Foo95].
%
%% References
% [Foo95]  G.M. Foody: "Cross-entropy for the evaluation of the accuracy 
%      of a fuzzy land cover classification with fuzzy ground data", ISPRS
%      Journal of Photogrammetry and Remote Sensing, 50:2-12, 1995.
%      <http://www.sciencedirect.com/science/article/pii/092427169590116V>
%
% [Zhang96]  J.X. Zhang: "A surface-based approcah to handling uncertainties
%      in an urban-orientated spatial database", PhD Thesis, University of
%      Edinburgh, 1996.
%      See also: "A surface approach to the handling of uncertainties in an
%      integrated spatial database environment", Cartographica: International
%      Journal for Geographic Information and Geovisualization, 33(1):23-31,
%      1996.
%      <http://utpjournals.metapress.com/content/y6036594410t1n3t/>
%
% [ZK97]  J.X. Zhang and R.P. Kirby: "An evaluation of fuzzy approaches to
%      mapping land cover from aerial photographs", ISPRS Journal of 
%      Photogrammetry and Remote Sensing, 52(5):193-201, 1997.  
%      <http://www.sciencedirect.com/science/article/pii/S0924271697000269>
%
% [Bur96]  P.A. Burrough: "Natural objects with indeterminate boundaries", 
%      in Geographic Objects with Indeterminate Boundaries, P.A. Burrough and
%      A.U. Frank eds., Taylor & Francis (London, UK), pp.3-28, 1996.
%      <http://www.mendeley.com/research/natural-objects-indeterminate-boundaries/>
%
% [ZK00]  J. Zhang and R.P. Kirby: "A comparison of alternative criteria
%      for defining fuzzy boundaries on fuzzy categorical maps", Geo-spatial
%      Information Science, 3(2):26-34, 2000.
%      <http://www.springerlink.com/content/p0k58753805p2101/>
% 
%% See also
% Called: 
% <FUZZBOUND_BASE.html |FUZZBOUND_BASE|>.

%% Function implementation
function fmap = fuzzbound(fmemb,varargin)

%%
% parsing parameters

error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(1, 1, nargout, 'struct'));

% mandatory parameter
if ~isnumeric(fmemb)
    error('fuzzbound:inputerror','a matrix is required in input'); 
end

% optional parameters
p = createParser('FUZZBOUND');   
% principal optional parameters
p.addParamValue('thmax', [], @(x)isempty(x) || (isscalar(x) && x>=0. && x<=1));
p.addParamValue('thcon', [], @(x)isempty(x) || (isscalar(x) && x>=0. && x<=1));
p.addParamValue('thent', [], @(x)isempty(x) || (isscalar(x) && x>=0. && x<=1));

% parse and validate all input arguments
p.parse(varargin{:}); 
p = getvarParser(p);        

%%
% main calculation

fmap = fuzzbound_base(fmemb, p.thmax, p.thcon, p.thent);

end % end of fuzzbound
