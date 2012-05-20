function [d,dcum]=diffsoleil(a,b,varargin)
%DIFFSOLEIL calculates the difference between B matching A using scattergramsoleil
%   Syntax: d = diffsoleil(A,B [,'property1',value1])
%       A,B: intensity images
%   Pair properties related to anisodiff2D (tangential smoothing)
%             num_iter:  2
%              delta_t:  1/7
%                kappa:  10
%               option:  2
%   Pair properties for rescaling
%                depth: 256
%              prctile: [.1 99.9]
%
%
% Reference: 
% P. A. Bromiley and N. A. Thacker and P. Courtney, Non-parametric image subtraction using grey level scattergrams, Image and Vision Computing, 2002, 20, 9-10.
%
% Example:2
% sample = samplesoleil(db,sampleengine(20)); % Series 20
% topleftposition = round(64/1024.*[430 430]); % x,y
% dimension       = round(64/1024.*2.^nextpow2([240 230])); % width, height
% imval2 = sample.im; imval2 = imval2(topleftposition(2)-1+(1:dimension(2)),topleftposition(1)-1+(1:dimension(1)),:,:); % restrict to a region of interest
% imsiz=size(imval2); imval2 = reshape(imval2,[prod(imsiz(1:2)) imsiz(3)]);
% % [a,b]=diffsoleil(reshape(imval2(:,1),64,64),reshape(imval2(:,300),64,64))
% figure, imagesc(diffsoleil(reshape(imval2(:,1),imsiz(1),imsiz(2)),reshape(imval2(:,300),imsiz(1),imsiz(2))))
%
%
%   See also: anisodiff2D scattergramsoleil

% Soleil experiments SUN2011c - 15/02/12 - INRA\Olivier Vitrac - rev. 27/02/12

% Revision history
% 27/02/11 rename percentile as prctile

%% default
default_aniso2D = struct(...
    'num_iter',2,...
    'delta_t',1/7,...
    'kappa',10,...
    'option',2 ...
    );
default_rescale = struct(...
    'depth',256,...
    'filtzero',20,...
    'prctile',[.1 99.9] ...
    );

% argcheck
if nargin<2, error('2 arguments are required'), end
siz = size(a);
if ~matcmp(siz,size(b)), error('A and B must be of the same size'); end
[options_rescale,remaining] = argcheck(varargin,default_rescale);
if (options_rescale.depth<256 && options_rescale.depth>4000), error('depth must be ranged between 64 and 4000'); end
if length(options_rescale.prctile)~=2, error('percentile must be a 1x2 vector'); end
options_rescale.depth = round(options_rescale.depth);
options_rescale.filtzero = round(min(options_rescale.depth/4,options_rescale.filtzero));
options_aniso2D = argcheck(remaining,default_aniso2D);

%% rescaling
options_rescale.prctile = sort(options_rescale.prctile);
range = prctile([a(:);b(:)],options_rescale.prctile);
if options_rescale.depth<=256
    a8 = uint8( (options_rescale.depth-1) * (a-range(1))/(range(2)-range(1)) );
    b8 = uint8( (options_rescale.depth-1) * (b-range(1))/(range(2)-range(1)) );
else
    a8 = uint16( (options_rescale.depth-1) * (a-range(1))/(range(2)-range(1)) );
    b8 = uint16( (options_rescale.depth-1) * (b-range(1))/(range(2)-range(1)) );
end

%% scattergram
s = scattergramsoleil(a8,b8,options_aniso2D,'depth',options_rescale.depth,'printon',false);

%% difference
scum = zeros(options_rescale.depth,1,class(a));
d = zeros(siz,class(a));
[~,gbmax] = max(s);
gbmax_likely = min((1:options_rescale.depth),gbmax);
gbmax_smooth = min(gbmax_likely,filtzero(gbmax_likely,options_rescale.filtzero));
for ga=1:options_rescale.depth % gray level
    %   cumulative difference
    u = (max(ga,gbmax(ga)):options_rescale.depth)'; % gray index matching gray range
    if length(u)>1, scum(ga) = trapz(u,s(u,ga)); end
    % local difference
    ua = find(a8==ga); % position index matching gray range
    d(ua) = max(0,b8(ua)-gbmax_likely(ga)); %  max(0,b8(ua)-min(ga,gbmax(ga)));
end

%% outputs
d = d * (range(2)-range(1))/(options_rescale.depth-1);
if nargout>1
    dcum = scum(a8+1) * (range(2)-range(1))/(options_rescale.depth-1);
end