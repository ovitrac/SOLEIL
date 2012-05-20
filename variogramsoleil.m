function [braw,b2,start] = variogramsoleil(a,varargin)
% VARIOGRAMSOLEIL generates a sampling variogram of time-series on uint8 frames
%   Syntax: b = variogramsoleil(a [,'property1',value1','property2','value2',...])
%           [b,b2,start] = variogramsoleil(...)
%   INPUTS
%        a: nrows x nwidth x nframes array of class uint8 (nframes > 5) 
%   anisodiff2D properties/values (tangential smoothing): to be used ONLY for calculating b2 and start
%             num_iter:  10
%              delta_t:  1/7
%                kappa:  2
%               option:  2
%                 mask:  [] (indexes of a)
%              printon:  true
%
% OUTPUTS
%     braw: nframesx256 array variogram with values ranged between 0 and 1 (not subjected to anisodiff2D since it use onappropriated boundary conditions)
%       b2: nframesx256 second derivative of b (based on anisodiff2D)
%    start: nframesx1 starting index of the right distribution
%
% EXAMPLE 1
%   db = loaddbsoleil;
%   out = moviesoleil(db,'series',51,'ind',190:300,'background',[],'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil','videoon',false,'noplot','lowmemory');
%   [v,v2,vstart] = variogramsoleil(permute(out.rgb3(:,:,1,:),[1 2 4 3]));
%   figure, imagesc(v(:,2:end)), figure, plot(vstart), figure, plot(min(vstart)-1:255,v(:,min(vstart):end)'), xlabel('intensity'), ylabel('frequency')
%   figure, plot(0:size(v,1)-1,sum(v(:,min(vstart):end),2)), xlabel('frame index'), ylabel('filling fraction')
%
% EXAMPLE 2
%   out2 = moviesoleil(db,'series',20,'ind',110:245,'background',[],'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil','videoon',false,'noplot','lowmemory');
%   [v,v2,vstart] = variogramsoleil(permute(out2.rgb3(:,:,1,:),[1 2 4 3]));
%   figure, imagesc(v(:,2:end)), figure, plot(vstart), figure, plot(min(vstart)-1:255,v(:,min(vstart):end)'), xlabel('intensity'), ylabel('frequency')
%   figure, plot(0:size(v,1)-1,sum(v(:,min(vstart):end),2)), xlabel('frame index'), ylabel('filling fraction')
%
% See also: scattergramsoleil loaddbsoleil moviesoleil anisodiff2D diffsoleil

% Soleil experiments SUN2011c - 26/03/12 - INRA\Olivier Vitrac and JMV - rev. 26/03/12

% Revision history
% 26/03/12 RC by OV: implement v2, start - do not apply anisodiff2D to v data, finalize example1, add example 2
% 09/04/12 by JMV: modify the definition of ibest by taking min(length(is,5)) instead of 5 (some cases have lower length than 5)
% 19/04/12 add mask

% default
default = struct(...
    'num_iter',10,...
    'delta_t',1/7,...
    'kappa',2,...
    'option',2,...
    'depth',256,...
    'mask',[],...
    'printon',true ...
    );
depth = 256;
scale = (0:depth-1)';

% argcheck
if nargin<1, error('one argument is required'), end
if ~strcmp(class(a),'uint8'), error('input must be uint8 class'), end
[height,width,nframes] = size(a);
if ndims(a)<3 || ndims(a)>3, error('a must be ny x nx x nt array'), end
if nframes <=5, error('more than 5 frames are required'), end
o = argcheck(varargin,default);
if ~isempty(o.mask), nmask = length(o.mask); end

% lookup frequencies for each frame
% loop for on time
braw = zeros(nframes,depth); % depth = 256
screen = ''; t0 = clock;
for iframe = 1:nframes
    if o.printon, screen = dispb(screen,'VARIOGRAMSOLEIL samples %d/%d',iframe,nframes); end
    tmp = a(:,:,iframe);
    if isempty(o.mask)
        braw(iframe,:) = histc(double(tmp(:)),scale)/(height*width);
    else
        braw(iframe,:) = histc(double(tmp(o.mask)),scale)/nmask;
    end
end
if o.printon, dispb(screen,'VAROGRAMSOLEIL saming completed in %0.4g s, tangential filtering...',etime(clock,t0)); end


% 2nd derivative of filteres braw values
if nargout>1
    % anisodiff2D (same default parameter than scattergramsoleil.m)
    b = anisodiff2D(braw,o.num_iter,o.delta_t,o.kappa,o.option,o.printon);
    % second order derivative in intensity (not in time)
    b2 = ndf(scale,b',2)';
    if nargout>2 % find begining of the second distribution (largest increasing peak with positive start)
        start = zeros(nframes,1);
        for iframe = 1:nframes
            [p,~,a]=monotone(b2(iframe,:),'+');
            [~,is] = sort(a./abs(b2(iframe,p))','descend'); % curvature normalized by its primitive
            %[~,is] = max(a./abs(b2(iframe,p))'); % for efficiency instead of sort
            [~,ibest] = max(a(is(1:min(5,length(is))))); % search through the first five (more robust) instead of taking the first
            start(iframe) = p(is(ibest)); % largest curvature and peak
        end
    end
end