function [background,Ib,thres] = synchrotronpattern(I,p,f,q)
%SYNCHROTRONPATTERN retrieve the synchrotron pattern from an image series
%   [background,Ib,thres] = synchrotronpattern(I,p,f,q)
%       I: heightxwidthxn image array
%       p: csaps regularizing parameter (default = 1e-5)
%       f: f or [fy fx] median filter (default = 5)
%       q: q, [q 100-q], [qmin qmax] threshold to match 0-255 (default = 0.1)


%CONFOCAL 2.0 - 16/11/11 - INRA/Olivier Vitrac - rev.

% default
p_default = 1e-5;
f_default = 5;
q_default = 0.1;

% arg check
if ~ismember(ndims(I),[2 3]), error('only intensity images or image series are accepted'), end
if nargin<2, p = []; end
if nargin<3, f = []; end
if nargin<4, q = []; end
if isempty(p), p = p_default; end
if isempty(f), f = f_default; end
if isempty(q), q = q_default; end
if length(q)<2, q = [q 100-q]; end
if all(q<1), q = q*100; end
q = sort(q(1:2),'ascend');
[width,height,n] = size(I);
if ~isa(I,'double') || ~isa(I,'single'), I = single(I); end
screen = '';

% background extraction
t0 = clock;
x = {1:height 1:width};
screen = dispb(screen,'SYNCHROTRON PATTERN: initialization...');
Im = mean(I,3);
if ~isa(Im,'double'), Im = double(Im); end
screen = dispb(screen,'SYNCHROTRON PATTERN: background smoothing...');
background = csaps(x,Im,p,x);

% background removal
if nargout>1
    Ib = zeros([width,height,n],class(I));
    for i=1:n
        t1=clock;
        Ib(:,:,i) = medfilt2(I(:,:,i),[f f])./background;
        screen = dispb(screen,'SYNCHROTRON PATTERN: background removed in image %d/%d in %0.4g s...',i,n,etime(clock,t1));
    end
end

% statistical analyis
if nargout>2
    thres = zeros(n,2);
    screen = dispb(screen,'synchrotron pattern: background smoothing...');
    for i=1:n
        tmp = Ib(:,:,i);
        thres(i,:) = prctile(tmp(:),q);
    end
    thresmin = min(thres(:,1));
    thresmax = max(thres(:,2));
    Ib = 255 * (Ib - thresmin) / (thresmax - thresmin);
end

dispb(screen,'SYNCHROTRON PATTERN completed for %d images in %0.5g s',n,etime(clock,t0));