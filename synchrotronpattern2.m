function [background,Ib,threshout] = synchrotronpattern2(I,varargin)
%SYNCHROTRONPATTERN retrieves the synchrotron pattern from an image series using a wavelet decomposition
%   background = synchrotronpattern(I [,'property1',value1,'property2,value2,...])
%   [background,Ib,thresh] = synchrotronpattern(...)
%
%   INPUTS
%       I: heightxwidthxn image array
%       List of available pair properties:
%           background: frames taken into account to get the background (default = [])
%                      empty value mean all frames
%               class: output of classsoleil (default=[])
%             details: 'v', 'h' or 'd' for vertical, horizontal and diagonal details respectively
%                      (see detcoef2 for details)
%        detailslevel: wavelet level decomposition (see wmaxlev for details)
%             iframes: indices of frames to be normalized (default = [])
%                      empty value mean all frames
%              interp: interpolation method for coefficients (default = 'lanczos3'), see imresize
%                      available methods: 'nearest', 'box', 'triangle', 'cubic', 'lanczos2', 'lanczos3'
%               level: wavelet level decomposition (see wmaxlev for details)
%          medianfilt: median filter width (default = []), set to 1, [] or NaN to remove median filtering
%             prctile: threshold percentiles for uint8 conversion (default = [0.1 99.9])
%             wavelet: wavelet type (default = 'db1')
%                      TIP: wavemngr('read',1) list available wavelets and their shorthands to be used
%
%           Keyword 'noscaling' force results not to be scaled
%
%   OUTPUTS
%        background: heightxwidth background image
%                Ib: heightxwidthxn image array as I but removed background
%            thresh: nx2 array of threshold values (for each image)

%CONFOCAL 2.0 - 01/01/12 - INRA/Olivier Vitrac - rev. 10/03/12

% Revision history
% 02/01/12 release candidate
% 03/01/12 threshold calculated but not applied and uint8 removed : user can change it out of this script
% 03/01/12 add background, ind
% 04/01/12 rename background as reference and ind as iframes
% 12/01/12 add class, add noscaling, rename thres as thresh
% 19/02/12 fix permute height and width
% 23/02/12 approximated image and original one have the same mean
% 26/02/12 update help
% 10/03/12 fix recursion to account for update in argcheck (first defined) - see archeck, rev. 28/02/12

% default
default = struct(...
    'iframes',[],... index of images to keep
    'background',[],... index of images used only for background determination
    'class',[],...
    'wavelet','db1',...
    'level',5,...
    'interp','lanczos3',...
    'prctile',[0.1 99.9],...
    'medianfilt',[],...
    'details','',...
    'detailslevel',[] ...
    );
kwlist = 'noscaling';

% argcheck
if nargin<1, error('one argument is required'); end
if ~ismember(ndims(I),[2 3]), error('only intensity images or image series are accepted'), end
options = argcheck(varargin,default,kwlist);
if ~isa(I,'double') || ~isa(I,'single'), I = single(I); end

% check parameter values
[height,width,n] = size(I);
if length(options.level)>1, error('level must be a scalar'), end
maxlevel = wmaxlev([height width],options.wavelet);
if options.level>maxlevel
    dispf('WARNING:\t the level value (%d) exceeds the maximum level value (%d) for decomposing an image %dx%d with ''%s'' wavelets',...
        options.level,maxlevel,height,width,options.wavelet);
    options.level = min(maxlevel,options.level);
end
if length(options.prctile)<2, options.prctile = [options.prctile 100-options.prctile]; end
options.prctile = [min(options.prctile) max(options.prctile)];
medianfilton = ~isempty(options.medianfilt) && ~any(isnan(options.medianfilt(1))) && options.medianfilt(1)~=1;
if medianfilton, options.medianfilt = options.medianfilt(1); end
detailson = ~isempty(options.details);
if detailson 
    if ~ischar(options.details), error('details must be a string'); end
    if ~ismember(options.details,{'v' 'h' 'd'}), error('the detail value ''%s'' is incorrect (accepted values = v, h or d, see detcoef2)',options.details), end
end
if isempty(options.detailslevel), options.detailslevel = options.level; end
options.detailslevel = min(options.detailslevel,options.level(end));
options.detailslevel = options.detailslevel(:)';
if isempty(options.background), options.background = 1:n; end
if isempty(options.iframes), options.iframes = 1:n; end

%% recursion
if ~isempty(options.class)
    classes = unique(options.class);
    nclasses = length(classes);
    [options.iframes,jframes] = sort(options.iframes);
    nframes = numel(options.iframes);
    Ib0 = zeros([height,width,nframes],class(I));
    thresh = zeros(nframes,length(options.prctile));
    for iclass = 1:nclasses
        isclass = options.class==classes(iclass); % true it belongs to the good class
        jcurrentframes = jframes(isclass(options.iframes));
        dispf('----> study based on class %d/%d',iclass,nclasses)
        [background,Ib0(:,:,jcurrentframes),thresh(jcurrentframes,:)] = synchrotronpattern2(I(:,:,jcurrentframes),'class',[],'background',[],'iframes',[],'noscaling',options);
    end
    thresmin = min(thresh(:,1));
    thresmax = max(thresh(:,end));
    if ~options.noscaling
        Ib0 = 255 * (Ib0 - thresmin) / (thresmax - thresmin);
    end
    if nargout>1, Ib=Ib0; end
    if nargout>2, threshout = thresh; end
    return
end

%% main

% time counter
t0 = clock;
screen = '';

% stack average if required
if n>1
    screen = dispb(screen,'SYNCHROTRON PATTERN: stack averaging between frames #%d and #%d...',min(options.background),max(options.background));
    Im = mean(I(:,:,options.background),3);
else
    Im = I;
end

% Downsampling using wavelets (common to all background extraction methods)
screen = dispb(screen,'SYNCHROTRON PATTERN: background extraction...');
%if ~isa(Im,'double'), Im = double(Im); end
[C,S] = wavedec2(Im,options.level,options.wavelet);
A = appcoef2(C,S,options.wavelet); % do this: A = reshape(C(1:S(1,1)*S(1,2)),S(1,:));
A = mean(Im(:)) * A / mean(A(:)); % approximated image and original one must have the same mean

% Background based on wavelets
background = imresize(A,[height, width],options.interp);
if detailson
    dtmp = zeros([height, width],class(background));
    for dlevel = options.detailslevel
        dtmp = dtmp + imresize(detcoef2(options.details,C,S,dlevel),[height, width],options.interp);
    end
    d = (dtmp+background)./background;
    background = min(background,background.*d);
end

% Background removal (if asked)
if nargout>1
    nframes = numel(options.iframes);
    Ib = zeros([height,width,nframes],class(I)); %Ib = zeros([width,height,n],'uint8');
    thresh = zeros(nframes,length(options.prctile));
    for i=1:nframes
        t1=clock;
        if medianfilton, 
            tmp = medfilt2(I(:,:,options.iframes(i)),[options.medianfilt options.medianfilt])./background;
        else
            tmp = I(:,:,i)./background; % not options.iframes(i) as I contains nframes frames
        end
        Ib(:,:,i) = tmp;
        thresh(i,:) = prctile(tmp(:),options.prctile);
        screen = dispb(screen,'SYNCHROTRON PATTERN: background removed in frame %d/%d in %0.4g s...',i,nframes,etime(clock,t1));
    end
    thresmin = min(thresh(:,1));
    thresmax = max(thresh(:,end));
    if ~options.noscaling
        Ib = 255 * (Ib - thresmin) / (thresmax - thresmin);
    end
end

% third output
if nargout>2, threshout = thresh; end

dispb(screen,'SYNCHROTRON PATTERN completed for %d images in %0.5g s',n,etime(clock,t0));


