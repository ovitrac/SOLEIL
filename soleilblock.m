function I = soleilblock(rgb,blocks,varargin)
%SOLEILBLOCK read and assemble in rectangular rgb images
%   syntax: I = soleilblock(rgb,blocks[,'property1',value1,'property2',value2,...])
%
% INPUTS:
%          ind: cell array frames to display
%       blocks: block definition as used in LEGBLOCKIMAGE ([] = sub-block)
%               example: {2 [1;1]}
%       method: method used to rescale images
%          sep: structure coding for horizontal and vertical separators between images
%               default = struct('x',0,'y',0)
%         sepx: horizontal separator (default = 0, used if "sep" is not defined)
%         sepy: vertical separator (default = 0, used if "sep" is not defined)
%
% OUTPUT:
%            I: merged image according to blocks
%
% EXAMPLE: I=soleilblock(rgb,{[1 1;1 1]},'ind',[2 3 1]); figure, imshow(I)
%
%
%   see also: LEICAREADBLOCK, LEGBLOCKIMAGE


% CONFOCAL 1.1 - 14/03/09 - INRA/Olivier - rev.

%revision history
% 13/03/12 first release adapted from leicareadblock
%          if elements to display < elements to merge, empty frames are added

% definitions
o_default = struct(...
    'method','bicubic',...
    'colsep',zeros(1,1,3,'uint8'),...
    'sepx'  ,0,... %pixels
    'sepy'  ,0,... %pixels
    't0'    ,clock,...
    'sep'   ,[],...
    'ind'   ,[]...
    );

%arg check
o = argcheck(varargin,o_default,'case');
if nargin<2, error('two arguments are required (series and blocks)'); end
if ndims(rgb)<3 || size(rgb,3)~=3, error('input must be a rgb image or a series of rgb image coded as a 4-D array'), end
if ~iscell(blocks), error('blocks must be a cell'), end
[frameheight,framewidth,~,nframes] = size(rgb);
if isempty(o.ind), o.ind = 1:nframes; end
n = length(o.ind);
if isempty(o.sep), o.sep = struct('x',o.sepx,'y',o.sepy); end

disp('SOLEILBLOCK...')

% check all sub-blocks and the whole assembly (they must be rectangle)
nsubB  = 0; %# number of subblocks (also rectangular)
sizB   = size(blocks);
widthB  = zeros(sizB(1),1);
heightB = zeros(sizB(2),1);
for i = 1:sizB(1)
    for j = 1:sizB(2)
        w = sum(blocks{i,j},2);
        h = sum(blocks{i,j},1);
        nsubB = nsubB + numel(blocks{i,j});
        if (~all(w==w(1))) || (~all(h==h(1))), error('the subblock blocks{%d,%d} must be rectangle',i,j), end
        heightB(j) = heightB(j) + h(1);
        widthB(i)  = widthB(i)  + w(1);
    end
end
if (~all(widthB==widthB(1))) || (~all(heightB==heightB(1))), error('all blocks must be rectangle'), end
heightB = heightB(1);
widthB  = widthB(1);
if nsubB<n, error('the number of defined blocks (%d) is lower than the number of images (%d)',nsubB,n); end
% if elements to display < elements to merge, empty frames are added
if nsubB>n,
    warning('Elements to display do not fill a rectangle: parts of the merged image are displayed in white')
    rgb(:,:,:,end+1) = 256*ones([frameheight,framewidth,3]); 
    if size(o.ind,1) == 1, o.ind = [o.ind nframes+1*ones(1,nsubB-n)];
    elseif size(o.ind,2) == 1, o.ind = [o.ind;nframes+1*ones(1,nsubB-n)'];
    end
    n=length(o.ind);
end


% image info
disp('     check all images')
Heightmax = frameheight;
Widthmax = framewidth;

% define the size of each image
imsiz = repmat(struct('xmag',[],'ymag',[],'Height',[],'Width',[]),n,1);
assembly = cell(sizB);
iB = 0;
for i = 1:sizB(1)
    for j = 1:sizB(2)
        for is = 1:size(blocks{i,j},1)
            for js = 1:size(blocks{i,j},2)
                iB=iB+1;
                if iB<=n
                    imsiz(iB).xmag   = blocks{i,j}(is,js); %/widthB;
                    imsiz(iB).ymag   = blocks{i,j}(is,js); %/heightB;
                    imsiz(iB).Width  = round(imsiz(iB).xmag * Widthmax);
                    imsiz(iB).Height = round(imsiz(iB).ymag * Heightmax);
                    assembly{i,j}(is,js) = iB;
                end
            end
        end
    end
end

% check whether images are not truncated
xmagmin   = min([imsiz.xmag]);
ymagmin   = min([imsiz.ymag]);
Widthmin  = min([imsiz.Width]);
Heightmin = min([imsiz.Height]);
for i=1:n
    imsiz(i).Width = round(imsiz(i).Width * (imsiz(i).xmag/xmagmin)/(imsiz(i).Width/Widthmin));
    imsiz(i).Height = round(imsiz(i).Height * (imsiz(i).ymag/ymagmin)/(imsiz(i).Height/Heightmin));
end

% preallocating (for speed)
disp('     preallocate')
I = repmat({repmat(o.colsep,[Heightmax+2*o.sep.x,Widthmax+2*o.sep.y,1])},n,1); % uint8

% load all images
for i=1:n % for each image
    requestedHeight = imsiz(i).Height;
    requestedWidth = imsiz(i).Width;
    imtmp = rgb(:,:,:,o.ind(i));
    if (frameheight~=Heightmax) || (framewidth~=Widthmax)
        disp(sprintf('     ->update the size of this frame (method=''%s'')',o.method))
        I{i}(o.sep.x+1:end-o.sep.x,o.sep.y+1:end-o.sep.y,:) = imresize(imtmp,[Heightmax,Widthmax],o.method);
    else
        I{i}(o.sep.x+1:end-o.sep.x,o.sep.y+1:end-o.sep.y,:) = imtmp;
    end
    if (size(I{i},1)~=requestedHeight+2*o.sep.y) || (size(I{i},2)~=requestedWidth+2*o.sep.x)
        disp(sprintf('     ->resize all channels (method=''%s'')',o.method))
        I{i} = imresize(I{i},[requestedHeight+round(2*o.sep.y*requestedHeight/(size(I{i},1)-2*o.sep.y)),...
                              requestedWidth+round(2*o.sep.x*requestedWidth/(size(I{i},2)-2*o.sep.x))],...
                        o.method);
    end
end

% assemble all images into a single image according to blocks specifications
disp('     assembling')
J = cell(sizB);
for i = 1:sizB(1)
    for j = 1:sizB(2)
        tmp = cell(size(blocks{i,j}));
        tmp(1:numel(assembly{i,j})) = I(assembly{i,j}(:));
        J{i,j} = cell2mat(tmp);
    end
end
I = cell(sizB(1),1);
for i = 1:sizB(1)
    I{i} = cell2mat(J(i,:));
end
I = cell2mat(I);
disp(sprintf('SOLEILBLOCK end in %0.4g s',etime(clock,o.t0)))