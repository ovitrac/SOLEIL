function hout = legblockimage(blocks,hax,obj,sep,linewidth,vshift)
%LEGBLOCKIMAGE legend an image created by assembling rectangular blocks
%   syntax: legblockimage(blocks [,hax,obj,sep,linewith,vshift])
%           hp = legblockimage(...)
%
%       blocks: mxn cell array coding for rectangle images (see example)
%               Axes are assumed to image (ydir = 'reverse')
%           very simple example:
%           A: one image, B: another image with the same size as A, C = imresize(A,2) or C = [A A; A A];
%           imshow([[A;B] C])
%           corresponding blocks is {[1;1] 2}
%
%       hax: current axes (default = [])
%
%       obj: object to plot (implemented objects are: scale, text and frame)
%           examples of object constructors
%             scale object: struct('line',100,'text','\color[rgb]{1 0 1}100 µm','position',2,'color',[1 1 1])
%              text object: struct('text',{'\color[rgb]{1 0 1}a' '\color[rgb]{1 0 1}b' '\color[rgb]{1 0 1}c' '\color[rgb]{1 0 1}d'}
%             frame object: struct('frame',[10 10 200 200],'ind',1,'color',[1 1 1]) --> 200x200 frame
%           Comments:
%           color and position can be omitted
%           use TEX commands to change text properties
%
%           position: 0 (center), 1 (north east), 2 (north west), 3 (south west), 4 (south east)
%
%       sep: distance to edge (value between 0 and 1) (default = 0.2)
%
%       linewidth: default = 3
%
%       vshift: vertical text shift (in pixels) (default = 10)
%
%       hp: handles of plotted objects

% CONFOCAL 1.1 - 10/03/09 - INRA/Olivier - rev. 13/03/09

%revision history
%13/03/09 add frame object

% Definitions
pos_default = 2;
color_default = [1 1 1];
sep_default = 0.2;
linewidth_default = 3;
vshift_default = 10;

%arg check
if nargin<2, hax = []; end
if nargin<3, obj = []; end
if nargin<4, sep = []; end
if nargin<5, linewidth = []; end
if nargin<6, vshift = []; end
if isempty(hax), hax=gca; end
if isempty(sep), sep = sep_default; end
if isempty(linewidth), linewidth= linewidth_default; end
if isempty(vshift), vshift = vshift_default; end
if ~iscell(blocks), error('blocks must be a cell'), end
nsubB  = 0;
sizB   = size(blocks);
widthB  = zeros(sizB(1),1);
heightB = zeros(sizB(2),1);
for i = 1:sizB(1)
    for j = 1:sizB(2)
        w = sum(blocks{i,j},2);
        h = sum(blocks{i,j},1);
        nsubB = nsubB + numel(blocks{i,j});
        if (~all(w==w(1))) || (~all(h==h(1))), error('the sublock blocks{%d,%d} must be rectangle',i,j), end
        heightB(j) = heightB(j) + h(1);
        widthB(i)  = widthB(i)  + w(1);
    end
end
if (~all(widthB==widthB(1))) || (~all(heightB==heightB(1))), error('the whole blocks must be rectangle'), end
heightB = heightB(1);
widthB  = widthB(1);
xax = get(hax,'xlim'); dxax = diff(xax);
yax = get(hax,'ylim'); dyax = diff(yax);

% find position of each image
xpos = zeros(nsubB,1);
ypos = zeros(nsubB,1);
width = zeros(nsubB,1);
height = zeros(nsubB,1);
iB = 0;
yB = yax(1);
for i = 1:sizB(1)
    xB = xax(1);
    for j = 1:sizB(2)
        y = yB;
        for is = 1:size(blocks{i,j},1)
            x = xB;
            for js = 1:size(blocks{i,j},2)
                iB=iB+1;
                xpos(iB)   = x;
                ypos(iB)   = y;
                width(iB)  = dxax * blocks{i,j}(is,js)/widthB;
                height(iB) = dyax * blocks{i,j}(is,js)/heightB;
                x = x + width(iB);
            end
            y = y + height(iB);
        end
        xB = x;
    end
    yB = y;
end


% text (to be modified)
hp = [];
if isempty(obj)
    for iB = 1:nsubB
        hp(iB,1)=text(xpos(iB)+width(iB)/2,ypos(iB)+height(iB)/2,sprintf('\\color[rgb]{1 1 1}%d',iB),'fontsize',48,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
elseif isstruct(obj) && isfield(obj,'text')
    for iB = 1:nsubB
        iobj = min(length(obj),iB);
        if isfield(obj(iobj),'position') && ~isempty(obj(iobj).position), pos = obj(iobj).position; else pos=pos_default;  end
        if isfield(obj(iobj),'color') && ~isempty(obj(iobj).color), col = obj(iobj).color; else col=color_default;  end
        switch pos
            case 1, xposleg = xpos(iB)+(1-sep)*width(iB); yposleg = ypos(iB)+sep*height(iB); xdirection = -1; vpos = 'top'; vdec = vshift;
            case 2, xposleg = xpos(iB)+sep*width(iB); yposleg = ypos(iB)+sep*height(iB); xdirection = +1; vpos = 'top'; vdec = vshift;
            case 3, xposleg = xpos(iB)+sep*width(iB); yposleg = ypos(iB)+(1-sep)*height(iB); xdirection = +1; vpos = 'bottom'; vdec = -vshift;
            case 4, xposleg = xpos(iB)+(1-sep)*width(iB); yposleg = ypos(iB)+(1-sep)*height(iB); xdirection = -1; vpos = 'bottom'; vdec = -vshift;
            otherwise, xposleg = xpos(iB)+width(iB)/2; yposleg = ypos(iB)+height(iB)/2; xdirection = 0; vpos = 'middle';  vdec = 0;
        end
        if isfield(obj,'line')
            hp(iB,1)=text(xposleg+xdirection*obj(iobj).line/2,yposleg+vdec,obj(iobj).text,'horizontalalignment','center','verticalalignment',vpos);
            hp(iB,2)=line(xposleg+[0 xdirection]*obj(iobj).line,yposleg * [1 1],'linestyle','-','linewidth',linewidth,'color',col);
        else
            hp(iB,1)=text(xposleg,yposleg,obj(iobj).text,'horizontalalignment','center','verticalalignment',vpos);
        end
    end 
elseif isstruct(obj) && isfield(obj,'frame')
    for iobj = 1:length(obj)
        if isfield(obj(iobj),'color') && ~isempty(obj(iobj).color), col = obj(iobj).color; else col=color_default;  end
        if ~isfield(obj(iobj),'ind') || isempty(obj(iobj).ind), obj(iobj).ind = iobj; end
        if obj(iobj).ind>nsubB, error('the image indexed %d does not exist',obj(iobj).ind), end
        iB = obj(iobj).ind;
        hp(iobj) = line(  xpos(iB)+obj(iobj).frame(1)+[0 0 1 1 0]'*obj(iobj).frame(3),...
                          ypos(iB)+obj(iobj).frame(2)+[0 1 1 0 0]'*obj(iobj).frame(4),...
                          'linestyle','-','linewidth',linewidth,'color',col);
    end
end

if nargout, hout = hp; end