function [hframes,hplot,seriesproperties]=montagesoleil(rgb,varargin)
%MONTAGESOLEIL montage alternative for SOLEIL results (diplay selected frames and a global kinetics for the red channel)
%   Syntax: montagesoleil(rgb [,'property1','value1','property2','value2'...])
%           [hframes,hplot] = montagesoleil(...)
%
% INPUTS
%           rgb: Width x Height x 3 x nframes matrix (3 is the number of channels)
%           IMPLEMENTED PAIR/PROPERTIES (property:value)
%                   ind: index of frames to display as images
%              realtime: nframesx1: time scale at zero when oil is added in region of interest, to make kinetic from effective oil deposit on visualized cells
%                indall: all frames chosen for kinetic (not mandatory if realtime is an input)
%                    t0: time when oil is added (not mandatory if realtime is an input)
%           pixellength: width of one pixel (in microns) (default = 1270/1024)
%         scalebarwidth: length of the scale bar relative to the frame width (default = 1/5)
%                        effective length calculated as 10*(round(o.roi(3)*o.bardiv/10)) to a have a multiple of 10
%                series: name of the series (default = 'unknown')
%             textcolor: set the color of the legends (time and scale) (default = [0.8 0.8 0.8])
%             linecolor: set the color of the scale line (default = [0.8 0.8 0.8] to get white color)
%             xsubplots: x input for subplots() (default = [])
%             ysubplots: y input for subplots (default = [])
%        sclararization: anonymous function that returns a scalar for each frame, its outputs are plotted in hplot
%                        (default = @(x) squeeze(mean(mean(x(:,:,1,:),2),1)))
%                  sepx: number of (white) pixels added to one frame for horizontal separation, left AND right (default = 10)
%                  sepy: number of (white) pixels added to one frame for vertical separation, bottom AND up (default = 10)
%               hframes: user override axes handles to plot frames (default = [])
%                 hplot: user override axes handle to plot the kinetics (default = [])
%               polygon: nverticesx2 array or npolygonx1 cell array containing nverticesx2 arrays (default = [])
%                        set polygon(s) (e.g. cells) on to which intensities are calculated
%                  plot: [] (default) displays frames, filling ratio from variogram
%                        'debug'      displays frames, filling ratio from variogram, intensity ratio from percentiles and average
%                        'none'       no display
%
% ADVANCED INPUTS (recusive usage of montagesoleil)
%                 object: nobjectx1 array structure with fields (matching db(iseries).batch(ibatch).object supplied by loaddbsoleil())
%                        ind: index used to display frames, relatively to indall (to set before calling montagesoleil: ind_montagesoleil = ind - indall(1) + 1)
%                   realtime: time scale in seconds corresponding to indexes indall
%                        roi: region of interest used to display cropped frames and to define polygon area
%                    polygon: [xtopleft,ytopleft,width,height] array defining an area of a cell filled by oil
%                             polygon is defined relatively to roi, from the top left corner
%
% for information (to be deleted)
%     series=17; indall = [41:50 52:2:70 75:5:100 110:10:220]; ind =  [41 45 50  54 60 66   70 80 120];
%     out = moviesoleil(db,'series',series,'ind',indall,'background',[],'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil',defaultanalysis,'noplot','lowmemory');
%     [~,~,filename]=montagesoleil(out.rgb3,'ind',find_multiple(ind,indall),'realtime',indall-indall(1)+1,'series',titre(series,indall));
%     print_pdf(600,filename,pwd,'nocheck');
%
%
% OUTPUTS
%               hframes: valid axes handles to display frames (default = [])
%                 hplot: valid axes handle for plotting the intensity kinetic
%      seriesproperties: inputs
%                        seriesfilename: name of the experiment (number of series, frames min and max, method used (invasion or diffsoleil)
%                        I: merged images cropped according to roi and displayed
%                        fillingratio: ratio of cell area filled by oil, using variogram
%                        variogram: frequency of red channel intensity plotted for intensities from 0 to 255 and versus frame index
%                        blocks: montage used in soleilblock to merge frames

%   See also: subplots, montage, soleilblock, legendblock, legblockimage
%
%
%
% EXAMPLE 1:
% test = zeros([256 256 3 6]); for i=1:5, for j = 1:3, test(:,:,j,i)= rand(256); end, end
% montagesoleil(test,'ind',2:2:6,'series','test','realtime',1:2:5,'indred',2:6)
%
% EXAMPLE 2:
% db = loaddbsoleil;
% out = moviesoleil(db,'series',20,'ind',110:245,'background',[],'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil','videoon',false,'noplot','lowmemory');
% montagesoleil(out.rgb3,'ind',[1 41 71 126],'realtime',[510:645],'series',20,'polygon',[561 615 655 589 504 474;493 508 568 628 619 544]')
%
% EXAMPE 3:
% series = 18; indall = 50:99; ind = [88:89 96:99]; background = 38; roi = [590 180 256 256]; matrix = 'CELL layer 1';
% outrgb = zeros([1024 1024 3 length(indall)]);
% out = moviesoleil(db,'series',series,'ind',indall,'background',[],'image','rgb3','roi',[],'synchrotronpatternengine','homomorphicfilersoleil',defaultanalysis,'noplot','rgb_rescale',true);
% outrgb(:,:,2,:) = out.rgb3(:,:,2,:);
% out = moviesoleil(db,'series',18,'image','rgb3','roi',[],'ind',indall,'background',background,'rgb3_num_iter',8,defaultanalysis,'noplot','rgb_rescale',true);
% outrgb(:,:,1,:) = out.rgb3(:,:,1,:);
% outrgbroi = outrgb(roi(2):roi(2)+roi(4)-1,roi(1):roi(1)+roi(3)-1,:,:);
% [~,~,filename]=montagesoleil(uint8(outrgbroi),'ind',ind-indall(1)+1,'realtime',indall-ind(1)+1,'series',sprintf('%s/Series=%d - frame=%d:%d',matrix,series,indall([1 end])));

% 12/03/12 - INRA\Jean-Michael Vauvre - Olivier Vitrac - 10/05/12

% Revision history
% 12/03/12 first release
% 13/03/12 debug scaletot and add "ind" as variable not mandatory [default = size(rgb,4)] for easier use when rgb has already selected frames
% 14/03/12 add surfacered
% 15/03/12 add title, seriesfilename
% 29/03/12 add fillingratio, add example 2
% 16/04/12 add polygon to constructor
% 18/04/12 implement polygon
% 19/04/12 add polygon as mask for variogram, add variable plot
% 07/05/12 add object
% 08/05/12 add indall and t0 as alternative to calculate realtime (calculation outside the function or inside)
% 08/05/12 add matrix to clarify result (matrix used: cell layer or batter)
% 10/05/12 check of the green part (comments), debugging name of series (object -> object(iobject))
% 14/05/12 add sepx and sepy as variables to keep it in output

% default
o_default = struct(...
    'object',[],...
    'series','unknown',...
    'ind',[],...
    'indall',[],...
    't0',[],...
    'realtime',[],...
    'pixellength',1270/1024,...
    'prctile',[25 50 75],...
    'scalebarwidth',1/5,...
    'textcolor',[1 1 1], ...
    'linecolor',[0.8 0.8 0.8],...
    'xsubplots',[],...
    'ysubplots',[],...
    'hframes',[],...
    'hplot',[],...
    'scalarization',@(x) squeeze(mean(mean(x(:,:,1,:),2),1)),...
    'sepx',10,...
    'sepy',10,...
    'polygon',[],...
    'plot',[] ...
);
defaultpaperposition = [0.6345    6.3452   20.3046   15.2284]; % to be updated
batchfields = {'ind' 'indall' 't0' 'roi' 'polygon'}; % required batch fields with objects (see loaddbsoleil for details)
montagefields = setdiff(batchfields,{'indall','t0','roi'});
defaultoutput = struct('seriesfilename','','I',[],'fillingratio',[],'variogram',[],'blocks',[],'inputs',[]);

% argcheck (part 1)
o = argcheck(varargin,o_default,'case','nostructexpand',1);
object = o.object;
o = argcheck(rmfield(o,'object'),o_default,'case');
o.object = object;
if nargin<1, error('one argument is required'); end
if ndims(rgb)<3 || size(rgb,3)~=3, error('input must be a rgb image or a series of rgb image coded as a 4-D array'), end
[frameheight,framewidth,~,nframes] = size(rgb);
numpixels = frameheight*framewidth;
if isnumeric(o.series), o.series = sprintf('Series %02d', o.series); end
if isempty(o.ind), o.ind = 1:nframes; end
if isempty(o.realtime) && isempty(o.object)
    if ~isempty(o.indall) && ~isempty(o.t0), o.realtime = o.indall - o.t0;
    else error('realtime is a mandatory input, use indall and t0 as alternative'),
    end
end
if ~isempty(o.polygon) && ~iscell(o.polygon), o.polygon = {o.polygon}; end

% Recursion
if ~isempty(o.object) && isstruct(o.object)
    if ~all(isfield(o.object,batchfields)), error('object must be built with loaddsoleil()'), end
    nobject = length(o.object);
    oglobal = rmfield(o,[{'object'} montagefields]);
    outtmp = repmat(defaultoutput,nobject,1);
    [hframestmp, hplottmp] = deal(cell(nobject,1));
    for iobject=1:nobject
        if numel(o.object(iobject).roi)~=4, error('bad roi values for object %d/%d (must of length 4)',iobject,nobject), end
        % name of series: definition of the matrix where oil is observed
        if o.object(iobject).matrix == 1, matrix = 'layer1';
        elseif o.object(iobject).matrix == 2, matrix = 'layer2';
        elseif o.object(iobject).matrix == 3, matrix = 'layer2fast';
        elseif o.object(iobject).matrix == 4, matrix = 'batter';
        else matrix = 'undefined matrix';
        end
        [hframestmp{iobject},hplottmp{iobject},outtmp(iobject)] = montagesoleil(....
            rgb(o.object(iobject).roi(2):o.object(iobject).roi(2)+o.object(iobject).roi(4)-1,o.object(iobject).roi(1):o.object(iobject).roi(1)+o.object(iobject).roi(3)-1,:,:),...
            'series',sprintf('%s - %d of %d - %s',o.series,iobject,nobject,matrix),...
            oglobal, ...
            'indall',o.object(iobject).indall,...
            't0',o.object(iobject).t0,...
            'ind',find_multiple(o.object(iobject).ind,o.object(iobject).indall),... % tranlation to fit sampled "indall" values
            'realtime',o.object(iobject).indall-o.object(iobject).t0,... % translation to put time at zero just before oil is added and at one when oil is added
            'polygon',o.object(iobject).polygon ...
            );
    end
    if nargin>0, hframes = hframestmp; end
    if nargin>1, hplot = hplottmp; end
    if nargin>2, seriesproperties = outtmp; end
    return
end

% prefetch outputs
npolygon = max(1,size(o.polygon,1));
if nargout>0, hframes = NaN(npolygon,1); end
if nargout>1, hplot = NaN(npolygon,1); end
if nargout>2, seriesproperties = repmat(defaultoutput,npolygon,1); end

% creation of the block structure to merge displayed images
nframestodisplay = length(o.ind);
if isempty(o.ysubplots) && ~isempty(o.xsubplots)
    if length(o.xsubplots)>1, error('xsubplots must be a scalar in this context'), end
    o.ysubplots = ceil(nframestodisplay/o.xsubplots);
elseif isempty(o.xsubplots) && ~isempty(o.ysubplots)
    if length(o.ysubplots)>1, error('ysubplots must be a scalar in this context'), end
    o.xsubplots = ceil(nframestodisplay/o.ysubplots);
elseif length(o.xsubplots)==1 && length(o.ysubplots)==1 % subplot syntax instead of subplots
    if length(o.xsubplots)>1, error('xsubplots must be a scalar in this context'), end
    if length(o.ysubplots)>1, error('ysubplots must be a scalar in this context'), end
else
    if nframestodisplay <=3, o.xsubplots = 1; else o.xsubplots = ceil(nframestodisplay^(1/3)); end
    o.ysubplots = ceil(nframestodisplay/o.xsubplots);
end
blocks = {repmat(ones(1,o.xsubplots),o.ysubplots,1)};

% creation of the merged image
if ~strcmp('none',o.plot)
    I=soleilblock(rgb(:,:,:,o.ind),blocks,'sepx',o.sepx,'sepy',o.sepy,'colsep',uint8(255*ones(1,1,3)));
end

% Do for each polygon
numprctile = length(o.prctile);

for ipolygon=1:npolygon
    
    % define the indexes of pixels that are located in the cell area (if no cell, default = all pixels)
    if ~isempty(o.polygon)
        [Yim,Xim] = ind2sub([frameheight,framewidth],1:numpixels);
        inpol = find(inpolygon(Xim(:),Yim(:),o.polygon{ipolygon}(:,1),o.polygon{ipolygon}(:,2)));
        ninpol = length(inpol);
    else
        inpol = 1:numpixels;
        ninpol = numpixels;
    end    

    % plot
    if strcmp('debug',o.plot)
        % filling ratio in area delimited by a polygon
        intensitypolygon = zeros(nframes,npolygon);
        for iframe = 1:nframes
            tmp = rgb(:,:,1,iframe);
            intensitypolygon(iframe,ipolygon) = sum(tmp(inpol))/(255*ninpol);
        end

        % Calulation of ratio of pixels above a percentile value (25, 50 and 75)
        surfacered = zeros(length(o.realtime),numprctile);
        tmp = squeeze(double(rgb(:,:,1,:)));
        tmpinpol = zeros(nframes,ninpol);%tmpinpol = zeros(nframes*ninpol,1);
        for i=1:nframes
            tmpbis = tmp(:,:,i);
            tmpinpol(i,:) = tmpbis(inpol);%tmpinpol(1+(i-1)*ninpol:i*ninpol) = tmpbis(inpol);
        end
        ptiles = prctile(tmpinpol(:),o.prctile);
        for i=1:nframes
            tmp = tmpinpol(1+(i-1)*ninpol:i*ninpol);
            for j=1:numprctile
                surfacered(i,j) = length(find(tmp(:)>ptiles(j)))/ninpol;
            end
        end
    end

    % filling ratio
    [v,~,vstart] = variogramsoleil(permute(rgb(:,:,1,:),[1 2 4 3]),'mask',inpol);
    fillingratio = zeros([1 size(v,1)]);
    for i=1:size(v,1),
        for j=min(vstart):255
            tmp = zeros(size(min(vstart):255));
            tmp(j-min(vstart)+1) = v(i,j+1)*(j-min(vstart))/(255-min(vstart));
        end, fillingratio(i) = sum(tmp);
    end

    % plots
    if ~strcmp('none',o.plot)
        % create the figure if required
        if isempty(o.hframes) || isempty(o.hplot)
            hfig = figure;
            formatfig(hfig,'figname',sprintf('%s - frames %d to %d',o.series,o.ind(1),o.ind(end)),'paperposition',defaultpaperposition) 
        else
            hfig = gcf;
        end

        % creation of hframes axes (if required)
        if isempty(o.hframes)
            hplaceholder = subplots([1.2 0.4 1.2 0.3],1,0,0,'alive',1);
            o.hframes = subplots(1,1,0,0,'position',hplaceholder,'strict')'; % row filling
            % delete(o.hframes(nframestodisplay+1:end))
            % o.hframes = o.hframes(1:nframestodisplay);
            delete(hplaceholder)
        elseif ~all(ishandle(o.hframes))
            error('some handles in hframes are not valid')
        end
        % part removed as image is now merged:
        % if length(o.hframes)<nframestodisplay, error('insufficient number of axes (%d) to display %d frames',length(o.hframes),nframestodisplay), end

        % creation of hplot (if required)
        if isempty(o.hplot)
            hplaceholder = subplots([1.5 0.4 1.2 0.3],[0.2 0.6 0.2],0,0,'position',hfig,'alive',8,'strict')';    
            o.hplot = subplots(1,1,0,0,'position',hplaceholder,'strict');
            delete(hplaceholder)
        elseif ~ishandle(o.hplot)
            error('invalid axes handle for hplot')
        end


        % display frames
        subplot(o.hframes)
        imshow(I)
        title(regexprep(o.series,'_','\\_'),'fontsize',12,'fontweight','bold')

        % plot the kinetics
        % % OBSOLETE VERSION
        % % mean intensity vs. time
        % scalarkinetic = o.scalarization(rgb);
        % subplot(o.hplot(1)), hold on, box off
        % plot(o.realtime,scalarkinetic,'k-','linewidth',2)
        % plot(o.realtime(o.ind),scalarkinetic(o.ind),'ro','markerfacecolor','r'), axis tight
        % % surface filling vs. time
        % hplotprop = getprop(o.hplot);
        % ax2 = axes(hplotprop,'Color','none','YAxisLocation','right','YTick',0:0.2:1,'YTickLabel',0:0.2:1,'YLim',[0 1]);
        % color = cbrewer('seq','YlGnBu',numprctile+1);
        % for i=1:numprctile
        %     plot(o.realtime(o.ind),surfacered(o.ind,i),'o','color',color(i+1,:))
        %     line(o.realtime,surfacered(:,i),'color',color(i+1,:),'markerfacecolor','r')
        % end
        
        % % NEW VERSION
        % filling ratio as calculated from variogram
        subplot(o.hplot(1)), hold on, box off
        plot(o.realtime ,fillingratio,'k-','linewidth',2)
        plot(o.realtime(o.ind),fillingratio(o.ind),'ro','markerfacecolor','r'), axis tight, set(gca,'YTick',0:0.2:1,'YTickLabel',0:0.2:1,'YLim',[0 1])
        
        if strcmp('debug',o.plot)
            % surface filling vs. time (percentiles)
            hplotprop = getprop(o.hplot);
            ax2 = axes(hplotprop,'Color','none','YAxisLocation','right','YTick',0:0.2:1,'YTickLabel',0:0.2:1,'YLim',[0 1]);
            color = cbrewer('seq','YlGnBu',numprctile+1);
            for i=1:numprctile
                plot(o.realtime(o.ind),surfacered(o.ind,i),'o','color',color(i+1,:))
                line(o.realtime,surfacered(:,i),'color',color(i+1,:),'markerfacecolor','r')
            end
            % surface filling vs. time (total of intensity / number of pixels in polygon area)
            ax3 = axes(hplotprop,'Color','none','YAxisLocation','right','YTick',0:0.2:1,'YTickLabel',0:0.2:1,'YLim',[0 1],'visible','off'); %#ok<NASGU>
            for i=1:npolygon
                plot(o.realtime,intensitypolygon,'r-')
                plot(o.realtime(o.ind),intensitypolygon(o.ind),'ro') %axis tight
            end
        end
        % cpos = get(ax3,'position');
        % ax3display = axes(hplotprop,'position',[cpos(1)+cpos(3)+.09 cpos(2) 0.0001 cpos(4)],'XTick',[],'XColor',[0.7 0.7 0.7],'Ycolor',[1 0 0],'Color','none','YAxisLocation','right','YLim',[0 1]);

        % print scale
        scalebarwidth = 10*round(framewidth*(o.scalebarwidth/10)); % bar: scalebar size displayed in pixel unit
        scale = 10*round(size(rgb,2)*o.pixellength*(o.scalebarwidth/10)); % scale: scalebar size diplayed in micron unit
        subplot(o.hframes),
        tmp = arrayfun(@(i) sprintf('\\color[rgb]{%0.1g %0.1g %0.1g}\\fontsize{11}%s',o.textcolor(1),o.textcolor(2),o.textcolor(3),sprintf('           + %d s',i)),o.realtime(o.ind),'UniformOutput',false)';
        legblockimage(blocks,[],struct('text',tmp,'position',2),0.05)
        tmp = arrayfun(@(i) '',ones(1,nframestodisplay-1),'UniformOutput',false)';
        legblockimage(blocks,[],struct('text',{sprintf('\\color[rgb]{%0.1g %0.1g %0.1g}\\fontsize{11}%d \\mum',o.textcolor(1),o.textcolor(2),o.textcolor(3),scale);tmp}','line',scalebarwidth,'step',0.5,'position',3,'color',o.linecolor),0.1)
        if ~isempty(o.polygon)
            points = cell2mat(o.polygon);
            for i=1:o.xsubplots
                for j=1:o.ysubplots
                    line(points([1:end,1],1)+10+(20+framewidth)*(i-1),points([1:end 1],2)+10+(20+frameheight)*(j-1),'color',o.linecolor)
                end
            end
        end

        % set xlabel and ylabel for hplot
        %formatax(o.hplot,'fontsize',10,'box','off')
        subplot(o.hplot), xlabel({'time (s)'},'fontsize',12), ylabel({'filling ratio'},'fontsize',12) %, ylabel({'surface filling'},'fontsize',12,'Parent',ax2)
        if strcmp('debug',o.plot)
            ylabel({'percentiles and' 'oil intensity ratio in cell'},'fontsize',12,'Parent',ax2)
        end
    end

% populate outputs
    if nargout>0, hframes(ipolygon) = o.hframes; end
    if nargout>1, hplot(ipolygon) = o.hplot; end
    if nargout>2,
        seriesproperties(ipolygon).seriesfilename = regexprep(o.series,{',|;' '*' '=' ':' '\\' '/' ' ' '_+-+' '-+_+' '_+'},{'_' 'x' '_' '-' '_' '_' '_' '_' '_' '_'});
        seriesproperties(ipolygon).I = I;
        seriesproperties(ipolygon).fillingratio = fillingratio;
        seriesproperties(ipolygon).variogram = v;
        seriesproperties(ipolygon).blocks = blocks;
        seriesproperties(ipolygon).inputs = o;
    end
    
end % next ipolygon