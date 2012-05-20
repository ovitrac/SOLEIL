function out = samplesoleil(db,varargin)
%SAMPLESOLEIL samples images from a database for rapid analizes
%   syntax: sample = moviesoleil(db,options) with options = struct('param1',value1,'param2',value2,...)
%   syntax: sample = moviesoleil(db,'param1',value1,'param2',value2,...)
%   syntax: sample = moviesoleil(db,options,'param3',value3,'param4',value4,...)
%
% INPUTS
%--------
%   db: nexperimentx1 structure array as defined by loadadbsoleil()
%
%   Main properties (parameter: default value)
%          series: mx1 array of experience indices
%        imagefmt: '*.tif'
%           class: uint8
%      resolution: 128 (vales < 16, force wavelets instead of imresize)
%         prctile: [0.1 99.9]
%            step: 1
%
% OUTPUTS
%----------
%   sample = mx1 structure array with fields
%        idx: db index
%       file: nframesx1 cell rray of filenames
%         im: [resolutionxresolutionxnframes of class]
%        ind: frame indices, empty for all (default=[])
%     thresh: [2x2 single] row-wise thresholds
%     frames: frame indoces
%       step: scalar
%         id: file id
%        nfo: file info
%
%   Example 1:
%       db = loaddbsoleil('local',pwd)
%       an = loadansoleil('local',pwd,'anfile','analysis_results.ods','sheetname','scenario1')
%       sample = samplesoleil(db,'series',[an(1:2).series],'step',3);
%       figure(1), clf, hold on
%       plot([sample(1).frames],immean(sample(1).im,sample(1).thresh),'b-')
%       plot([sample(2).frames],immean(sample(2).im,sample(2).thresh),'r-')
%       figure(2), clf
%       cax = [min(sample(2).thresh([35  50],1)) max(sample(2).thresh([35  50],2))];
%       subplot(121), imagesc(imscale(sample(2).im(:,:,35),sample(2).thresh(35,:))), caxis(cax)
%       subplot(122), imagesc(imscale(sample(2).im(:,:,50),sample(2).thresh(50,:))), caxis(cax)
%       figure, monotone(immean(sample(1).im,sample(1).thresh))
%
%   Example 2:
%       % Series = 2, index of the frame used as example = 150
%       %% load database and sample a series
%       db = loaddbsoleil; %('dbfile','/home/olivier/SUN_interpretation/Telemos/description_results.ods');
%       sample = samplesoleil(db,'series',1,'step',1,'resolution',4); % change series number for each script
%       class = classsoleil(sample,6);
%       scale = (size(sample.im,1)-1)/(db(1).imfinfo.Width-1); downscaling = @(roi) [round((roi(1:2)-1)*scale)+1 round(roi(3:4)*scale)];
%       %% plots (update/add/remove rois if needed)
%       rois = {[1 1 1024 1024]
%               [295 500 101 251]
%               [430 625 120 201]
%               [820 1 201 371]  }; nrois = length(rois); colrois = jet(nrois);
%       roisdown = cellfun(@(roi) downscaling(roi),rois,'UniformOutput',false);
%       I = cellfun(@(roi) squeeze(mean(mean(sample.im(roi(2):(roi(2)+roi(4)-1),roi(1):(roi(1)+roi(3)-1),:),1),2)),roisdown,'UniformOutput',false); I = cellfun(@(val) val/val(1),I,'UniformOutput',false);
%       formatfig(figure,'paperposition',[0.33952 0.58872  20.305 28.5],'figname','example'), hs = subplots(1,[.8 1 .5],0,0.04); colormap('bone')
%       subplot(hs(1)), imshow(imread(char(db(1).files(150)))), cellfun(@(roi,col) rectangle('Position',roi,'EdgeColor',col,'FaceColor','none','LineWidth',2),rois,num2cell(colrois,2)), brighten(.8)
%       subplot(hs(2)), plotpub(I,'marker','none','color',col,'linewidth',2);
%       subplot(hs(3)), plot(class)
%       formatax(hs,'fontsize',12); xlabel('frame index','fontsize',12), subplot(hs(2)), ylabel('relative variation','fontsize',12); subplot(hs(3)), ylabel({'class of' 'background'},'fontsize',12), titles(hs,'','fontsize',14,'x',.95,'y',.95)



% See also: immean, loaddbsoleil, loadansoleil, moviesoleil, synchrotronpattern, synchrotronpattern2

% Soleil experiments SUN2011c SUN2011d - 08/01/12 - INRA\Olivier Vitrac - rev. 10/03/12

% revision history
% 08/01/12 release candidate
% 10/01/12 fix thres instead of thresh, fix normalization
% 12/01/12 enable wavelet
% 27/01/12 add example 2, send to linux
%          for ip in 19 20 21 22 23 24 36 27 28 111 112 113 114; do scp -p /home/olivier/codes/confocal/samplesoleil.m olivier@10.75.4.$ip:/home/olivier/codes/confocal/. ; done
%          scp -p /home/olivier/codes/confocal/samplesoleil.m Olivier@10.75.4.102:/cygdrive/d/Data/Olivier/INRA/Codes/confocal/.
%          scp -p /home/olivier/codes/confocal/samplesoleil.m Olivier@10.75.4.115:/cygdrive/c/Data/Olivier/INRA/Codes/confocal/.
% 09/03/12 fix same mean for coefficients of wavedec2() (same problem as in synchrotronpattern2, rev. 23/02/12)
% 10/03/12 add ind

% raw properties

% default properties
o_default = struct(...
     'series',[],...
     'imagefmt','*.tif',...
        'ind',[], ...
      'class','uint8',...
 'resolution',128,...
    'prctile',0.1,...
       'step',1, ...
    'wavelet','db1' ...
);

% argcheck
if nargin<1, error('one argument is at least required'); end
if ~isstruct(db) || isempty(db), error('db must be a valid database'); end
o = argcheck(varargin,o_default,'case');
o.class = lower(o.class);
if ~ismember(o.class,{'uint8','uint16','single','double'}), error('unrecognized class ''%s''',o.class); end
if isempty(o.series), error('a valid series is required'), end
if length(o.resolution)<2, o.resolution = o.resolution([1 1]); end, o.resolution = o.resolution(1:2);
if length(o.prctile)<2, o.prctile = [o.prctile(1) 100-o.prctile(1)]; end, o.prctile = o.prctile(1:2); o.prctile=o.prctile(:)';
imresizeon =  all(o.resolution>15);

% for each series
nseries = length(o.series);
out = repmat(struct('idx',[],'file','','im',[],'thresh',[],'frames',[],'step',[],'id','','nfo',''),nseries,1);
i = 0;
for idx = o.series
    i = i + 1;
    out(i).idx = idx;
    out(i).file = explore(o.imagefmt,db(out(i).idx).fullpath,[],'abbreviate');
    nfiles = length(out(i).file);
    if isempty(o.ind), ind = 1:o.step:nfiles; else ind = o.ind; end
    nind = length(ind);
    out(i).nfo = imfinfo(fullfile(out(i).file(1).path,out(i).file(1).file));
    out(i).id = lastdir(out(i).file(1).path);
    out(i).thresh = zeros(nind,2);
    out(i).frames = zeros(nind,1);
    if imresizeon
        out(i).im = zeros(o.resolution,o.class);
    end
    screen = ''; j = 0;
    for iframe=ind
        j = j + 1; completion = 100 * ((i-1)*nind + j)/(nseries*nind);
        screen=dispb(screen,'SAMPLESOLEIL series %d (%s): loading frame %d/%d [%0.4g%% completed]...',out(i).idx,out(i).id,iframe,nfiles,completion);
        tmp = single(imread(fullfile(out(i).file(iframe).path,out(i).file(iframe).file)));
        out(i).thresh(j,:) = prctile(tmp(:),o.prctile);
        out(i).frames(j) = iframe;
        if imresizeon
            tmp = (tmp-out(i).thresh(j,1))/(out(i).thresh(j,2)-out(i).thresh(j,1));
            switch o.class
                case 'uint8', tmp = uint8(256*tmp);
                case 'uint16', tmp = uint16(65536*tmp);
                case 'single', tmp = single(tmp);
                case 'double', tmp = double(tmp);
            end
            out(i).im(:,:,j) = imresize(tmp,o.resolution);
         else % use wavelet instead without scaling and single precision (instead of uint8)
             tmpmean = mean(tmp(:));
            [C,S] = wavedec2(single(tmp),min(wmaxlev(size(tmp),o.wavelet),o.resolution(1)),o.wavelet);
            tmp = appcoef2(C,S,o.wavelet);
            tmp = tmpmean * tmp / mean(tmp(:)); 
            if j==1, out(i).im = zeros(size(tmp),'single'); end
            out(i).im(:,:,j) = tmp;
        end
    end
    out(i).step = o.step;
end