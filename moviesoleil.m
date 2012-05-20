function out = moviesoleil(db,varargin)
%MOVIESOLEIL main function to load/analyze/visualize SUN2011c results (multipurpose/integrated function to simplify user scripts)
%   syntax: out = moviesoleil(db,options) with options = struct('param1',value1,'param2',value2,...)
%   syntax: out = moviesoleil(db,'param1',value1,'param2',value2,...)
%   syntax: out = moviesoleil(db,options,'param3',value3,'param4',value4,...)
%   syntax: out = moviesoleil(db,...,keyword/flag,...)
%
%*******************
%      INPUTS
%*******************
%
%          db: nexperimentx1 structure output of loadadbsoleil()
%       
%
%        >>MAIN PROPERTIES ('parameter name': value, also noted in examples paremeter=value)
%          ---------------
%                  series: mx1 array of experience indices
%
%                     ind: nx1 array of frames indices (ceonsecutive or not)
%
%                   image: char or string cell array to define images to diplay
%                          Recommended reconstruction value : 'rgb3' (most robust and last method)
%                           List of accepted values:
%                               'im': raw image
%                              'rgb': composite image based on a median threshold
%                                     with cells in green
%                                      with oil colored in red (upper half intensity scale)
%                             'rgb2': composite image based on the difference with a reference frame
%                                      cells in green with an intensity weighted by factorgreen (default = 0)
%                                     differences weighted by factorred (default=0) and factorblue (default=0)
%        recommended---->     'rgb3': as 'rgb2' but use diffsoleil instead of standard image substraction
%                                     This non-parametric method does not need explicit model of fading/bleaching.
%                                     'rgb3_num_iter' can be set to use for anisotropic filtering.
%                                     set 'rgb_rescale' to true to rescale the result of diffsoleil (using 'prctile' values)                                    
%                           'scaled': image with no background (is included as soon as rgb, rgb2 or rgb3 is calculated)
%                             'none': load images only (as im but without display)
%
%                     roi: region of interest to be defined as [xleft yleft width heigth]
% 
%synchrotronpatternengine: char to set the illumination correction method (see the function with similar name for details)
%                          Recommended correction values: 'synchrotronpattern2', 'homomorphicfiltersoleil' (when the previous one fails)
%                          List of accepted values:
%                       '' or 'none': no correction
%               'synchrotronpattern': low frequency components are retrieved by csaps (2D cubic spline)
%                                     NOTE that results are rescaled with percentiles (by default 0.1 and 99.9th) to 255
%  OK ---->    'synchrotronpattern2': low frequency components are retrieved by wavedec2 (2D DWFT)
%                                     NOTE that results are rescaled with percentiles (by default 0.1 and 99.9th) to 255
%                                     TIP: use 'noscaling' along with 'synchrotronpatternparam'
%  OK ---->'homomorphicfiltersoleil': illumination correction is performed in Fourier space (2D FFT) and applied along a Butterworth filter
%                                     NOTE that homomorphicfiltersoleil works a slightly differently from synchrotronpattern, synchrotronpattern2
%                                     Differences: i) no scaling is applied (it is not the intend)
%                                                 ii) when background=[], the previous frame is used as background reference
%                                                     differences between images are inferred as exp(cumsum(log(1+out.scaled),3)) instead as
%                                                     standard differences
%                                                iii) setting explicitely frame indices/classes (or even images) in background enables
%                                                     to return to the standard behavior using diffsoleil (rgb3) or standard substraction (rb2)
%                                                 iv) 'rgb_rescale' is set to true (due to the absence of scaling)
% 
%             'background': index of frames to be used as to determine and correct illumination (averaging over index) with synchrotronpatternengine
%                          When used with 'synchrotronpattern2' and 'homomorphicfiltersoleil', a by-class-approach is also available.
%                          1xnframes vector (output of classsoleil) giving the class index of each frame (the first frame of each class is used)
%                          An automatic background determination is enabled with the keyword 'backgroundauto'
%                          An empty value (default) mean all frames, except when synchrotronpatternengine is set to 'homomorphicfiltersoleil'
%                       
% 
%             'reference': index of frames to be used as reference for image substraction (averaging is applied if needed)
%                          when used with 'rgb2' an experimental cumulant is available using a by-class-approach:
%                          1xnframes vector (output of classsoleil) giving the class index of each frame
%
%      RGB intensities for each channel
%                      'factorgreen': scalar value (default = 1)
%                        'factorred': scalar value (default = 4)
%                       'factorblue': scalar value (default = 0)
%
%
%         --------------------------
%       >>PROPERTIES ACTING AS FLAGS (parameter: value) that control the behavior of MOVIESOLEIL
%         --------------------------
%
%                          'videoon': flag to grab a video (default = false)
%                          'pauseon': flag to pause display (default = false)
%                            'clear': internal flag to accumulate results (default = false), must not set by user
%
%
%                >>KEYWORDS: as FLAGS but they are assumed to be false (not set) when non-used
%                  -------- TIP: Keywords cannot be stored in the structure options
%
%                  'noplot' removes plots (faster)
%               'lowmemory' preserves memory (empty scale an im fields at each iteration)
%          'backgroundauto' enables an automatic extraction of background based on the most populated class of all sampled frames (when
%                           the most populated class is not included in ind, the second is used etc). Finally only the frames included in
%                           ind and of the proper class are considered to define an average image used as background. As a result, the same
%                           correction is applied to all images.
%         'backgroundauto2' as 'backgroundauto' but do a specific correction for each class of frames defined in ind.
%                           TIP 'backgroundauto2' is faster than 'backgroundauto' because classification is performed
%                           only on considered frames
%
%         --------------------------
%       >>ADVANCED PROPERTIES/VALUES (for non-standard applications)
%         --------------------------
%
%       --> PARAMETERS AFFECTING HOW FRAMES ARE LOADED IN MEMORY
%                  imagefmt: file pattern to find frames db.fullpath with EXPLORE (default = '*.tif')
% 
%       --> PARAMETERS AFFECTING RGB RECONSTRUCTION/ASSEMBLING (only 'rgb3' is recommended for production)
%                       rgb_rescale: when set to true (default = false), the difference result is rescaled between 0 and 255 (using 'perctile' values)
%                                     this flag have effect on 'rgb2' 'rgb3' ('rgb' is not altered by this parameter)
%       When  'rgb' is used as image
%                rgbthreshpercentile: percentiles to scale intensities (default = [1 4])
%       When  'rgb2' is used as image
%                 smoothrgb: set as true to use gaussian filter on rgb2 images
%                            use inputs: 'hsize' (default = 100)
%                                        'sigma'    (default = 30)
%       When  'rgb3' is used as image (diffsoleil() properties mainly)
%                      rgb3_num_iter: 2x1 array to apply anisotropic filtering to source and final image before substraction with diffsoleil            
%                            delta_t: anisodiff2D parameter (default = 1/7)
%                              kappa: anisodiff2D parameter (default = 10)
%                             option: anisodiff2D parameter (default = 2)
%                              depth: scattergramsoleil parameter (default = 256)
%                            prtcile: scattergramsoleil parameter (default = [.1 99.9])
%
%
%       --> ILLUMINATION CORRECTION PARAMETERS (only 'synchrotronpattern2' and 'homomorphicfilersoleil' are recommended for production)
%       The property 'synchrotronpatternparam' must be followed by a structure with fields that depend on the method set in synchrotronpatternengine
%       General syntax: moviesoleil(...,'synchrotronpatternparam',struct('parameter1',value1,'parameter2',value2,...),...)
%       When 'synchrotronpattern' is used as synchrotronpatternengine (depreciated, old method)
%                                  p: csaps regularizing parameter (default = 1e-5)
%                                  f: f or [fy fx] median filter (default = 5)
%                                  q: q, [q 100-q], [qmin qmax] threshold to match 0-255 (default = 0.1)
%       When 'synchrotronpattern2' is used as synchrotronpatternengine (note that only some properties are public)
%                            details: 'v', 'h' or 'd' for vertical, horizontal and diagonal details respectively (see detcoef2 for details)
%                       detailslevel: wavelet level decomposition (see wmaxlev for details)
%                            interp: interpolation method for coefficients (default = 'lanczos3'), see imresize
%                                    available methods: 'nearest', 'box', 'triangle', 'cubic', 'lanczos2', 'lanczos3'
%                             level: wavelet level decomposition (see wmaxlev for details)
%                        medianfilt: median filter width (default = []), set to 1, [] or NaN to remove median filtering
%                         prctile: threshold percentiles for uint8 conversion (default = [0.1 99.9])
%                         wavelet: wavelet type (default = 'db1')
%                                  TIP: wavemngr('read',1) list available wavelets and their shorthands to be used
%                       noscaling: if true force results not to be scaled
%       When 'homomorphicfilersoleil' is used as synchrotronpatternengine
%                           order: filter order (default = 2)
%                          cutoff: Lower cut off frequency (default = 10)
%                        passband: [0.0999 1.01]
%
%       --> FADING CORRECTION PARAMETERS (not required when 'rgb3' is used as image)
%                bleachcorrection: anonymous function to correct bleaching effects e.g. @(I0,t) I-I0*(slopevsI0(1) + slopevsI0(2)*I0)*t
%                            where slopevsI0 is a parameter in fading.m and t is the number of frame-1 i.e. beginning at 0 (default = [])
%
%                 >>VIDEO PARAMETERS (alternatively use grabmoviesoleil() to generate movies based on RGB outputs)
%                           videofps: frame rate per second (default = 10)
%                       videoquality: best = 100 (default = 40)
%                   videocompression: codec (default = 'none' on linux and 'divx' on windows)
%                                     install codec: http://www.xvid.org/Downloads.15.0.html
%
%
%*******************
%      OUTPUTS
%*******************
% Note that some fields (rgb, rgb2, rgb3) below can be empty according to the values of 'image'.
% scaled and im are empty when the keyword 'lowmemory' is used.
%
%   out = mx1 structure array with fields
%         id: automatic identifier (e.g. '4_darkbatter_permeation_stack')
%     series: db index
%        ind: frame indices
%     length: number of frames
%       file: [1xnframes struct]
%        nfo: [1x1 struct]
%         im: [HeightxWidthxnframes uint16]
%        min: min intensity value
%        max: max intensity value
%     thresh: [2x2 single] row-wise thresholds matching rgbthreshpercentile
%     median: median intensity value
%        rgb: [HeightxWidthx3xnframes uint8] rgb image --> (to be used with grabmoviesoleil, montage, imshow)
%       rgb2: [HeightxWidthx3xnframes uint8] rgb image --> (to be used with grabmoviesoleil, montage, imshow)
%       rgb3: [HeightxWidthx3xnframes uint8] rgb image --> (to be used with grabmoviesoleil, montage, imshow)
%     scaled: [HeightxWidthxnframes single] scaled image (rescaled to 255 except if synchrotronpatternengine is set to 'homomorphicfiltersoleil')
%
%
%   see also: loaddbsoleil, diffsoleil, grabmoviesoleil, classsoleil, scattergramsoleil, synchrotronpattern, synchrotronpattern2, homomorphicfilersoleil, montage, imshow, hist, histfit, histfitlog
%
%******************************
%     OBSOLETE EXAMPLE
%******************************
% PREHISTORIC EXAMPLE (series 38 - frames 20 to 31)
% ---------------------------------------------------
%     db = loaddbsoleil;
%     figure
%     out = moviesoleil(db,'series',38,'ind',20:31,'pauseon',true); % press a key to see all images
%     figure, montage(out.rgb) % the whole series with rgb colors
%     figure, montage(permute(out.scaled,[1 2 4 3])), caxis([0 255])
%
%
%**********************************************
%     PRODUCTION EXAMPLES: LOCAL OIL UPTAKE
%  Best resulats are obtained with:
%       synchrotronpatternengine = 'synchrotronpattern2' (intensity is normalized with the DWFT of a reference image)
%       image = 'rgb3' (the red channel is reconstructed by image substraction using diffsoleil)
%       For convenience, the keyword 'backgroundauto' can be used to set automatically the best background.
%       When 'backgroundauto' is not used, it is recommended to set explicetly background and reference frames.
%       Advanced users can combine samplesoleil et classsoleil to optimize their background and reference (not detailed here)
%  TIP: do not use a large value for factorred when 'rgb_rescale'=true
%**********************************************
%     db = loaddbsoleil;
% %fast test (without backgroundauto) - note that frames can be discontinuous
%     defaultanalysis = struct('factorgreen',1,'factorred',4,'factorblue',0,'synchrotronpatternengine','synchrotronpattern2','medfilt',5,'videoon',false);
%     out = moviesoleil(db,'series',227,'image','rgb3','roi',[],'ind',[1 25 50 100 200 300],'reference',25,'rgb3_num_iter',8,defaultanalysis,'noplot');
%     figure, montage(out.rgb3)
% %do analysis with ROI and 'backgroundauto' (without generating any movie)
%     out = moviesoleil(db,'backgroundauto','series',20,'image','rgb3','roi',[430 430 256 256],'ind',1:50:300,'rgb3_num_iter',8,defaultanalysis,'noplot');
%     figure, montage(out.rgb3)
%
%**********************************************
%     PRODUCTION EXAMPLES: OIL INVASION
%  TIP: Usually oil invasion leads to a non-uniform background correction between 'dry' and 'wet' states.
%       As a result, synchrotronpatternengine='synchrotronpattern2' tends to give poor results.
%       Better results are obtained with synchrotronpatternengine='homomorphicfilersoleil';
%       Two methods are recommended to calculate the red channel (oil):
%               1) factorizing in Fourier Space A/B*B/C*C/D... as log(A)-log(B)+log(B)-log(C)+log(C)-log(D)
%                  homomorphic filtering ensures that only low-frequency details are removed (main lighting)
%                  high-frequency details (phase) are kept and acculated in a conservative way (low risk of divergence)
%                 ==> this mode is activating by setting backgroung to []
%                  PRO: the gives high contrast and is fast (no image substraction)
%                  CON: The method works best when it is applied to consecutive frames
%                       The reference channel (green) is based on the homomorphic filtering of the first image
%                  NOTE that rgb_rescale is aplplied automatically.
%               2) calculating (B/A-A/A)+(C/A-B/A)-(D/A-C/A)+... as (exp(log(B+1))-exp(log(A+1)))+(exp(log(C+1))-exp(log(B+1)))+...
%                  this method compares each frame to a same frame or set of frames and use diffsoleil to perform image substractions
%                  PRO: The method computes absolute differences (and not local variations), more sensitive
%                       The method can be applied to non-consecutive frames (e.g. [58:68 100 150 200 250])
%                       The background (green channel) is updated.
%                  CON: The method uses heavy resources.
%                  ATTENTION: be-sure to apply rgb_rescale=true
%**********************************************
% SERIES 120 (BATTER): ind = 128:250; synchrotronpatternengine='homomorphicfilersoleil'; background=[];
%     db = loaddbsoleil;
%     defaultanalysis = struct('factorgreen',1,'factorred',4,'factorblue',0,'synchrotronpatternengine','synchrotronpattern2','medfilt',5,'videoon',false);
%  %method 1 on consecutive frames - note that background = [] and image = 'rgb3' (even if diffsoleil is not used)
%     out = moviesoleil(db,'series',120,'ind',58:100,'background',[],'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil',defaultanalysis,'noplot');
%     figure, montage(out.rgb3(:,:,:,1:30));
%  %method 2 on the same frames (for comparison) - note that rgb_rescale=true
%     out = moviesoleil(db,'series',120,'ind',58:100,'background',58,'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil',defaultanalysis,'noplot','rgb_rescale',true)
%     figure, montage(out.rgb3(:,:,:,1:30));
%  %method 2 on non-consecutive frames (very fast)
%     out = moviesoleil(db,'series',120,'ind',[58:68 100 150 200 250],'background',58,'image','rgb3','synchrotronpatternengine','homomorphicfilersoleil',defaultanalysis,'noplot','rgb_rescale',true)
%     figure, montage(out.rgb3)
%
%
%**********************************************
%     HOW TO CREATE A MOVIE
%**********************************************
% rgb3 analysis - series 54
%     db = loaddbsoleil;
%     defaultanalysis = struct('factorgreen',1,'factorred',4,'factorblue',0,'synchrotronpatternengine','synchrotronpattern2','medfilt',5,'videoon',false);
%     out = moviesoleil(db,'series',54,'image','rgb3','roi',[430 511 340 340],'ind',13:1:40,'background',13:15,'rgb3_num_iter',2,defaultanalysis,'noplot');
%     figure, montage(out.rgb3(:,:,:,2:4:end))
%     grabmoviesoleil(out.rgb3,'id','RGB3_series54','overwrite')

% Soleil experiments SUN2011c - 17/11/11 - INRA\Olivier Vitrac - rev. 19/03/12

% revision history
% 17/11/11 release candidate
% 23/11/11 fix image to 'none' to prevent graphical outputs
% 04/01/12 rename indref, reference as background (for compatibility with synchrotronpattern2)
% 04/01/12 add rgb2, indback replaces indref, indref for reference, add keywords 'noplot', 'lowmemory'
% 04/01/12 workaround for XVID compression error, add Series%03d to movie names
% 05/01/12 fix erravi
% 09/01/12 add smoothrgb as option of rgb2: use gaussian filter (add hsize and sigma), help update
% 12/01/12 add factorgreen, help update
% 12/01/12 background based on classes
% 13/01/12 fix old behavior without class
% 16/01/12 manage several references as cumulants
% 17/01/12 add refpast, refprctile
% 10/02/12 add bleachcorrection
% 17/02/12 add rgb3 based on diffsoleil
% 18/02/12 add roi, rgb3_num_iter
% 19/02/12 updated help, clear sample
% 21/02/12 add translate o.roi to the top left if required
% 23/02/12 retrieve correctly the value of synchrotronpatternparam as a structure
% 24/02/12 no fading correction when ntrend<2
% 26/02/12 help major revision, add homomorphicfilersoleil
% 27/02/12 new examples, add rgb3_rescale
% 28/02/12 2nd release candidate with tested examples
% 10/03/12 add backgroundauto2
% 19/03/12  fix "convert to rgb2 or rgb3" 'line 536) when o.image is empty (when image is set to 'none')

% raw properties
synchrotronpatternparam = struct('p',1e-4,'f',5,'q',0.1);

% default properties
o_default = struct(...
    'clear',false,...
    'image',{{'im' 'rgb' 'rgb2' 'rgb3' 'scaled'}},...
    'imagefmt','*.tif',...
    'ind',[], ...
    'background',[],...
    'bleachcorrection',[],...
    'factorgreen',1,...
    'factorred',4,...
    'factorblue',0,...
    'hsize',100,...
    'outputfolder','output',...
    'pauseon',false,...
    'reference',[],...
    'refpast',10,...
    'refprctile',50,...
    'resolution',4,...
    'rgbthreshpercentile',[1 4],...
    'series',[],...
    'sigma',30,...
    'synchrotronpatternengine','synchrotronpattern',... synchrotronpattern, synchrotronpattern2, '' or 'none'
    'synchrotronpatternparam',synchrotronpatternparam,...
    'videoon',false,...          true to grab a move
    'videofps',10,...            frame/s
    'videoquality',40,...        image quality (not used with codec divx)
    'videocompression','none', ... default
    'delta_t',1/7,...
    'kappa',10,...
    'option',2, ...
    'depth',256,...
    'prctile',[.1 99.9],...
    'roi',[],...
    'rgb3_num_iter',0,...     'num_iter',2,...
    'rgb_rescale',false ...
);
kw = {'noplot' 'lowmemory' 'smoothrgb','backgroundauto','backgroundauto2'};
result = struct('id',[],'idx',[],'ind',[],'length',[],'extlength',[],'file',[],'nfo',[],'im',[],'min',[],'max',[],'scaledthres',[],'thresh',[],'median',[],'rgb',[],'scaled',[],'rgb2',[],'rgb3',[]);
[red,green,blue] = deal(1,2,3);

% argcheck
o = argcheck(varargin,o_default,kw,'case');
if ~isstruct(o.synchrotronpatternparam)
    tmp = argcheck(varargin,o_default,kw,'case','nostructexpand');
    if isstruct(tmp.synchrotronpatternparam), o.synchrotronpatternparam = tmp.synchrotronpatternparam; end
end
    
if strcmpi(o.synchrotronpatternengine,'synchrotronpattern')
    o.synchrotronpatternparam = argcheck(o.synchrotronpatternparam,synchrotronpatternparam);
elseif strcmpi(o.synchrotronpatternengine,'none'), o.synchrotronpatternengine = '';
end
% not needed: ind0 = intersect(ind,ind0); %indref is an index of ind

% customization
switch localname
    case 'WSLP-OLIVIER2', override.videocompression = 'xvid'; % Codec: http://www.xvid.org/, Utility: http://virtualdub.sourceforge.net/
    case 'LP-MOL5',       override.videocompression = 'xvid';
    otherwise, override = [];
end
% argcheck
o = argcheck(override,o,'','case');
if ischar(o.image) && strcmpi(o.image,'none'), o.image = ''; end
if ~isempty(o.image)
    if ~iscell(o.image), o.image = {o.image}; end
    if ~iscellstr(o.image), error('image must be a char cell array or a string'), end
end
nplots = length(o.image);
if o.noplot, nplots = 0; end

% autobackground
sampledone = false;
if (o.backgroundauto) || (o.backgroundauto2 && isempty(o.ind))
    o.backgroundauto2 = false;
    dispf(['\n' repmat('-',1,50)])
    dispf('\nbackground will be automatically set with CLASSSOLEIL...\nTIP: use ''backgroundauto2'' to restrict the analysis to supplied frames.') 
    if length(o.series)>1, error('backgroundauto is available only for unique series'), end
    sample = samplesoleil(db,'series',o.series,'step',1,'resolution',o.resolution); sampledone = true;
    nclass = max(3,ceil(sqrt(size(sample.im,3))/2.5));
    dispf('\nbackground set with CLASSSOLEIL (using %d classes)...',nclass) 
    currentclass = classsoleil(sample,nclass);
    nframesinclass = arrayfun(@(c) length(find(currentclass==c)),1:nclass); % number of frames in the two first classes (possibly extend 1:2 to 1:3)
    [nframesinclass,sortedclasses] = sort(nframesinclass,'descend');
    foundclass = false; iclass = 0;
    while ~foundclass && iclass<nclass
        iclass = iclass + 1;
        if ~isempty(intersect(o.ind,find(currentclass==sortedclasses(iclass)))), foundclass = true; end
    end
    dispf('\nbackground set with CLASSSOLEIL: class %d/%d is found optimal (ranked #%d) including %d frames between %d and %d...',...
        sortedclasses(iclass),nclass,iclass,nframesinclass(iclass),find(currentclass==sortedclasses(iclass),1,'first'),find(currentclass==sortedclasses(iclass),1,'last')) 
    o.background = o.ind(currentclass(o.ind)==sortedclasses(iclass));
    dispf('\nbackground set to: [ %s]',sprintf('%d ',sort(o.background)))
    dispf(['\n' repmat('-',1,50)])
elseif o.backgroundauto2 && ~isempty(o.ind)
    dispf(['\n' repmat('-',1,50)])
    dispf('\nbackground will be automatically set with CLASSSOLEIL (using supplied frames only)...') 
    if length(o.series)>1, error('backgroundauto is available only for unique series'), end
    sample = samplesoleil(db,'series',o.series,'ind',o.ind,'step',1,'resolution',o.resolution); sampledone = true;
    % figure, figure, montage(permute(uint8(255*(sample.im-min(sample.thresh(:,1)))/(max(sample.thresh(:,2))-min(sample.thresh(:,1)))),[1 2 4 3]))
    nclass = max(3,ceil(sqrt(length(o.ind))/2.5));
    o.background = classsoleil(sample,nclass);
end

% keep memory
if sampledone, clear sample, end

% argcheck - continued
useclasson = ~isempty(findduplicates(o.background)); %true if background is the output of classsoleil
if useclasson, classidx = o.background; o.background = []; else classidx = []; end
if ~isempty(o.roi) && length(o.roi)~=4, error('roi must a 1x4 vector [xleft ybottom width height]'), end
if length(o.rgb3_num_iter)<2, o.rgb3_num_iter = o.rgb3_num_iter([1 1]); end

% result
if isempty(o.series), error('no series is defined'); end
if any(o.series>length(db)), error('invalid series numbers, please check'), end
if o.clear, out = result; incr = 0; else out = repmat(result,length(o.series),1); incr = 1;  end

% output folder
if ~exist(o.outputfolder,'dir'), mkdir(o.outputfolder), end

% for each series
i = 0;
for idx = o.series
    i = i + incr;
    
    % series
    out(i).idx = idx;
    
    % common
    out(i).file = explore(o.imagefmt,db(out(i).idx).fullpath,[],'abbreviate');
    nfiles = length(out(i).file);
    if isempty(o.ind),        out(i).ind = 1:nfiles;
    elseif iscell(o(i).ind),  out(i).ind = o.ind{i};
    else                      out(i).ind = o.ind;
    end    
    if isempty(o.background);
        if useclasson
            if o.backgroundauto2, background = o.ind;
            else background = 1:nfiles; end
        else
            background = out(i).ind;
        end
    elseif iscell(o.background), background = o.background{i};
    else                         background = o.background;
    end
    if isempty(o.reference);     reference = out(i).ind(1);
    elseif iscell(o.background), reference = o.reference{i};
    else                         reference = o.reference;
    end
    ind = union(union(out(i).ind,background),reference);
    [~,ind0]    = intersect(ind,out(i).ind); %ind0 is now an index of ind
    [~,indback] = intersect(ind,background); %indback is now an index of ind
    [~,indref]  = intersect(ind,reference); %indref is now an index of ind
    out(i).file = out(i).file(ind);
    out(i).length   = length(out(i).ind);
    out(i).extlength = length(ind);
    out(i).nfo = imfinfo(fullfile(out(i).file(1).path,out(i).file(1).file));
    % Use Roi instead
    if ~isempty(o.roi)
        imsize = [out(i).nfo.Width out(i).nfo.Height];
        outscreen = find(o.roi([1 2])+o.roi([3 4])-1>imsize);
        o.roi(outscreen) = o.roi(outscreen)-(o.roi(outscreen)+o.roi(2+outscreen) - 1 -imsize(outscreen));
        if any(o.roi(1:2)<1), error('the roi [%d %d %d %d] exceeds the image size. Please, decrease the ROI size',o.roi(1),o.roi(2),o.roi(3),o.roi(4)), end
        out(i).nfo.Width = o.roi(3); out(i).nfo.Height = o.roi(4);
    end
    out(i).id = lastdir(out(i).file(1).path);
    out(i).im = zeros(out(i).nfo.Height,out(i).nfo.Width,out(i).extlength,sprintf('uint%d',out(i).nfo.BitDepth));
    outputmovie =  fullfile(o.outputfolder,sprintf('Series%03d_%s_%03d-%03d.avi',idx,out(i).id,min(out(i).ind),max(out(i).ind)));
    
    % check whether the output already exist
    if exist(outputmovie,'file')
        warning('The output of ''%s'' already exists',out(i).id) %#ok<WNTAG>
        fileinfo(outputmovie)
        dispf('Remove the output file above first.\n\tuse: delete(''%s'')',outputmovie)
        o.videoon = false;
    end
        
    % common (suite)
    screen = '';
    for iframe=1:out(i).extlength
        screen=dispb(screen,'MOVIESOLEIL series %d (%s): loading frame %d/%d...',out(i).idx,out(i).id,iframe,out(i).extlength);
        tmp = imread(fullfile(out(i).file(iframe).path,out(i).file(iframe).file));
        if ~isempty(o.roi), tmp = tmp(o.roi(2):(o.roi(2)+o.roi(4)-1),o.roi(1):(o.roi(1)+o.roi(3)-1)); end
        out(i).im(:,:,iframe) = tmp;
        if ~isempty(o.bleachcorrection), out(i).im(:,:,iframe) = out(i).im(:,:,iframe) + o.bleachcorrection(iframe); end
        % to avoid the problem of "bad" first frame (when oils is added at t>t0), out(i).im(:,:,o.reference) can replace out(i).im(:,:,1) ?
    end

    % convert to rgb
    if ismember('rgb',o.image)
        screen = dispb(screen,'MOVIESOLEIL series %d (%s): RGB conversion based on distribution of intensities...',out(i).idx,out(i).id);
        out(i).min = min(out(i).im(:));
        out(i).max = max(out(i).im(:));
        tmp = single(out(i).im(:));
        out(i).thresh = reshape(prctile(tmp,[o.rgbthreshpercentile 100-o.rgbthreshpercentile]),length(o.rgbthreshpercentile),2);
        out(i).median = median(tmp( (tmp>=min(out(i).thresh(:))) & (tmp<=max(out(i).thresh(:))) ));
        out(i).rgb = zeros(out(i).nfo.Height,out(i).nfo.Width,3,out(i).length,'uint8');
        % green component
        greencomponent = tmp; greencomponentmax = out(i).median; %prctile(green((green>out(i).thresh(1)) & (green<out(i).median)),99);
        greencomponent(greencomponent>out(i).median)=greencomponentmax;
        greencomponent = 128*(greencomponent-out(i).thresh(2,1))./(greencomponentmax-out(i).thresh(2,1)); greencomponent(greencomponent>255)=255; greencomponent(greencomponent<0)=0;
        out(i).rgb(:,:,2,:) = permute(reshape(uint8(greencomponent),out(i).nfo.Height,out(i).nfo.Width,out(i).length),[1 2 4 3]);
        % red component
        redcomponent = tmp; redcomponent(redcomponent<out(i).median)=out(i).thresh(2,1); redcomponentmin = prctile(redcomponent(redcomponent>out(i).thresh(2,1)),1);
        redcomponent = 220*(redcomponent-redcomponentmin)./(out(i).thresh(1,2)-redcomponentmin); redcomponent(redcomponent>255)=255; redcomponent(redcomponent<0)=0;
        out(i).rgb(:,:,1,:) = permute(reshape(uint8(redcomponent),out(i).nfo.Height,out(i).nfo.Width,out(i).length),[1 2 4 3]);
    end
    
   
    % scaled intensity with background removed
    if ~isempty(o.synchrotronpatternengine) && ~strcmpi(o.synchrotronpatternengine,'none')
        dispb(screen,'MOVIESOLEIL series %d (%s): background removal...',out(i).idx,out(i).id);
        if strcmpi(o.synchrotronpatternengine,'synchrotronpattern')
            [~,out(i).scaled,out(i).scaledthres] = synchrotronpattern(out(i).im,...
                o.synchrotronpatternparam.p,synchrotronpatternparam.f,synchrotronpatternparam.q);
        elseif strcmpi(o.synchrotronpatternengine,'synchrotronpattern2')
            if useclasson
                [~,out(i).scaled,out(i).scaledthres] = synchrotronpattern2(out(i).im,o.synchrotronpatternparam,'background',indback,'class',classidx(indback)); % applied to all images
            else
                [~,out(i).scaled,out(i).scaledthres] = synchrotronpattern2(out(i).im,o.synchrotronpatternparam,'background',indback);
            end
        elseif strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil')
            if useclasson
                [out(i).scaled,out(i).scaledthres] = homomorphicfilersoleil(out(i).im,classidx,o.synchrotronpatternparam);
            else
                if isempty(o.background)
                    [out(i).scaled,out(i).scaledthres] = homomorphicfilersoleil(out(i).im,[],o.synchrotronpatternparam);
                else
                    [out(i).scaled,out(i).scaledthres] = homomorphicfilersoleil(out(i).im,indback,o.synchrotronpatternparam);
                end
            end
        else error('improper value ''%s'' for property ''synchrotronpatternengine'' [''synchrotronpattern'',''synchrotronpattern2'',''homomorphicfilersoleil'']',o.synchrotronpatternengine)
        end
        ntrend = size(out(i).scaled,3);
        if ntrend>1 %&& ~strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil')
            trend = squeeze(mean(mean(out(i).scaled,1),2));
            ratetrend = polyfit((1:ntrend)',trend,1); screen='';
            if ratetrend(1)<0
                for itrend=1:ntrend
                    screen = dispb(screen,'WARNING::\t fading has been detected, the linear trend is removed for frame %d/%d',itrend,ntrend);
                    out(i).scaled(:,:,itrend) =  ratetrend(2) + out(i).scaled(:,:,itrend)-polyval(ratetrend,itrend);
                end
            end
        end
    end
    
   % convert to rgb2 or rgb3 (method using diff soleil or exp(cumsum(log(1+out.scaled),3)) if homomorphicfilersoleil with no background
   if ~isempty(o.image) && (ismember('rgb2',o.image) || ismember('rgb3',o.image))
       if ismember('rgb3',o.image), rgboutput = 'rgb3'; diffmethod='DIFFSOLEIL'; else rgboutput = 'rgb2'; diffmethod='substraction'; end
       if isempty(out(i).scaled), error('please set a synchrotronpatternengine first'); end
       if (length(indref)==1) %|| ~useclasson % standard behavior
           % initialization
           if strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil') && isempty(o.background) % product factorization
               outref = zeros(size(out.scaled(:,:,1)));
           else
               outref = mean(out(i).scaled(:,:,indref),3); % standard behavior
           end
           if strcmp(rgboutput,'rgb3') && o.rgb3_num_iter(1)>0, outref = anisodiff2D(outref,o.rgb3_num_iter(1),[],[],[],false); end
           out(i).(rgboutput) = zeros(out(i).nfo.Height,out(i).nfo.Width,3,out(i).length,'uint8');
           % dscaling if requis
           dscalingflag =  (o.rgb_rescale || (strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil') && isempty(o.background)));
           if dscalingflag, jset = {[1 out(i).length] 1:out(i).length}; dscalingcalibration=true; dscale = [0 0];
           else jset ={1:out(i).length}; dscalingcalibration=false; dscalemin = 0; dscale = 1;
           end
           for jlist = jset
               iscale = 0;
               for j=jlist{1} %1:out(i).length
                   if dscalingcalibration, screen = dispb(screen,'MOVIESOLEIL series %d (%s): RGB conversion based on %s (calibration %d/%d)...',out(i).idx,out(i).id,diffmethod,j,out(i).length); iscale = iscale+1;
                   else screen = dispb(screen,'MOVIESOLEIL series %d (%s): RGB conversion based on %s (%0.3g %%)...',out(i).idx,out(i).id,diffmethod,100*j/out(i).length);
                   end
                   if strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil') && isempty(o.background) % product factorization
                       if dscalingcalibration
                           d = exp(sum(log(1+out(i).scaled(:,:,1:ind0(j))),3)); % as below but with one step integration
                       else
                           screen = dispb(screen,'MOVIESOLEIL series %d (%s): RGB conversion based on homomorphicfilersoleil CUMULANTS (%0.3g %%)...',out(i).idx,out(i).id,100*j/out(i).length);
                           outref = outref + log(1+out(i).scaled(:,:,ind0(j))); % diffsoleil is not used
                           d = exp(outref);
                       end
                   else % standard behavior
                       if strcmp(rgboutput,'rgb2')
                           d= out.scaled(:,:,ind0(j)) - outref;
                       else % rgb3
                           if o.rgb3_num_iter(2)==0
                               d = diffsoleil(outref,out.scaled(:,:,ind0(j)),o);
                           else
                               d = diffsoleil(outref,anisodiff2D(out.scaled(:,:,ind0(j)),o.rgb3_num_iter(2),[],[],[],false),o);
                           end
                       end
                   end % homomorphicfilersoleil
                   if dscalingcalibration
                       dscale(iscale) = prctile(d(:),o.prctile(iscale));
                   else
                       d = dscale*(d-dscalemin);
                       if strcmpi(o.synchrotronpatternengine,'homomorphicfilersoleil') && isempty(o.background) % product factorization
                               if j==1 % a background is reconstructed
                                   back = homomorphicfilersoleil(out(i).im(:,:,ind0(j)),out(i).im(:,:,ind0(j)),o.synchrotronpatternparam,'dwft',true);
                                   dback = prctile(back(:),o.prctile([1 end]));
                                   out(i).(rgboutput)(:,:,green,1) = 255*o.factorgreen * (back-dback(1))/(dback(2)-dback(1));
                               else
                                   out(i).(rgboutput)(:,:,green,j) = out(i).(rgboutput)(:,:,green,1);
                               end
                       else % not homomorphicfilersoleil
                           if o.rgb_rescale
                               back = out(i).scaled(:,:,ind0(j));
                               dback = prctile(back(:),o.prctile([1 end]));
                               out(i).(rgboutput)(:,:,green,j) = 255*o.factorgreen * (back-dback(1))/(dback(2)-dback(1));
                           else
                               out(i).(rgboutput)(:,:,green,j) = uint8(o.factorgreen * out(i).scaled(:,:,ind0(j)));
                           end
                       end % homomorphicfilersoleil
                       if o.smoothrgb
                           %figure, imagesc( sign(d).*abs(d).^(1.3)./abs(dgauss))
                           dgauss = abs(imfilter(d,fspecial('gaussian',o.hsize,o.sigma)));
                           dtmp = sign(d).*sqrt(dgauss.*abs(d)); %dgauss.*d/mean(mean(d,1),2);
                           out(i).(rgboutput)(:,:,red,j)    = uint8(o.factorred*dtmp);
                           out(i).(rgboutput)(:,:,blue,j)   = uint8(o.factorblue*dtmp);
                       else
                           out(i).(rgboutput)(:,:,red,j)   = uint8(o.factorred*d);
                           out(i).(rgboutput)(:,:,blue,j)  = uint8(o.factorblue*d);
                       end % smooth rgb
                   end % dscaling
               end % next j
               if dscalingcalibration
                   dscalemin = dscale(1); dscale = 255/(dscale(2)-dscale(1)); dscalingcalibration = false;
               else  dscalemin = 0; dscale = 1;
               end
           end % next jlist
       else % several references % behavior added 16/01/11 (could be merge with the default behavior once validated)
           screen='';
           indref = min(max(indref,min(ind0)),max(ind0));
           indref = unique([ind0(1);indref(:);ind0(end)])'; nref = length(indref);
           out(i).rgb2 = zeros(out(i).nfo.Height,out(i).nfo.Width,3,out(i).length,'uint8');
           cumulant = zeros(out(i).nfo.Height,out(i).nfo.Width); lastiref=1;
           outref = mean(out(i).scaled(:,:,max(1,indref(1)-o.refpast):indref(1)),3);
           for j=1:out(i).length
               iref = find(indref<=ind0(j),1,'last'); % reference
               if isempty(iref)
                   d = zeros(out(i).nfo.Height,out(i).nfo.Width);
               else
                   if iref>lastiref
                       okd = (d>0) | (d>prctile(d(:),o.refprctile));
                       cumulant(okd)=cumulant(okd)+abs(d(okd));
                       lastiref=iref;
                       outref = mean(out(i).scaled(:,:,max(1,indref(iref)-o.refpast):indref(iref)),3);
                   end
                   screen = dispb(screen,'MOVIESOLEIL series %d (%s): RGB conversion based on difference with reference %d/%d (%0.3g %%)...',out(i).idx,out(i).id,iref,nref,100*j/out(i).length);
                   d=out.scaled(:,:,ind0(j))-outref;
               end
               %d(d<0)=0; %too strict
               % add cumulant and denoise using wavelets
               [thr,sorh,keepapp] = ddencmp('den','wv',d+cumulant);
               dc = wdencmp('gbl',d+cumulant,'sym4',2,thr,sorh,keepapp);
               %[C,L] = wavedec2(d+cumulant,2,'sym8');
               %THR = wthrmngr('dw2dcompLVL','scarcehi',C,L,2);
               %dc = wdencmp('lvd',d+cumulant,'sym8',2,THR,'s');
               %
               out(i).rgb2(:,:,green,j) = uint8(o.factorgreen * out(i).scaled(:,:,ind0(j)));
               out(i).rgb2(:,:,red,j)   = uint8(o.factorred*(dc));
               out(i).rgb2(:,:,blue,j)  = uint8(o.factorblue*(dc));
           end % next j
       end % if
   end
    
    % grab video frames
    if nplots % something to plot
        fig = gcf;
        formatfig(fig,'Units','pixels','position',[300 200 nplots*out(i).nfo.Width out(i).nfo.Height],'color','k','menubar','none','figname',out(i).id);
        hs = zeros(nplots,1);
        for iplot=1:nplots
            hs(iplot) = axes('Units','pixels','position',[(iplot-1)*out(i).nfo.Width 0 out(i).nfo.Width out(i).nfo.Height]);
        end
        id = regexprep(out(i).id,'\_','\\_');
        if o.videoon,
            aviobj = avifile(outputmovie,'fps',o.videofps,'quality',o.videoquality,'compression',o.videocompression); %#ok<TNMLP>
        end
        for iframe=1:out(i).length % for each frame
            for iplot = 1:nplots
                subplot(hs(iplot))
                switch lower(o.image{iplot})
                    case 'im', imshow(out(i).im(:,:,ind0(iframe))), caxis(out(i).thresh(1,:));
                    case 'rgb',  if ~isempty(out(i).rgb), imshow(out(i).rgb(:,:,:,iframe)), end
                    case 'rgb2', if ~isempty(out(i).rgb2), imshow(out(i).rgb2(:,:,:,iframe)), end
                    case 'rgb3', if ~isempty(out(i).rgb3), imshow(out(i).rgb3(:,:,:,iframe)), end
                    case 'scaled'
                        if nplots>1
                            imshow(out(i).scaled(:,:,ind0(iframe))), caxis([0 255]);
                        else
                            set(fig,'colormap',jet(4096))
                            imagesc(out(i).scaled(:,:,ind0(iframe)))
                            caxis([0 255])
                        end
                    otherwise
                        error('''%s'' is an unknown image type',o.image{iplot})
                end
            end % next iplot
            text(out(i).nfo.Width-5,out(i).nfo.Height-10,sprintf('\\color{white}%s: frame %0.3d/%d',id,ind(ind0(iframe)),nfiles),'fontsize',12,'horizontalalignment','right')
            if o.videoon
                screensize = get(0,'ScreenSize'); posfig = get(fig,'position');
                set(fig,'position',[round((screensize(3:4)-posfig(3:4))/2) posfig(3:4)]), drawnow % to prevent: Failed to set stream format. Try specifying a different codec.
                try
                    aviobj = addframe(aviobj,getframe(fig));
                catch erravi
                    dispf('ERROR generating frame %d/%d',iframe,out(i).length)
                    rethrow(erravi)
                end
            elseif o.pauseon && (out(i).length>1)
                screen = dispb(screen,'\tpress a key (use CTRL+C to stop it)'); pause, screen = dispb(screen,'next frame');
            end
        end % next frame
        if o.videoon, ok=close(aviobj); fileinfo(outputmovie), end %#ok<NASGU>
        
       
    end % if nplots
    
    if o.lowmemory, out(i).scaled=[]; out(i).im=[]; end
    
end % next idx
