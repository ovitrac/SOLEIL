function [B,threshout] = homomorphicfilersoleil(A,I,varargin)
%HOMOMORPHICFILTERSOLEIL homomorphic filtering to remove heterogeneous illumination (Butterworth filtering is applied if illumination is unknwon)
%   B = homomorphicfilersoleil(A,I [,'property1','value1',...])
%   [B,thresh] = homomorphicfilersoleil(...)
%       A: (height x width x nframes) intensity image (convert RGB image to grayscale first with rgb2gray)
%       I: (height x width) intensity image to be used as illumination (default = [])
%          when a series of image is used as height x width x nI array, averaging is applied I = mean(I,3)
%          if I is empty, 2 different behaviors are applied:
%               - if nframes==1, an internal method is applied to guess I (see [4])
%               - if nframes==1, A(:,:,max(1,i-1)) is used as reference for A(:,:,i) 
%       list of properties/values:
%           order: filter order (default = 2)
%          cutoff: Lower cut off frequency (default = 10)
%        passband: [0.0999 1.01]
%         prctile: [.1 99.9]
%
%   Additional filtering based on DWFT (check that important frequential information is not removed)
%            dwft: true to activate DWFT on I (default = false)
%                  ATTENTION: setting dwft to true modifies significantly the energy spectrum (exp(cumsum(log(1+out(1).scaled),3)) will diverge)
%         wavelet: DWFT (default = 'db1')
%           level: DWFT level (default = 5)
%          interp: interpolation method for coefficients (default = 'lanczos3'), see imresize (this step remove important high freq)
%
%       
%       Keywords: 'keepclass' to have A and B of the same class
%                 'dwft'    apply a DWFT on I
%
%       B: restored image (default type = 'double')
%
% References:
% [1] http://en.wikipedia.org/wiki/Homomorphism, http://en.wikipedia.org/wiki/Homomorphic_filtering
% [2] http://www.ensta-paristech.fr/~bazeille/fr/TDS07.pdf
% [3] http://appnote.avrportal.com/appnotes/Image-processing/HOMOMORPHIC-PROCESSING-AND-ITS-APPLICATION-TO-IMAGE-ENHANCEMENT-.pdf
% [4] http://www.mathworks.com/matlabcentral/fileexchange/21357-homomorphic-filtering
% [5] image of lena (used in [3]): http://sipi.usc.edu/database/download.php?vol=misc&img=4.2.04
%

% Soleil experiments SUN2011c - 25/02/12 - INRA\Olivier Vitrac - rev. 26/02/12

% revision history
% 26/02/12 add recursion
% 27/02/12 release candidate, ready for moviesoleil (recursion with I=[])
% 28/02/12 fix isindexI for scalars

% default
default = struct('order',2,'cutoff',10,'passband',[0.0999 1.01],'wavelet','db1','level',5,'interp','lanczos3','prctile',[.1 99.9],'dwft',false);
keywordlist = {'keepclass' 'dwft'};

% argcheck
if nargin<1, error('one argument is required'); end
if nargin<2, I = []; end
if ndims(I)>3, error('RGB images are not accepted, I must be a 2D or 3D-array'); end
if ndims(A)>3, error('RGB images are not accepted, A must be a 2D or 3D-array'); end
if ~isempty(I), I = double(I); end, if size(I,3), I = mean(I,3); end
[height,width,nframes] = size(A); classA = class(A); A = double(A);
isindexI = (numel(I)==length(I)); %&& (length(I)==nframes);
if ~isindexI && ~isempty(I) && ~matcmp([height width],size(I)), error('A and I must be of the same size'); end
options = argcheck(varargin,default,keywordlist,'property'); % promote inheriting of property dwft
if length(options.passband)~=2, error('passband must be a 1x2 vector'); end
if options.passband(2)<options.passband(1), options.passband = options.passband([2 1]); end
options.prctile = options.prctile(:)';

% recursion with 3 cases: a) with classes, b) with a set of reference background (aveaging), c) using the previous frame as I
if nframes>1
    if options.keepclass, B = zeros(size(A),classA); end
    screen = dispb('','HOMOMORPHICFILTER initializing...');
    thresh = zeros(nframes,length(options.prctile));
    t0 = clock; 
    if isindexI && ~isempty(I) % I contains a frame index (case a)
        if ~isempty(findduplicates(I)) % duplicates found, I contains classes
            I = I(:)'; [~,jI] = unique(I,'first'); nclass = length(jI); iclass = 0;
            for c = I(sort(jI)) % keep the initial order
                iclass = iclass + 1;
                ind = find(I==c); nind = length(ind);
                for i=1:nind
                    t1 = clock;
                    tmp = homomorphicfilersoleil(A(:,:,ind(i)),A(:,:,ind(1)),options); B(:,:,ind(i)) = tmp;
                    screen = dispb(screen,'HOMOMORPHICFILTER SOLEIL filtered frame %d/%d of class %d (%0.3g)/%d in %0.5g s',i,nind,iclass,c,nclass,etime(clock,t1));
                    thresh(ind(i),:) = prctile(tmp(:),options.prctile);
                end
            end
        else % frame index are not classes, averaging is applied (case b)
            J = mean(A(:,:,I),3); nI = length(I);
            for i=1:nframes
                t1 = clock;
                if ismember(i,I) && nI<10 % 'dwft' must be forced to true (true is added at the end since 'property' is used)
                    tmp = homomorphicfilersoleil(A(:,:,i),J,options,'dwft',true); B(:,:,i) = tmp;
                    screen = dispb(screen,'HOMOMORPHICFILTER SOLEIL filtered frame %d/%d in %0.5g s (background against itself with dwft=true)',i,nframes,etime(clock,t1));
                else
                    tmp = homomorphicfilersoleil(A(:,:,i),J,options); B(:,:,i) = tmp;
                    screen = dispb(screen,'HOMOMORPHICFILTER SOLEIL filtered frame %d/%d in %0.5g s (common background)',i,nframes,etime(clock,t1));
                end
                thresh(i,:) = prctile(tmp(:),options.prctile);
            end
        end
    elseif isempty(I) % no index, the previous frame is used as background (case c)
        for i=1:nframes
            t1 = clock;
            tmp = homomorphicfilersoleil(A(:,:,i),A(:,:,max(1,i-1)),options); B(:,:,i) = tmp;
            screen = dispb(screen,'HOMOMORPHICFILTER SOLEIL filtered frame %d/%d in %0.5g s (background=%d)',i,nframes,etime(clock,t1),max(1,i-1));
            thresh(i,:) = prctile(tmp(:),options.prctile);
        end
        
    else error('recursion fails due to unrecognized background index')
    end
    dispb(screen,'HOMOMORPHICFILTER SOLEIL filtered %d frames in %0.5g s',nframes,etime(clock,t0));
    if nargout>1, threshout = thresh; end
    return
end


% logscale and equivalent DFT
logA = log(1 + A);
fft_logA = fft2(logA); % DFT, do not apply fftshift() as H is symmetric

% Luminance
H = ((options.passband(2)-options.passband(1)).*highpassbutterworth(options.cutoff,height,width,options.order)) + options.passband(1);
if isempty(I)
    % separation of reflectance and luminance based on high pass Butterworth
    filter_fft_logA = (1-H).*fft_logA;
    B = exp(abs(ifft2(filter_fft_logA))); % restored image
else        
    if options.dwft % experimental (to be validated for production, it seems robust)
        [C,S] = wavedec2(I,options.level,options.wavelet);
        aI = appcoef2(C,S,options.wavelet);
        aI = mean(I(:)) * aI / mean(aI(:));
        I = imresize(aI,[height width],options.interp);
    end
    fft_logI = fft2(log(1 + I)); % figure, mesh(20*log10(abs(fftshift(fft_logI))))
    diffA = H.*(fft_logA - fft_logI);
    B = abs(exp(ifft2(diffA))) - 1; % relative variation
end


% class conversion
if options.keepclass
    switch classA
        case 'uint8',  B = uint8(B);
        case 'uint16', B = uint16(B);
        case 'uint32', B = uint32(B);
        case 'single', B = single(B);
        case 'double', B = double(B);
        otherwise, error('class ''%s'' is yet not implemented')
    end
end

if nargout>1, threshout = prctile(B,options.prctile); end

end % end function

%% subfunctions
% Butterworth high pass filter
function H = highpassbutterworth(cutoff,height,width,order)
H=zeros(height,width);
for i=1:height
    for j=1:width
        dist=(((i-height/2).^2+(j-width/2).^2)).^(.5);
        H(i,j)=1/(1+((cutoff/dist)^(2*order)));
    end
end
end

%-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
% the functions below are not used but stored for future evolutions or future readings (see references)
%-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

%Butterworth Bandpass Filter
function filtered_image = butterworthbpf(I,d0,d1,n) %#ok<DEFNU>
% Butterworth Bandpass Filter
% This simple  function was written for my Digital Image Processing course
% at Eastern Mediterranean University taught by
% Assoc. Prof. Dr. Hasan Demirel
% for the 2010-2011 Spring Semester
% for the complete report:
% http://www.scribd.com/doc/51981950/HW4-Frequency-Domain-Bandpass-Filtering
%
% Written By:
% Leonardo O. Iheme (leonardo.iheme@cc.emu.edu.tr)
% 23rd of March 2011
%
%   I = The input grey scale image
%   d0 = Lower cut off frequency
%   d1 = Higher cut off frequency
%   n = order of the filter
%
% The function makes use of the simple principle that a bandpass filter
% can be obtained by multiplying a lowpass filter with a highpass filter
% where the lowpass filter has a higher cut off frquency than the high pass filter.
%
% Usage BUTTERWORTHBPF(I,DO,D1,N)
% Example
% ima = imread('grass.jpg');
% ima = rgb2gray(ima);
% filtered_image = butterworthbpf(ima,30,120,4);

f = double(I);
[nx ny] = size(f);
f = uint8(f);
fftI = fft2(f,2*nx-1,2*ny-1);
fftI = fftshift(fftI);

subplot(2,2,1)
imshow(f,[]);
title('Original Image')
subplot(2,2,2)
fftshow(fftI,'log')
title('Image in Fourier Domain')
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);

for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Create Butterworth filter.
        filter1(i,j)= 1/(1 + (dist/d1)^(2*n));
        filter2(i,j) = 1/(1 + (dist/d0)^(2*n));
        filter3(i,j)= 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end
% Update image with passed frequencies.
filtered_image = fftI + filter3.*fftI;

subplot(2,2,3)
fftshow(filter3,'log')
title('Filter Image')
filtered_image = ifftshift(filtered_image);
filtered_image = ifft2(filtered_image,2*nx-1,2*ny-1);
filtered_image = real(filtered_image(1:nx,1:ny));
filtered_image = uint8(filtered_image);

subplot(2,2,4)
imshow(filtered_image,[])
title('Filtered Image')
end


%% fftshow
function fftshow(f,type)
% Usage: FFTSHOW(F,TYPE)
%
% Displays the fft matrix F using imshow, where TYPE must be one of
% 'abs' or 'log'. If TYPE='abs', then then abs(f) is displayed; if
% TYPE='log' then log(1+abs(f)) is displayed. If TYPE is omitted, then
% 'log' is chosen as a default.
%
% Example:
% c=imread('cameraman.tif');
% cf=fftshift(fft2(c));
% fftshow(cf,'abs')
%
if nargin<2,
    type='log';
end
switch type
    case 'log'
        fl = log(1+abs(f));
        fm = max(fl(:));
        imshow(im2uint8(fl/fm))
    case 'abs'
        fa=abs(f);
        fm=max(fa(:));
        imshow(fa/fm)
    otherwise
        error('TYPE must be abs or log.');
end;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the code below is not used at all (but stored here as it gives very impressive results (recreate functions in separated files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Matlab Source Code of [3]
% clear all;
% close all;
% I = imread(fullfile(find_path_toolbox('confocal'),'lena_4.2.04.tiff'));
% G = rgb2gray(I);
% GM = im2double(G);
% GN = 1-GM;
% GN = im2uint8(GN);
% t = 256;
% ill = illumination_v1(t);
% product = GM.*ill;
% figure;
% imshow(GM);
% title('Original Image');
% figure;
% imshow(product);
% title('Corrupted Image');
% figure;
% kk = imhist(product);
% imhist(product);
% title('Histogram of the corrupted image');
% figure;
% imhist(kk);
% title('Histogram after equalization');
% figure;
% jj = histeq(product);
% imshow(jj);
% title('Image restored after histogram equalization');
% log_GM = log(GM + 1);
% log_product = log(product + 1);
% fft_GM = (fftshift(fft2(log_GM)));
% figure;
% original_value = (20*log10(abs(fft_GM)));
% mesh(20*log10(abs(fft_GM)));
% title('Frequency characteristics original image');
% fft_product = (fftshift(fft2(log_product)));
% figure;
% mesh(20*log10(abs(fft_product)));
% title('Frequency characteristics of illuminated image');
% diff_ill = (fft_product - fft_GM);
% figure;
% mesh(20*log10(abs(diff_ill)));
% title('Frequency characteristics of illumination');
% filt_filt = ifft2(diff_ill);
% figure;
% mesh(20*log10(abs(filt_filt)));
% title('Frequency characteristics of filter');
% restored = (fft_product - diff_ill);
% figure;
% desired_value = (20*log10(abs(restored)));
% mesh(20*log10(abs(restored)));
% title('Frequency characteristics of restored image');
% restored_im = ifft2(restored);
% restored_image = (exp(restored_im)) - 1;
% figure;
% abs_restored_image = abs(restored_image);
% imshow(abs_restored_image);
% title('Restored Image');
% %metrics for testing the performance
% final = original_value - desired_value;
% sum_ave = det(final);
% %error image
% error = restored_image - GM;
% error1 = im2uint8(error);
% figure;
% imshow(error1);
% title('Error image');
% Function file
% function[ill] = illumination_v1(t);
% p = 1/t:1/t:1;
% for k = 1:t,
% ill(k,:) = p(:)';
% end
% figure;
% imshow(ill);
% title('Illumination Pattern');
% p = 1/t:1/t:1;
% for k = 1:t,
% ill(k,:) = p(:)';
% end
% ill = rot90(ill);
% figure;
% imshow(ill);
% title('illumination pattern');
% p = 1/t:1/t:1;
% for k = 1:t,
% il(k,:) = sin(9*p*pi + 0.2);
% end
% for k = 1:t,
% il_1(k,:) = sin(9*p*pi + 0.2);
% end
% il_1 = rot90(il_1);
% ill = il.*il_1;
% figure;
% imshow(ill);
% title('Illumination Pattern');
% p = 1/t:1/t:1;
% for k = 1:t,
% il(k,:) = p(:)';
% end
% figure;
% imshow(il);
% title('Illumination Pattern');
% ill = imrotate(il,-45,'bilinear','crop');
% imshow(ill);
% title('Rotated Illumination Pattern');
% dia = sqrt(2*(t^2));
% p = 1/dia:1/dia:1;
% for k = 1:dia,
% il(k,:) = p(:)';
% end
% size_il = size(il)
% figure;
% imshow(il);
% title('Illumination Pattern');
% J = imrotate(il,-45,'bilinear','crop');
% ill = imcrop(J,[0,0,t,t]);
% size_ill = size(ill)
% figure;
% imshow(ill);
% title('Rotated Illumination Pattern');
% p = 1/t:1/t:1;
% for k = 1:t,
% il(k,:) = exp(p*pi/10);
% end
% for k = 1:t,
% il_1(k,:) = exp(p*pi/10);
% end
% il_1 = rot90(il_1);
% ill = il.*il_1;
% figure;
% imshow(ill);
% title('Illumination Pattern');