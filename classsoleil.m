function [class,pdistclass] = classsoleil(sample,nclass)
%CLASSSOLEIL split time frames in clusters
%   default method based on KMEANS
%       syntax: class = samplesoleil(sample [,nclass])
%       sample: sample produced by samplesoleil (best results with resolution=4)
%       nclass: number of classes (default = 3)
%        class: cluster index
%   method based on PDIST (KMEANS result is also included)
%       syntax: [class,pdistclass] = samplesoleil(...)
%        class: as above
%   pdistclass: structure with fields class and linkage (note that internally nclass+1 is applied)
%   
% Example
%    db = loaddbsoleil('local','~/SUN_interpretation/Telemos/'); 
%    sample = samplesoleil(db,'series',13,'step',1,'resolution',4);
%    class = classsoleil(sample,3);
%    figure, hist(class,1:3)
%   [class,pdistclass] = classsoleil(sample,5)
%   figure, dendrogram(pdistclass.linkage,'colorthreshold','default')
%
%
% See also: samplesoleil, immean, loaddbsoleil, loadansoleil, moviesoleil, synchrotronpattern, synchrotronpattern2

% Soleil experiments SUN2011c SUN2011d - 12/01/12 - INRA\Olivier Vitrac - rev. 13/01/12

% revision history
% 12/01/12 fix definition of class, add
% 13/01/12 add link

    % default
    nclass_default = 3;

    % arg check
    if nargin<1, error('one argument is required'); end
    if ~isstruct(sample) || ~isfield(sample,'im') || isscalar(sample.im) || ndims(sample.im)~=3 || size(sample.im,3)<10
        error('sample must be created with SAMPLESOLEIL, e.g.: sample = samplesoleil(db,''series'',[an(1).series],''step'',3,''resolution'',4)');
    end
    if nargin<2, nclass = []; end
    if isempty(nclass), nclass = nclass_default; end

    % clustering
    siz = size(sample.im);
    t0 = clock;
    screen = dispb('','CLASSSOLEIL\t running on %dx%dx%d...',siz(1),siz(2),siz(3));
    tmp = reshape(sample.im,[siz(1)*siz(2) siz(3)])'; avgtmp = mean(tmp,2);
    class = sortclass(kmeans(tmp,nclass,'distance','sqEuclidean'),nclass);
    dispb(screen,'CLASSSOLEIL\t completed in %0.5g s',etime(clock,t0));


    % add link
    if nargout>1
        link = linkage(pdist(tmp,'euclidean'),'average');
        c = sortclass(cluster(link,'maxclust',nclass+1),nclass+1);
        pdistclass = struct('class',c,'linkage',link);        
    end

    % internal function
    % function sorting with increasing intensity
    function finalcluster = sortclass(rawcluster,ncluster)        
        avgc = zeros(ncluster,1);
        for i=1:ncluster
            avgc(i) = mean(avgtmp(rawcluster==i));
        end
        [~,iavgc] = sort(avgc);
        reorder(iavgc) = 1:ncluster;
        finalcluster = reorder(rawcluster);
    end

end % end function
