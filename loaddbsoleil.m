function db = loaddbsoleil(varargin)
%LOADDBSOLEIL load a SOLEIL database, based on SUN2011c (i.e. ODS file with useful fields, to be incremented)
%   syntax: db = loaddbsoleil()
%   syntax: db = loaddbsoleil('param1',value1,'param2',value2...)
%       list of implemented properties:
%              root: raw data location (parent folder that contains raw images)
%             local: working directory (interpretation folder)
%            dbfile: database file (default = 'description_results.ods')
%excesspathtoremove: relocation of bad paths (see code for details)
%         sheetname: working sheetname (default = 'id')
%       filepattern: *.tif
%      prefetchfile: (default ='')
%        noprefetch: flag for prefetch (default=false)
%                    can be used also as a keyword
%
% 	OUUTPUT db = structure with fields (example below):
%         folder1: 'SUN2011d'
%         folder2: 'day2'
%         folder3: '54PERM'
%         folder4: NaN
%         folder5: NaN
%         folder6: NaN
%      experiment: 54
%        typeflag: 2
%      matrixflag: 2
%          frames: 300
%            size: 1024
%     pixellength: 1.2402
%           zstep: 0
%        isseries: 1
%          isperm: 1
%         isstack: 0
%          iswalk: 0
%          format: 'permeation'
%         product: 'French-fries'
%         SUM2012: 1
%         comment: 'one_cell_filled'
%        fullpath: 'C:\Data\Olivier\INRA\Projects\SUN\SUN2011d\day2\54PERM\'
%         isvalid: 1
%     framecounts: 300
%           files: {300x1 cell}
%         imfinfo: [1x1 struct]
%           batch: [1x1 struct]
%          nbatch: 1
%
%   Output for montagesoleil() is stored as:
%       db(iseries).nbatch = number of natch
%       db(iseries).batch(ibatch).nobject = number of objects (polygons)
%       db(iseries).batch(ibatch).object(iobject) with fields (example below)
%                    ind: [250 251 252 253 254 255 256]
%                 indall: [1x67 double]
%             background: [234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249]
%               realtime: [-15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51]
%                    roi: [590 255 256 256]
%                polygon: [7x2 double]
%   
%
%
% See also: moviesoleil, montagesoleil

% Soleil experiments SUN2011c-SUN2011d - 17/11/11 - INRA\Olivier Vitrac - rev. 12/05/12

% Revision history
% 03/01/12 add prefetchfile, first release candidate since run2011d
% O3/01/12 parameters 'select' and 'priority' removed to avoid the change in the number and order of lines
%          this modification is set to work in parallel with another ods file for image analyses
% 04/01/12 concatenate an arbitrary number of folder1, folder2, folder3..foldern (regarless their order)
% 13/01/12 fix o_default.root for LP-MOL5
% 26/01/12 use loadodsprefetch instead of loadods
% 26/01/12 help updated, add fullpath in dbfile
% 26/01/12 add files and imfinfo
% 27/01/12 sort files
% 27/01/12 send to linux
%          for ip in 19 20 21 22 23 24 36 27 28 111 112 113 114; do scp -p /home/olivier/codes/confocal/loaddbsoleil.m olivier@10.75.4.$ip:/home/olivier/codes/confocal/. ; done
%          scp -p /home/olivier/codes/confocal/loaddbsoleil.m Olivier@10.75.4.102:/cygdrive/d/Data/Olivier/INRA/Codes/confocal/.
%          scp -p /home/olivier/codes/confocal/loaddbsoleil.m Olivier@10.75.4.115:/cygdrive/c/Data/Olivier/INRA/Codes/confocal/.
% 13/04/12 remove the last sep in the regexp defining the folder paths of series (solve bug with consecutive NaN that are not removed when no sep are present)
% 06/05/12 read and parse batch data (Olivier)
% 12/05/12 batchfields: remove realtime, add matrix t0 (removed by JMV in description_results.ods)
% 14/05/12 fix parsing (retrieved from 10/05/12)
%          tmp = cellfun(@(f) db(i).(sprintf('%s%d',f,nbatch)),batchfields,'UniformOutput',false); tmp = cell2struct(tmp,batchfields,2);
%          tmp = cellfun(@(f) db(i).(sprintf('%s%d',f,j)),batchfields,'UniformOutput',false); tmp = cell2struct(tmp,batchfields,2);

% default
o_default = struct(...
    'root','/home/olivier/SUN_rawdata_donotmodify/Telemos',... raw data location (machine specific)
    'local','',... working directory
    'excesspathtoremove','',... % workaround for either relocation or bad paths (keep it empty if OK with paths), to be fixed by JMV/AP
    'dbfile','description_results.ods',...
    'sheetname','id',... %'all'
    'suggest','suggest',... % worksheet used to put data that sort/select results according to priorities
    'filepattern','*.tif',...       %parameter removed 'priority',false,... %set 'priority' as 'true' to sort database
    'prefetchfile','',...
    'noprefetch',false...
    );
switch localname % user overrides for some machines
    case 'WSLP-OLIVIER2',                       o_default.root = 'D:\Data\Olivier\INRA\Projects\SUN\';
    case {'WSLP-OLIVIER3' 'LP-MOL6' 'LP-MOL7'}, o_default.root = 'C:\Data\Olivier\INRA\Projects\SUN\';
    case 'LP-MOL5',                             o_default.root = 'C:\Data\SUN_rawdata_donotmodify\Telemos\';
    case 'E075016',                             o_default.root = 'D:\Documents and Settings\jevauvre\My Documents\PhD\SUN_rawdata_donotmodify\Telemos\';
end
kwlist = 'noprefetch';
% batch definitions
nbatchmax = 10; % number maximal of batches
% batchfields = {'ind','indall','background','realtime','roi','polygon'}; % fields used in batch (followed in ODS file by an index ranging 1..nbatch)
% realtime removed 12/5/2012
batchfields = {'ind','indall','background','roi','polygon','matrix','t0'}; % fields used in batch (followed in ODS file by an index ranging 1..nbatch)
batchiscell = cell2struct(repmat({true},length(batchfields),1),batchfields,1); %batchiscell.polygon = true;

% argcheck
o = argcheck(varargin,o_default,kwlist,'case');
if isempty(o.prefetchfile)
    o.prefetchfile = sprintf('PREFETCH_%s',regexprep(o.dbfile,'.ods$','.mat','ignorecase'));
end
t0 = clock;

if isempty(o.local), [o.local,odsname,odsext]=fileparts(o.dbfile);  o.dbfile = sprintf('%s%s',odsname,odsext); end
if isempty(o.local),
    switch localname
        case 'WSLP-OLIVIER2',                       o.local = 'D:\Data\SUN_interpretation\Telemos\';
        case {'WSLP-OLIVIER3' 'LP-MOL6' 'LP-MOL7'}, o.local = 'C:\Data\SUN_interpretation\Telemos\';
        case 'LP-MOL5',                             o.local = 'C:\Data\SUN_interpretation\Telemos\';
        case 'E075016',                             o.local = 'D:\Documents and Settings\jevauvre\My Documents\PhD\SUN_interpretation\Telemos\';
        otherwise, o.local = '/home/olivier/SUN_interpretation/Telemos/';
    end
end

% PREFETCH / LOAD
% Before 26/01/12
% if ~o.noprefetch && exist(fullfile(o.local,o.prefetchfile),'file')
%     dispf('Load prefetchfile...')
%     fileinfo(fullfile(o.local,o.prefetchfile))
%     load(fullfile(o.local,o.prefetchfile))
% else
%     db = loadods(fullfile(o.local,o.dbfile),'sheetname',o.sheetname,'structarray',true,'forceboolean',true);
%     save(fullfile(o.local,o.prefetchfile),'db')
%     dispf('...prefetch file updated')
%     fileinfo(fullfile(o.local,o.prefetchfile))
% end
%
% After 26/01/12
db = loadodsprefetch(fullfile(o.local,o.dbfile),'sheetname',o.sheetname,'structarray',true,'forceboolean',true,'prefetchpath',o.local);

% full paths: concatenate an arbitrary number of folder1, folder2, folder3..foldern (regarless their order)
% NaN/ or NaN\ values (blank cells) are removed
screen = dispb('','LOADBSOLEIL: create full paths...');
nseries = length(db);
folderdepth = regexp(fieldnames(db),'^folder(\d+)$','tokens'); folderdepth = sort(str2double(uncell(folderdepth(~cellfun('isempty',folderdepth)))));
folders = cellfun(@(f){db.(f)},arrayfun(@(x) sprintf('folder%d',x),folderdepth,'UniformOutput',false),'UniformOutput',false);
fullpath = cell(nseries,1); sep = regexprep(filesep,'\\','\\\\'); folders = cat(1,folders{:});
for iseries=1:nseries
    fullpath{iseries} = fullfile(o.root,regexprep(sprintf(['%s' sep],folders{:,iseries}),{['^NaN' sep '|' sep 'NaN'],o.excesspathtoremove},{'',''},'ignorecase'));
end
[db.fullpath] = deal(fullpath{:});

% check that all files exist (workaround as some paths do no exist, to be fixed by JMV/AP)
screen = dispb(screen,'LOADBSOLEIL: check that all files exist...');
isvalid = true(length(db),1);
isvalid([db.isseries]) = cellfun(@(f)exist(f,'dir'),{db([db.isseries]).fullpath});
isvalid(~[db.isseries]) = cellfun(@(f)exist(f,'file'),{db(~[db.isseries]).fullpath});
if ~all(isvalid), warning('some foldernames are not valid, see below'), cellfun(@(f) dispf('\t ''%s'' is nonexistent or not a directory',f),{db(~isvalid).fullpath}), screen = ''; end %#ok<WNTAG>
isvalid = num2cell(isvalid,2);
[db.isvalid] = deal(isvalid{:});

% count the number of images for all series (workaround as some paths do no exist, to be fixed by JMV/AP)
screen = dispb(screen,'LOADBSOLEIL: counts the number of frames for each series...');
[db(~[db.isvalid]).framecounts] = deal(0);
[db([db.isvalid]).framecounts] = deal(1); % by default
validseries = [db.isseries] & [db.isvalid];
files = cellfun(@(f)sort(explore(o.filepattern,f,[],'fullabbreviate')),{db(validseries).fullpath},'UniformOutput',false);
[db(validseries).files] = deal(files{:});
framecounts = num2cell(cellfun(@(f)length(f),{db(validseries).files}));
[db(validseries).framecounts] = deal(framecounts{:});
for i=1:length(db), if validseries(i)==1, db(i).imfinfo = imfinfo(char(db(validseries(i)).files(1))); end, end
emptyseries = [db.framecounts]<1;
if any(emptyseries), warning('some (time/stack) series are not valid'), cellfun(@(f) dispf('\t ''%s'' does not contain frames',f),{db(emptyseries).fullpath}), screen = ''; end %#ok<WNTAG>

% batch identification/counts the maximum number of batches
ok = true; nbatch = 0;
while (nbatch<nbatchmax) && ok
    ok = all(isfield(db,cellfun(@(f) sprintf('%s%d',f,nbatch+1),batchfields,'UniformOutput',false)));
    if ok, nbatch = nbatch + 1; end
end
nbatchmax = nbatch;

% ------------- BATCH EXTRACTION AND PARSING by INRA\Olivier - 06/05/11 -------------
% batch identification/extracts batches
if nbatchmax>0
    screen = dispb(screen,'Prepare batch data...');
    for i=1:length(db)
        
        t1 = clock;
        % check all batch fields
        % set to [] any batchfield with NaN
        for f=batchfields
            for nbatch = 1:nbatchmax
                fb = sprintf('%s%d',f{1},nbatch);
                if isnumeric(db(i).(fb)) && ~isempty(db(i).(fb)) && isnan(db(i).(fb)), db(i).(fb)=[]; end
                if ischar(db(i).(fb)) && ~isempty(db(i).(fb)) % eval strings
                    try
                        db(i).(fb)=eval(db(i).(fb));
                    catch evalerr % generate an error if unable to evaluate the string
                        dispf('ERROR in evaluating field ''%s'' of db row %d/%d',fb,i,length(db))
                        dispf('\texpression=''%s''',db(i).(fb))
                        rethrow(evalerr)
                    end
                    if batchiscell.(f{1}) && ~iscell(db(i).(fb)), db(i).(fb) = {db(i).(fb)}; end % force cell if required
                end % endif string
            end % next nbatch value
        end % next batchfield
        ok = true; nbatch = 0;
        while (nbatch<nbatchmax) && ok
            ok =  ~all(cellfun(@(f) isempty(db(i).(sprintf('%s%d',f,nbatch+1))),batchfields));
            if ok, nbatch = nbatch + 1; end
        end
        
        % parse object data
        batch = repmat(struct('object',[],'nobject',0),nbatch,1); % prefetch output
        for j=1:nbatch % for each batch
            tmp = cellfun(@(f) db(i).(sprintf('%s%d',f,j)),batchfields,'UniformOutput',false); tmp = cell2struct(tmp,batchfields,2);
            npol = structfun(@length,tmp); npolmax = max(npol); npol = cell2struct(num2cell(npol,2),batchfields,1);
            object = repmat(cell2struct(repmat({[]},length(batchfields),1),batchfields,1),npolmax,1); % prefetch objects
            for f=batchfields % for each field in batch
                for k=1:npolmax % for each object
                    if k<=npol.(f{1}) % if field defined for current object
                        object(k).(f{1}) = tmp.(f{1}){k};
                    elseif npol.(f{1})>0 % propagate last object
                        object(k).(f{1}) = tmp.(f{1}){npol.(f{1})};
                    end 
                end % next object
            end % next field
            batch(j).object = object;
            batch(j).nobject = length(object);
        end
        db(i).batch = batch; % add batch
        db(i).nbatch = nbatch;
        screen = dispb(screen,'Batch data of %d/%d parsed in %0.5g s',i,length(db),etime(clock,t1));
        
    end % next db row
    
else % no batch
    [db(i).batch]  = deal([]);
    [db(i).nbatch] = deal(0);

end % if any batch defined in the database

% clean batch fields
flist = fieldnames(db);
for f=batchfields
    db = rmfield(db,flist(~cellfun(@isempty,regexp(flist,[f{1} '\d+']))));
end

% final
dispb(screen,'LOADBSOLEIL: read %d series in %0.3g s',nseries,etime(clock,t0));