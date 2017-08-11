function exinfo = setReceptiveFieldSize( exinfo )
% exinfo = setReceptiveFieldSize( exinfo )
% 
% searches the XPos and YPos files for each session and assigns the
% equivalent width of the receptive field in the x and y dimension to all
% corresponding entries
% 
% @CL

unit_id = unique([exinfo.id]);

for i =1:length(exinfo)
    exinfo(i).RFwx = nan;
    exinfo(i).RFwy = nan;
    exinfo(i).RFw = nan;
end

% loop through all unit recordings
for i = 1:length(unit_id)

    wY = nan; wX = nan;
    fprintf('working on unit : %1.1f \n', unit_id(i));
    
    %find all entries belonging to this session
    idx = find([exinfo.id]==unit_id(i));
    
    % XPos
    [fnames, fdir] = getExFileNames(unit_id(i));
    fnamesX = fnames( contains(fnames, 'XPos') );
    fnamesX = fnamesX( contains(fnamesX, 'c1') & contains(fnamesX, 'sortLH') );
    
    % if there was no RF experiment for this unit, than look for the
    % preceding unit recordings
    k = 1;
    while isempty(fnamesX) 
        [fnames, fdir] = getExFileNames(unit_id(i)-k);
        fnamesX = fnames( contains(fnames, 'XPos') );
        fnamesX = fnamesX( contains(fnamesX, 'c1') & contains(fnamesX, 'sortLH') );
        k = k+1;        
    end
    
    if ~isempty(fnamesX)
        for j =1:length(fnamesX)
            ex = loadCluster(fullfile(fdir, fnamesX{j}), 'loadlfp', false);
            wX(j) = getMarginalDist(ex.Trials, 'x0');
        end
    end
 
    % YPos
    [fnames, fdir] = getExFileNames(unit_id(i));
    fnamesY = fnames( contains(fnames, 'YPos') );
    fnamesY = fnamesY( contains(fnamesY, 'c1') & contains(fnamesY, 'sortLH') );

    k = 1;
    while isempty(fnamesY) 
        [fnames, fdir] = getExFileNames(unit_id(i)-k);
        fnamesY = fnames( contains(fnames, 'YPos') );
        fnamesY = fnamesY( contains(fnamesY, 'c1') & contains(fnamesY, 'sortLH'));
        k = k+1;        
    end
        
    if ~isempty(fnamesY)
        for j =1:length(fnamesY)
            ex = loadCluster(fullfile(fdir, fnamesY{j}), 'loadlfp', false);
            wY(j) = getMarginalDist(ex.Trials, 'y0');
        end
    end
    
    
    %========================================== assign receptive field size
    if isempty(fnamesY) && isempty(fnamesX) 
        
        if idx(1)>1 && (exinfo(idx(1)-1).date - exinfo(idx(1)).date) <0.8 
            wX = exinfo(idx(1)-1).RFwx;
            wY = exinfo(idx(1)-1).RFwy;
            
            for j = 1:length(idx)
                exinfo(idx(j)).RFwx = nanmean(wX);
                exinfo(idx(j)).RFwy = nanmean(wY);
            end
            
        else
            continue;
        end
    end
    
    for j = 1:length(idx)
        exinfo(idx(j)).RFwx = nanmean(wX);
        exinfo(idx(j)).RFwy = nanmean(wY);
        
        if isnan(exinfo(idx(j)).RFwy)
            exinfo(idx(j)).RFw = exinfo(idx(j)).RFwx;
            
        elseif isnan(exinfo(idx(j)).RFwx)
            exinfo(idx(j)).RFw = exinfo(idx(j)).RFwy;
            
        elseif ~isnan(exinfo(idx(j)).RFwx) && ~isnan(exinfo(idx(j)).RFwy)
            exinfo(idx(j)).RFw = ...
                sqrt(exinfo(idx(j)).RFwy.^2 + exinfo(idx(j)).RFwx^2);
            
        end
    end

           
    rf(i) = exinfo(idx(j)).RFw;
    ecc(i) = exinfo(idx(j)).ecc;
    clearvars wX wY
        
end



exinfo = setRFcorrected(exinfo, rf, ecc);
end

%%
function exinfo = setRFcorrected(exinfo, rf, ecc)
% computes the eccentricity independent receptive field width via linear
% regression 

ecc= ecc(~isnan(rf))'; rf = rf(~isnan(rf))';
tbl = table(ecc, rf);
lm = fitlm(tbl, 'rf~ecc');
lm.Coefficients.Estimate

for i =1 :length(exinfo)
   
    if ~isnan(exinfo(i).RFw)
        exinfo(i).RFw_corr = exinfo(i).RFw - ...
            (lm.Coefficients.Estimate(1)+exinfo(i).ecc * lm.Coefficients.Estimate(2));
    else
         exinfo(i).RFw_corr  = nan;
    end
end



end

%%
function w = getMarginalDist(trials, posname)
% equivalent width of the marginal distribution of spike rate to different
% bar position

trials = trials([trials.Reward]==1);
trials = trials([trials.(posname)]< 1000);


% bar positions
pos = unique([trials.(posname)]);

for i = 1:length(pos)
    tp = trials([trials.(posname)] == pos(i));
    meanspk(i) = nanmean([tp.spkRate]);
    sdnspk(i) = nanstd([tp.spkRate]);    
    nrep(i) = length(tp);
end


% if the mapping is not causing a selective response return with nan
if anova1([trials.spkRate], [trials.(posname)], 'off') > 0.08 
    % for debugging:
%         figure; errorbar(pos, meanspk, sdnspk);
%         text(pos, meanspk, num2str(nrep'));
    w = nan;
else
    
    % substract baseline firing rate
    meanspk = meanspk - min(meanspk);
    meanspk = meanspk/max(meanspk); % a redundant step...

    % area / height of the gaussian like curve
    A = sum(meanspk); h = max(meanspk);
    w = A / h;
    w = w* mean(diff(pos(pos<1000))); % normalize to the given unit
    
    
    % for debugging
%     figure;    plot(pos, meanspk, 'Displayname', posname); hold on
%     text(pos(3), meanspk(3), num2str(w));
end
end

%%
function nr = foldernr(session)
% prefixes zeros to the session number to get the correct foldername

s = num2str(session);   prefix = [];

if length(s)==1;        prefix = '000';
elseif length(s)==2;    prefix = '00';
elseif length(s)==3;    prefix = '0';       
end

nr = [prefix s];
end

function d = getD(ex)

if isfield(ex.Trials, 'x0')
    uq = unique([ex.Trials.x0]);
else
    uq = unique([ex.Trials.y0]);
end
d = diff(uq(1:2));
end


function [fnames, fdir] = getExFileNames(unitnr)

if mod(unitnr,1) == 0
    fdir = ['Z:\data\mango\' foldernr(unitnr)];
else
    fdir = ['Z:\data\kaki\' foldernr(unitnr-0.5)];
end

%get all filenames belonging to this session
fnames = dir( fdir );    fnames = {fnames.name};
end

