function exinfo = setPhaseSelTFexp( exinfo )
% exinfo = setPhaseSelTFexp( exinfo )
% 
% loops through all units oinf exinfo and tries to find experiments of them
% with varying temporal frequency in order to compute the phase selectivity
% (f1/f0).
%
%
% @CL



unit_id = unique([exinfo.id]); % all units within exinfo

% loop through all units
for i = 1:length(unit_id)
    
    fprintf('working on sessison : %1.1f \n', unit_id(i));
    %find all entries belonging to this session
    idx = find([exinfo.id]==unit_id(i));
    
    %get all tf files belonging to this session
    if exinfo(idx(1)).ismango
        fdir = ['Z:\data\mango\' foldernr(unit_id(i))];
    else
        fdir = ['Z:\data\kaki\' foldernr(unit_id(i))]; 
    end
    
    % get all filenames within the directory of this unit
    fnames = dir( fdir );    fnames = {fnames.name};
    fnames = fnames( contains(fnames, 'TF') ); % all experiments with variation of temporal frequency 
    
    if ~isempty(fnames)
        fnames = fnames( contains(fnames, 'c1') );
        f1f0 = getF1F0(fdir, fnames); % always plots the results
        
        % save the result plot 
        h = findobj('type', 'figure'); 
        set(h, 'Name', exinfo(idx(1)).figname);
        savefig(h , exinfo(idx(1)).fig_phasetf); 
        close(h);
    else
        f1f0 = nan;
    end
    
    % assign f1f0 to all experiments of this unit in exinfo
    for j = 1:length(idx)
        exinfo(idx(j)).tf_f1f0 = f1f0;
        exinfo(idx(j)).fig_phasetf = exinfo(idx(1)).fig_phasetf;
    end
end


function f1f0 = getF1F0(fdir, fnames)

% concatenate trials if the experiment was repeated
load(fullfile(fdir, fnames{1}), 'ex'); ex0 = ex;
for i =2:length(fnames)
    load(fullfile(fdir, fnames{i}), 'ex');
    ex0.Trials = [ex0.Trials, ex.Trials];
end

f1f0 = getPhaseSelectivity(ex0, 'plot', true);



function nr = foldernr(unitID)
% returns the file folder corresponding to the unit

s = num2str(floor(unitID));  % floor() accounts for kaki's id being +0.5
prefix = [];

if length(s)==1;        prefix = '000';
elseif length(s)==2;    prefix = '00';
elseif length(s)==3;    prefix = '0';       
end

nr = [prefix s];
