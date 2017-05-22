function exinfo = setPhaseSelTFexp( exinfo )



session = unique([exinfo.id]);

% loop through all sessions
for i = 1:length(session)
    
    fprintf('working on sessison : %1.1f \n', session(i));
    %find all entries belonging to this session
    idx = find([exinfo.id]==session(i));
    
    %get all tf files belonging to this session
    if exinfo(idx(1)).ismango
        fdir = ['Z:\data\mango\' foldernr(session(i))];
    else
        fdir = ['Z:\data\kaki\' foldernr(session(i)-0.5)];
    end
    
    fnames = dir( fdir );    fnames = {fnames.name};
    fnames = fnames( cellfun(@(x) ~isempty(strfind(x, 'TF')), fnames) );
    
    if ~isempty(fnames)
        fnames = fnames( cellfun(@(x) ~isempty(strfind(x, 'c1')), fnames) );
        f1f0 = getF1F0(fdir, fnames);
        
        h = findobj('type', 'figure'); 
        set(h, 'Name', exinfo(idx(1)).figname);
        savefig(h , exinfo(idx(1)).fig_phasetf); 
        close(h);
    else
        f1f0 = nan;
    end
    
    % assign f1f0
    for j = 1:length(idx)
        exinfo(idx(j)).tf_f1f0 = f1f0;
        exinfo(idx(j)).fig_phasetf = exinfo(idx(1)).fig_phasetf;
    end
end


function f1f0 = getF1F0(fdir, fnames)

load(fullfile(fdir, fnames{1}), 'ex'); ex0 = ex;
for i =2:length(fnames)
    load(fullfile(fdir, fnames{i}), 'ex');
    ex0.Trials = [ex0.Trials, ex.Trials];
end

f1f0 = getPhaseSelectivity(ex0, 'plot', true);



function nr = foldernr(session)

s = num2str(session);   prefix = [];

if length(s)==1;        prefix = '000';
elseif length(s)==2;    prefix = '00';
elseif length(s)==3;    prefix = '0';       
end

nr = [prefix s];
