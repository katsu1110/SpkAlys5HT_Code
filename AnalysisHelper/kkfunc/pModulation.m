function pmodulation = pModulation( ex0, ex2, exinfo, p_flag )
% pmodulation = pModulation( ex0, ex2, exinfo )
%
% this functions computes a distribution of spike counts for two
% experiments. the spikes were normalized for each stimulus
% conditions in the merged experiments, therefore allowing to factor out
% the stimulus dependent variation.
%
% pmodulation is a vector containing the p-values of the parametric and
% non-parametric test of two distributions.
%
%
%
% @CorinnaLorenz 12/7/2017
%


if nargin==3
    p_flag = false;
end


i_base = length(ex0.Trials);
ex_all = ex0;
fieldname0 = fieldnames(ex0.Trials);
fieldname2 = fieldnames(ex2.Trials);

k0 = ~ismember(fieldname0, fieldname2);
ex0.Trials = rmfield(ex0.Trials, fieldname0(k0));

k2 = ~ismember(fieldname2, fieldname0);
ex2.Trials = rmfield(ex2.Trials, fieldname2(k2));


if size(ex0.Trials, 1)>size(ex0.Trials, 2)
    ex_all.Trials = [ex0.Trials; ex2.Trials];
else
    ex_all.Trials = [ex0.Trials, ex2.Trials];
end
ex_all=znormex(ex_all, exinfo);

znormed0 = [ex_all.Trials(1:i_base).zspkcount];
znormed2 = [ex_all.Trials(i_base+1:end).zspkcount];
p_nonparam = ranksum(znormed0, znormed2);
[~,p_param] = ttest2(znormed0, znormed2);
pmodulation = [p_param,p_nonparam];


if p_flag
    figure;
    h0=histogram(znormed0,10, 'FaceColor', lines(1));hold on;
    h2=histogram(znormed2, 'FaceColor', getCol(exinfo));
    h2.BinWidth = min([h0.BinWidth, h2.BinWidth]);
    h0.BinWidth = min([h0.BinWidth, h2.BinWidth]);
    ylim = get(gca, 'YLim');
    plot(ones(1,2)*mean(znormed0), ylim, 'Color', h0.FaceColor, 'LineWidth', 2);
    plot(ones(1,2)*mean(znormed2), ylim, getCol(exinfo), 'LineWidth', 2);
    
    fname0=getFname(ex0);
    fname2=getFname(ex2);
    legend(fname0(15:end), fname2(15:end), 'mean', 'mean');
    xlabel('z scored spike count responses');
    
    title(sprintf('p_{param}=%1.2e, p_{nonparam}=%1.2e', ...
        p_param, p_nonparam));
    
    
end

end

