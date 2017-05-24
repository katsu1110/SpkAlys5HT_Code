function exinfo = evalBothEx(ex0, ex2, exinfo, p_flag, varargin)
%evalBothEx evaluates operations on both base (ex0) and drug (ex2) files



%% preferred/unpreferred stimuli occuring in both files
idx = ismember(exinfo.ratepar, exinfo.ratepar_drug) &  exinfo.ratepar<1000;
parb = exinfo.ratemn;

parb(~idx, :) = 0;
[~, exinfo.pfi] = max( parb );
exinfo.pfi_drug = find(exinfo.ratepar_drug == exinfo.ratepar( exinfo.pfi ));

parb(~idx) = 10^4;
[~, exinfo.upfi] = min( parb );
exinfo.upfi_drug = find(exinfo.ratepar_drug == exinfo.ratepar( exinfo.upfi ));


%% non-parametric change, i.e. relative change in the integrated area under the tc curve
exinfo = bootstrap_exinfo( exinfo );

%% gain change
if any(strcmp(varargin, 't2reg')) || any(strcmp(varargin, 'all'))
    [gslope, yoff, r2, bootstrp] = type2reg(exinfo, p_flag);
    
    exinfo.gslope = gslope(1);
    exinfo.yoff = yoff(1);
    exinfo.r2reg = r2(1);
    
    exinfo.yoff_rel = yoff(2);
    %exinfo.reg_bootstrp = bootstrp;
end

%% fitting parameters
if any(strcmp(varargin, 'tcfit')) || any(strcmp(varargin, 'all'))
    if isfield(exinfo.fitparam, 'others') && isfield(exinfo.fitparam.others, 'OR')
        
        %%% orientation
        for i = 1:length(exinfo.fitparam.others.OR)
            
            cont = exinfo.fitparam.others.OR(i).val.mn;
            cont = [cont; exinfo.ratemn(exinfo.ratepar>180)];
            drug = exinfo.fitparam_drug.others.OR(i).val.mn;
            drug = [drug; exinfo.ratemn_drug(exinfo.ratepar_drug>180)];
            [beta0, beta1, xvar] = ...
                perpendicularfit(cont, drug, var(drug)/var(cont));
            
            exinfo.fitparam.others.OR(i).yoff = beta0;
            exinfo.fitparam.others.OR(i).gslope = beta1;
            exinfo.fitparam.others.OR(i).r2reg = xvar;
        end
        
        %%% contrast
        cont = exinfo.fitparam.others.CO.val;
        drug = exinfo.fitparam_drug.others.CO.val;
        
        [beta0, beta1, xvar] = ...
            perpendicularfit(cont, drug, var(drug)/var(cont));
        
        exinfo.fitparam.others.CO.yoff = beta0;
        exinfo.fitparam.others.CO.gslope = beta1;
        exinfo.fitparam.others.CO.r2reg = xvar;
        
        plotCOTC(exinfo);
        
        
    elseif strcmp(exinfo.param1, 'sf')      %%% spatial frequency
        
        % first entry linear fit, second entry log scaled fit
        others_base = exinfo.fitparam.others;
        others_drug= exinfo.fitparam_drug.others;
        if others_base{1}.r2 > others_base{2}.r2 && others_base{1}.mu > 0
            exinfo.fitparam = others_base{1};
            exinfo.fitparam_drug = others_drug{1};
            exinfo.fitparam.type = 'linear';
        else
            exinfo.fitparam = exinfo.fitparam.others{2};
            exinfo.fitparam_drug = exinfo.fitparam_drug.others{2};
            exinfo.fitparam.type = 'log';
        end
        exinfo.fitparam_drug.others = others_drug;
        exinfo.fitparam.others = others_base;
        
    end
    
    
    %% pfreferred orientation fit evaluation
    
    if strcmp(exinfo.param1, 'or');
        if isfield(exinfo.fitparam, 'OR')
            
            for i = 1:length(exinfo.fitparam.OR)
                
                if abs(exinfo.fitparam.OR(i).mu - exinfo.fitparam_drug.OR(i).mu) > 150
                    if exinfo.fitparam.OR(i).mu < exinfo.fitparam_drug.OR(i).mu
                        exinfo.fitparam.OR(i).mu = exinfo.fitparam.OR(i).mu +180;
                    elseif exinfo.fitparam_drug.OR(i).mu < exinfo.fitparam.OR(i).mu
                        exinfo.fitparam_drug.OR(i).mu = exinfo.fitparam_drug.OR(i).mu +180;
                    end
                end
            end
            
        else
            
            if abs(exinfo.fitparam.mu - exinfo.fitparam_drug.mu) > 150
                if exinfo.fitparam.mu < exinfo.fitparam_drug.mu
                    exinfo.fitparam.mu = exinfo.fitparam.mu +180;
                elseif exinfo.fitparam_drug.mu < exinfo.fitparam.mu
                    exinfo.fitparam_drug.mu = exinfo.fitparam_drug.mu +180;
                end
            end
        end
    end
    
end


%% wave form
if any(strcmp(varargin, 'wavewdt')) || any(strcmp(varargin, 'all'))
    
    if isfield(ex0.Trials, 'Waves') && isfield(ex2.Trials, 'Waves')
        ind0 = ~cellfun(@isempty, {ex0.Trials.Waves});
        ind2 = ~cellfun(@isempty, {ex2.Trials.Waves});
        exinfo.wdt = waveWidth( exinfo, vertcat(ex0.Trials(ind0).Waves),...
            vertcat(ex2.Trials(ind2).Waves), p_flag );
    else
        fprintf(['missing wave file' exinfo.fname '\n']);
    end
end




end


function  plotCOTC(exinfo)

% tuning curve
h = figure('Name', exinfo.figname);

ratemn = exinfo.fitparam.others.CO.val;
ratemn_drug = exinfo.fitparam_drug.others.CO.val;


[co, idx] = sort(exinfo.fitparam.others.CO.co);
ratemn = ratemn(idx);
[co_drug, idx_drug] = sort(exinfo.fitparam_drug.others.CO.co);
ratemn_drug = ratemn_drug(idx_drug);


%%%%%%%
subplot(1,2,1)
%baseline
plot(co, ratemn, 'o', 'MarkerEdgeColor', getCol(exinfo)); ho % raw
plot(exinfo.fitparam.others.CO.x, ...
    exinfo.fitparam.others.CO.y, '--', 'Color', getCol(exinfo)); ho % fit
%drug
plot(co_drug, ratemn_drug, 'o', 'MarkerEdgeColor', getCol(exinfo), ...
    'MarkerFaceColor', getCol(exinfo)); ho % raw
plot(exinfo.fitparam_drug.others.CO.x, ...
    exinfo.fitparam_drug.others.CO.y, '-', 'Color', getCol(exinfo)); ho % fit

legend('baseline', exinfo.drugname, 'Orientation', 'horizontal', ...
    'Location', 'Northoutside');
set(gca, 'XScale', 'log'); xlabel('contrast'); ylabel('spk/frame')
xlim([min(co) max(co)]);


subplot(1,2,2);
scatter(ratemn, ratemn_drug, 40, getCol(exinfo), 'MarkerFaceColor', getCol(exinfo)); ho

plot(ratemn, ...
    ratemn.*exinfo.fitparam.others.CO.gslope+exinfo.fitparam.others.CO.yoff, ...
    'Color', getCol(exinfo));

eqax; crossl; unity; hold on
xlabel('Baseline'); ylabel(exinfo.drugname);

title(sprintf('contrast: gain %1.2f, offset: %1.2f, r2 %1.3f', ...
    exinfo.fitparam.others.CO.gslope,...
    exinfo.fitparam.others.CO.yoff,...
    exinfo.fitparam.others.CO.r2reg));
savefig(h, exinfo.fig_phase);
close(h)
end



