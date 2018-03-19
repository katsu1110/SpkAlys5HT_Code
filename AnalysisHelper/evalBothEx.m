function exinfo = evalBothEx(ex0, ex2, exinfo, p_flag, varargin)
%evalBothEx evaluates operations on both base (ex0) and drug (ex2) files
% and performed the comparative analyses
%
% @CL

% ==================== preferred/unpreferred stimuli occuring in both files
idx = ismember(exinfo.ratepar, exinfo.ratepar_drug) &  exinfo.ratepar<1000;
parb = exinfo.ratemn;

parb(~idx, :) = 0;
[~, exinfo.pfi] = max( parb );
exinfo.pfi_drug = find(exinfo.ratepar_drug == exinfo.ratepar( exinfo.pfi ));

parb(~idx) = 10^4;
[~, exinfo.upfi] = min( parb );
exinfo.upfi_drug = find(exinfo.ratepar_drug == exinfo.ratepar( exinfo.upfi ));


% =============================================================== bootstrap
exinfo = bootstrap_exinfo( exinfo );

% ====================================================== type-II regression
if any(strcmp(varargin, 't2reg')) || any(strcmp(varargin, 'all'))
    [gslope, yoff, r2, bootstrp] = type2reg(exinfo, p_flag);
    
    exinfo.gslope = gslope(1);
    exinfo.yoff = yoff(1);
    exinfo.r2reg = r2(1);
    
    exinfo.yoff_rel = yoff(2);
    %exinfo.reg_bootstrp = bootstrp;
end

% ====================================================== fitting parameters
if any(strcmp(varargin, 'tcfit')) || any(strcmp(varargin, 'all'))
    if isfield(exinfo.fitparam, 'others') && isfield(exinfo.fitparam.others, 'OR')
        % RC experiments with varying or and co only
        
        % compute the type-II regression for orientation tuning curves for each contrast condition ... 
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
        
        
        % ... and for contrast tuning curves averaged across all orientations
        cont = exinfo.fitparam.others.CO.val;
        drug = exinfo.fitparam_drug.others.CO.val;
        
        [beta0, beta1, xvar] = ...
            perpendicularfit(cont, drug, var(drug)/var(cont));
        
        exinfo.fitparam.others.CO.yoff = beta0;
        exinfo.fitparam.others.CO.gslope = beta1;
        exinfo.fitparam.others.CO.r2reg = xvar;
        
        plotCOTC(exinfo);
        
        
    elseif strcmp(exinfo.param1, 'sf')      %%% spatial frequency
        
        % find the best fit - 
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
    
    
    %%% pfreferred orientation fit evaluation
    % 
    % in case the preferred orientations (mu) are more than 150 degrees
    % apart they are at different extremes of the circular values and their difference could be misleading.
    % in this case, I allow mu to be in the range of [0 210].
    if strcmp(exinfo.param1, 'or')
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


% ======================================================= waveform analysis
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

% ======================================================= pupil GLM
if any(strcmp(varargin, 'pupil')) || any(strcmp(varargin, 'all'))
    len_tr0 = length(ex0.Trials);
    len_tr2 = length(ex2.Trials);    
    
    % structure initialization
    exinfo.psglm.tc = [];
    exinfo.psglm.inter_table = nan(2,2);
    exinfo.psglm.b = [];
    exinfo.psglm.perf = nan(1,7);
    
    % store pupil size time-course
    % -- drug -- label_tr -- TC ---
    nmin = 1000;
    for i = 1:len_tr0
        if length(ex0.Trials(i).pupil_raw) < nmin
            nmin = length(ex0.Trials(i).pupil_raw);
        end
    end
    for i = 1:len_tr2
        if length(ex2.Trials(i).pupil_raw) < nmin
            nmin = length(ex2.Trials(i).pupil_raw);
        end
    end
    
    psmat0 = zeros(len_tr0 , 2 + nmin);
    for i = 1:len_tr0
        psmat0(i,2) = ex0.Trials(i).n_stm;
        psmat0(i,3:end) = ex0.Trials(i).pupil_raw(1:nmin);
    end
    psmat2 = ones(len_tr2 , 2 + nmin);
    for i = 1:len_tr2
        psmat2(i,2) = ex2.Trials(i).n_stm;
        psmat2(i,3:end) = ex2.Trials(i).pupil_raw(1:nmin);
    end
    
    exinfo.psglm.tc = [psmat0; psmat2];
        
    % mean firing rate in the baseline condition    
    if exinfo.isRC==1
        bfr = mean([ex0.Trials.spkRate]);
    end
    
    % use only the last stimulus in the DG experiments
    if ~exinfo.isRC
        ex0.Trials = ex0.Trials(psmat0(:,2)==4);
        ex2.Trials = ex2.Trials(psmat2(:,2)==4);
        len_tr0 = sum(psmat0(:,2)==4);
        last_tr2 = sum(psmat2(:,2)==4);
    end
    
    % make matrices for GLM
    mat0 = zeros(len_tr0, 5);
    mat2 = ones(len_tr2, 5);
    
    % baseline
    for i = 1:len_tr0        
        
        % mean firing rate in this stimulus
        try
            if exinfo.isRC==1
                mat0(i,2) = bfr;
            else
                mat0(i,2) = exinfo.ratemn(exinfo.ratepar(1:end-1)==ex0.Trials(i).(exinfo.param1));
            end
        catch
            continue
        end
        
        % firing rate in this trial
        mat0(i,1) = ex0.Trials(i).spkRate;
        
        % drug condition
        mat0(i,3) = 0;
        
        % pupil size
        if exinfo.isRC
            mat0(i,4) = ex0.Trials(i).pupil_val;
        else
            mat0(i,4) = mean(ex0.Trials(i).pupil_raw);
        end
        
        % interaction ... mat0(i,5) = 0
    end
    
    % drug 
    for i = 1:len_tr2        
        
        % mean firing rate in this stimulus
        try
            if exinfo.isRC==1
                mat2(i,2) = bfr;
            else
                mat2(i,2) = exinfo.ratemn(exinfo.ratepar_drug(1:end-1)==ex2.Trials(i).(exinfo.param1));
            end
        catch
            continue
        end
        
        % firing rate in this trial
        mat2(i,1) = ex2.Trials(i).spkRate;
        
        % drug condition
        mat2(i,3) = 1;
        
        % pupil size
        if exinfo.isRC
            mat2(i,4) = ex2.Trials(i).pupil_val;
        else
            mat2(i,4) = mean(ex2.Trials(i).pupil_raw);
        end
        
        % interaction
        mat2(i,5) = mat2(i,4);
    end
        
    % data matrix
    mat = [mat0; mat2];
    mat = mat(~any(isnan(mat),2),:); 
    mat(:,3:end) = zscore(mat(:,3:end));

    % stimulus only    
    exinfo.psglm.perf(1) = 1 - (var(abs(mat(:,2) - mat(:,1)))/var(mat(:,1)));
    
    % interaction table
    %%%%%%%%%%%%%%%%%%%%%%%%
    % FR %% 5HT %%% base %%%
    % S-ps      %%%      %%%
    % L-ps      %%%      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    med0 = median(mat0(:,4));
    med2 = median(mat2(:,4));
    exinfo.psglm.inter_table(1,1) = mean(mat2(mat2(:,4) < med2,1));
    exinfo.psglm.inter_table(2,1) = mean(mat2(mat2(:,4) > med2,1));
    exinfo.psglm.inter_table(1,2) = mean(mat0(mat0(:,4) < med0,1));
    exinfo.psglm.inter_table(2,2) = mean(mat0(mat0(:,4) > med0,1));
    
    % normalize the table by the grand mean
    exinfo.psglm.inter_table = exinfo.psglm.inter_table/...
        mean(mean(exinfo.psglm.inter_table));
    
    % model fit (for gain change)
    for i = 1:3
        % fitting
        x_init = zeros(1, i);
        c = @(x) cost(x, mat, i);    
        options = optimset('MaxFunEvals',3000);
        [b1] = fminsearch(c, x_init, options);
        b = nan(1, 6);
        b(1:length(b1)) = b1;
        exinfo.psglm.b = [exinfo.psglm.b; b];
        
        % variance explained
        mout = mymodel(b1, mat, i);
        exinfo.psglm.perf(i+1) = 1 - (var(abs(mout - mat(:,1)))/var(mat(:,1)));
    end

    % model fit (for additive change)
    resid = abs(mat(:,1) - mout);
    for i = 1:3
        % fitting
        mdl = fitglm(mat(:,3:2+i), resid, ...
            'Distribution','normal', 'link', 'identity', 'Intercept', false);
        b2 = mdl.Coefficients.Estimate';
        b = nan(1, 6);
        b(4:3+length(b2)) = b2;
        exinfo.psglm.b = [exinfo.psglm.b; b];
        
        % variance explained
        mout_add = predict(mdl, mat(:,3:2+i));
        exinfo.psglm.perf(i+4) = 1 - (var(abs(mout + mout_add - mat(:,1)))/var(mat(:,1)));
    end
    
end

end

%% subfunctions

function  plotCOTC(exinfo)
% plots the contrast tuning curve from RC experiments
% only occurs if both or and co varied in the experiment
% 

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

% model fit for GLM with pupil size =====================
function mout = mymodel(x, mat, incl)
% fr(t) = fr(s)*exp(sum of weighted predictors)
switch incl
    case 1
        mout = mat(:,2).*exp(x(1)*mat(:,3));
    case 2
        mout = mat(:,2).*exp(x(1)*mat(:,3) + x(2)*mat(:,4));
    case 3
        mout = mat(:,2).*exp(x(1)*mat(:,3) + x(2)*mat(:,4) + x(3)*mat(:,5));
end
end

function f = cost(x, mat, incl)
mout = mymodel(x, mat, incl);
f = 0;
for i = 1:size(mat,1)
    if mout(i) > 0
        f = f - (mat(i,1)*log(mout(i)) - mout(i));
    end
%     f = f + sqrt((mat(i,1) - mout(i))^2);
end
end





