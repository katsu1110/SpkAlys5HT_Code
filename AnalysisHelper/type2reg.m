function [b1, b0, xvar, bootstrp]  = type2reg(exinfo, p_flag)
% type2reg returns the offset and slope for the best fit linear regression
% that is subject to both errors (x and y).


% minimum response for restricted regression fit (revision)
spkmin = 0;


% % only take those conditions that are present in both experiments
s1 = exinfo.ratepar;  %or/sf/co/... values used in the experiment
s2 = exinfo.ratepar_drug;

cont = exinfo.ratemn( ismember( s1, s2 ) );
drug = exinfo.ratemn_drug( ismember( s2, s1 ));


%-------------------------- calculate perpendicular distance minimizing fit
[b0, b1, xvar] = ...
    perpendicularfit(cont, drug, var(drug)/var(cont));

%-------------------------------------- repeat for normalized tuning curves
if any(drug<0)||any(cont<0)
    drug_norm = drug + abs(min([drug; cont])); cont_norm = cont +abs(min([drug; cont]));
else
    drug_norm = drug; cont_norm = cont;
end
drug_norm = drug_norm./max(cont_norm); cont_norm = cont_norm./max(cont_norm);
[b0(2), b1(2), xvar(2)] = ...
    perpendicularfit(cont_norm, drug_norm, var(drug_norm)/var(cont_norm));



%-------------------------- bootstrap the parameters using the resampled tc
if ~exinfo.isRC
    idx_nonullresp = cont> spkmin | drug> spkmin; % everything below 10% of the smallest response
    
    cont_wo_null = cont_norm(idx_nonullresp);
    drug_wo_null = drug_norm(idx_nonullresp);
    
    [b0(3), b1(3), xvar(3)] = perpendicularfit(cont_wo_null, drug_wo_null, ...
        var(drug_wo_null)/var(cont_wo_null));
    
    [bootstrp.b0, bootstrp.b1, bootstrp.xvar] = ...
        bootRegII( exinfo.rate_resmpl(ismember( s1, s2 )), ...
        exinfo.rate_resmpl_drug(ismember( s2, s1 )), spkmin );
    
    % exclude data with less than 4 data points via the regression fit
    if length(cont_wo_null)<4 || length(drug_wo_null)<4
        xvar(3) = 0;
    end
else
    b0(3) = 0;
    b1(3) = 1;
    xvar(3) =0;
    bootstrp = [];
end

%-------------------------------------------------- alternative computation
% triangular fit
% [Lt, beta1, xvarXY] = trianglefit(x, y) ;
% beta0 = mean(y)-m*mean(x);

%--------------------------------------------------------- plot the results
if p_flag
    
    %%% raw tuning curves
    h = figure('Name', exinfo.figname);
    subplot(2,2,1);
    scatter(cont, drug, 100, getCol(exinfo), 'filled',...
        'MarkerFaceAlpha', 0.5); hold on;
    if any(s1)>1000 && any(s2)>1000
        scatter(cont(end), drug(end), 100, 'g', 'filled',...
            'MarkerFaceAlpha', 0.5); hold on;
    end
    eqax;
    xlim_ = get(gca, 'xlim');
    plot(xlim_, xlim_.*b1(1) + b0(1), getCol(exinfo), 'LineWidth', 2); hold on;
    
    xlabel('baseline');
    ylabel(exinfo.drugname);
    
    title(sprintf( 'gain = %1.2f, add = %1.2f, \n r2=%1.2f, n=%1d', ...
        b1(1), b0(1), xvar(1), length(cont)),'FontSize', 8);
    unity; axis square; box off;
    
    ylim_ = get(gca, 'YLim');
    hold on
    plot([mean(cont) mean(cont)], ylim_, 'Color', [0.5 0.5 0.5 0.5])
    plot(xlim_, [mean(drug) mean(drug)], 'Color', [0.5 0.5 0.5 0.5])
    
    %%% normalized tuning curves
    subplot(2,2,2);
        scatter(cont_norm, drug_norm, 100, getCol(exinfo), 'filled',...
        'MarkerFaceAlpha', 0.5); hold on;
    if any(s1)>1000 && any(s2)>1000
        scatter(cont_norm(end), drug_norm(end), 100, 'g', 'filled',...
            'MarkerFaceAlpha', 0.5); hold on;
    end
    eqax;
    xlim_ = get(gca, 'xlim');
    plot(xlim_, xlim_.*b1(2) + b0(2), getCol(exinfo), 'LineWidth', 2); hold on;
    
    xlabel('baseline');
    ylabel(exinfo.drugname);
    
    title(sprintf( 'gain = %1.2f, norm. offset  = %1.2f, \n r2=%1.2f, n=%1d', ...
        b1(2), b0(2), xvar(2), length(cont_norm)),'FontSize', 8);
    unity; axis square; box off;    
    
    %%% normalized tuning curves without null responses
    subplot(2,2,3);
        scatter(cont_wo_null, drug_wo_null, 100, getCol(exinfo), 'filled',...
        'MarkerFaceAlpha', 0.5); hold on;
    eqax;
    xlim_ = get(gca, 'xlim');
    xlim_ = [0 xlim_(2)];
    set(gca, 'XLim', xlim_, 'YLim', xlim_); 
    plot(xlim_, xlim_.*b1(3) + b0(3), getCol(exinfo), 'LineWidth', 2); hold on;
    
    xlabel('baseline');
    ylabel(exinfo.drugname);
    
    title(sprintf( 'w/o zeros \n gain = %1.2f, norm. offset = %1.2f, \n r2=%1.2f, n=%1d', ...
        b1(3), b0(3), xvar(3), length(cont_wo_null)),'FontSize', 8);
    unity; axis square; box off;

    
    subplot(4,2,6);
    histogram(bootstrp.b0); ylabel('bootstrap'); xlabel('offset (5 and 95%CI)');
    hold on;
    plot(ones(2,1)*prctile(bootstrp.b0, 5), [0 200], 'k:')
    plot(ones(2,1)*prctile(bootstrp.b0, 95), [0 200], 'k:')
    crossl;
    
    subplot(4,2,8);
    histogram(log(abs(bootstrp.b1))); ylabel('bootstrap'); xlabel('slope (5 and 95%CI)');
    hold on;
    plot(ones(2,1)*prctile(log(abs(bootstrp.b1)), 5), [0 200], 'k:')
    plot(ones(2,1)*prctile(log(abs(bootstrp.b1)), 95), [0 200], 'k:')
    crossl;
    
    % save results
    savefig(h, exinfo.fig_regl);
    close(h);
    
end

if any(isnan(xvar))
    xvar(isnan(xvar)) = 0;
    b0(isnan(xvar)) = -inf; b1(isnan(xvar)) = 0;
end

end




function [boot_b0, boot_b1, boot_xvar] = bootRegII(resamples_B, resamples_D, spkmin)
% use the resampled stimulus responses to bootstrap the regression fit.
% factor out the null responses to avoid the iceberg effect confounding the
% resulting parameters.


%loop through the stimuli and average each resampled data set
%the result matrices are bootstrapped tuning curves grouped in columns
for j = 1:length(resamples_B)
    tc_res_b(j,:) = mean(resamples_B{j}, 1);
    tc_res_d(j,:) = mean(resamples_D{j}, 1);
end

% initalize variables
boot_b0 = nan(size(tc_res_b, 2), 1);
boot_b1 = boot_b0;
boot_xvar = boot_b0;

% loop throught the resamples
for k = 1:size(tc_res_b, 2)
    
    base = tc_res_b(:, k); 
    drug = tc_res_d(:, k); 
    

    
    % remove all zero-responses
    idx_nullresp = drug< spkmin | base< spkmin;
    
    drug = drug(~idx_nullresp);
    base = base(~idx_nullresp);

    % normalize
    drug = drug./max(base);
    base = base./max(base);
    
    if length(drug)>3 && length(base)>3
    
        [boot_b0(k), boot_b1(k), boot_xvar(k)] = ...
            perpendicularfit(base, drug, var(drug)/var(base));
    else
        continue;
    end
end

end


function [Lt, m, expvar] = trianglefit(x, y)
% triangular area calculation after Barker et al 1980

% center data around 0 ( easier calculation )
x2 = x - mean(x);
y2 = y - mean(y);

% estimated b1 via triangle minimization
m = sign(sum( y2.*x2 )) * sqrt( sum( y2.^2 ) / sum( x2.^2 ) );


ssx = sum( x2.^2 );
ssy = sum( y2.^2 );

% area of triangles
Lt = - sign(sum( x2.* y2 )) * sum( x2.* y2 ) + sqrt(ssx.*ssy);

expvar = 1-(Lt / sqrt(ssx*ssy));

end
