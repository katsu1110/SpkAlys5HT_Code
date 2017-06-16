function addTitle( dat )
% compute the statistics for dat and write tham as title to the plot
% 
% dat is the result structure created in createUnitPlot.m
% 
% 
% @CL


i_5HT = dat.is5HT; % logical indexing idnicating 5HT drug experiments in the dat structure

if any(i_5HT) && any(~i_5HT)
    
    % 5HT experiments
    xdat1 = dat.x(i_5HT)';     ydat1 = dat.y(i_5HT)';
    
    mnx1 = mean(xdat1);         mny1 = mean(ydat1);
    medx1 = median(xdat1);         medy1 = median(ydat1);
    
    % NaCl experiments
    xdat2 = dat.x(~i_5HT)';    ydat2 = dat.y(~i_5HT)';
    
    mnx2 = mean(xdat2);         mny2 = mean(ydat2);
    medx2 = median(xdat2);         medy2 = median(ydat2);
    
    
    % transform relative metrices to logarithmic space to use the tests
    % that assume the normal distribution
    if contains(dat.xlab, 'gain') || contains(dat.xlab, 'nonparam')
        xdat = log(dat.x(i_5HT));  
    elseif contains(dat.ylab, 'gain') || contains(dat.ylab, 'nonparam')
        ydat = log(dat.y(i_5HT));
    end
    
    % pearson's and spearman's correlation coefficients between x and y
    [rho, p] = corr(xdat1, ydat1);
    [rhos, ps] = corr(xdat1, ydat1, 'type', 'Spearman');
    
    [rho2, p2] = corr(xdat2, ydat2);
    [rhos2, ps2] = corr(xdat2, ydat2, 'type', 'Spearman');
    
    % paired t-test, parametric and non-parametric, for x and y
    [~, pttest] = ttest(xdat1, ydat1);
    psign = signrank(xdat1, ydat1);
      
    % results transformed into a long string that is used as title
    title(sprintf(['5HT(n=%1.0f) \t r_{p}=%1.2f p=%1.3f \t r_{s}=%1.2f p=%1.3f   \\mu_x = %1.2f \\mu_y = %1.2f   med_x = %1.2f med_y = %1.2f \n '...
        'NaCl(n=%1.0f) \t r_{p}=%1.2f p=%1.3f \t r_{s}=%1.2f p=%1.3f \\mu_x = %1.2f \\mu_y = %1.2f   med_x = %1.2f med_y = %1.2f\n'...
        'paired ttest p=%1.3f, signrank test p=%1.3f' ],...
        length(xdat1),rho, p, rhos, ps, mnx1, mny1, medx1, medy1, ...
        length(xdat2), rho2, p2, rhos2, ps2,  mnx2, mny2, medx2, medy2, ...
        pttest, psign), 'FontSize', 8, 'FontWeight', 'bold');
    
    
elseif any(i_5HT) || any(~i_5HT)
    % if there is no control or no 5HT condition
    
    [rho, p1] = corr(dat.x', dat.y');
    [rhos, ps1] = corr(dat.x', dat.y', 'type', 'Spearman');
    
    [~, p2] = ttest(dat.x', dat.y');
    p_wil = ranksum(dat.x', dat.y');
    
    title(sprintf(['n = %1.0f \t r_{p} =%1.2f p=%1.3f r_{s} =%1.2f p=%1.3f \n'...
        '\t paired-t p=%1.3f \t wilk p=%1.3f'] ,...
        [length(dat.x), rho, p1, rhos, ps1, p2, p_wil]), 'FontSize', 8, 'FontWeight', 'bold');
end

end

