function [beta1, beta0, xvar, bootstrp]  = type2reg(exinfo, p_flag)
% type2reg returns the offset and slope for the best fit linear regression
% that is subject to both errors (x and y).


% % only take those conditions that are present in both experiments
s1 = exinfo.ratepar;  %or/sf/co/... values used in the experiment
s2 = exinfo.ratepar_drug;

cont = exinfo.ratemn( ismember( s1, s2 ) );
drug = exinfo.ratemn_drug( ismember( s2, s1 ));


% calculate perpendicular distance minimizing fit
[beta0, beta1, xvar, bootstrp] = ...
    perpendicularfit(cont, drug, var(drug)/var(cont));


if any(drug<0)||any(cont<0)
drug2 = drug + abs(min([drug; cont])); cont2 = cont +abs(min([drug; cont]));
else
   drug2 = drug; cont2 = cont; 
end
drug2 = drug2./max(cont2); cont2 = cont2./max(cont2);
[beta0(2), beta1(2), xvar(2)] = ...
    perpendicularfit(cont2, drug2, var(drug2)/var(cont2));

% triangular fit
% [Lt, beta1, xvarXY] = trianglefit(x, y) ;
% beta0 = mean(y)-m*mean(x);

if p_flag
    
    h = figure('Name', exinfo.figname);
    subplot(1,2,1);
    scatter(cont, drug, 100, getCol(exinfo), 'filled',...
        'MarkerFaceAlpha', 0.5); hold on;
    if any(s1)>1000 && any(s2)>1000
        scatter(cont(end), drug(end), 100, 'g', 'filled',...
            'MarkerFaceAlpha', 0.5); hold on;
    end
    eqax;
    xlim_ = get(gca, 'xlim');
    plot(xlim_, xlim_.*beta1(1) + beta0(1), getCol(exinfo), 'LineWidth', 2); hold on;
    
    xlabel('baseline');
    ylabel(exinfo.drugname);
    
    title(sprintf( 'gain = %1.2f, add = %1.2f, \n r2=%1.2f', ...
        beta1(1), beta0(1), xvar(1)),'FontSize', 8);
    unity; axis square; box off;
    
    ylim_ = get(gca, 'YLim');
    hold on
    plot([mean(cont) mean(cont)], ylim_, 'Color', [0.5 0.5 0.5 0.5])
    plot(xlim_, [mean(drug) mean(drug)], 'Color', [0.5 0.5 0.5 0.5])
    
    %%% normalized
    subplot(1,2,2);
        scatter(cont2, drug2, 100, getCol(exinfo), 'filled',...
        'MarkerFaceAlpha', 0.5); hold on;
    if any(s1)>1000 && any(s2)>1000
        scatter(cont2(end), drug2(end), 100, 'g', 'filled',...
            'MarkerFaceAlpha', 0.5); hold on;
    end
    eqax;
    xlim_ = get(gca, 'xlim');
    plot(xlim_, xlim_.*beta1(2) + beta0(2), getCol(exinfo), 'LineWidth', 2); hold on;
    
    xlabel('baseline');
    ylabel(exinfo.drugname);
    
    title(sprintf( 'normalized gain = %1.2f, add = %1.2f, \n r2=%1.2f', ...
        beta1(2), beta0(2), xvar(2)),'FontSize', 8);
    unity; axis square; box off;
    
    ylim_ = get(gca, 'YLim');
    hold on
    plot([mean(cont) mean(cont)], ylim_, 'Color', [0.5 0.5 0.5 0.5])
    plot(xlim_, [mean(drug) mean(drug)], 'Color', [0.5 0.5 0.5 0.5])
    
    savefig(h, exinfo.fig_regl);
    close(h);
    
end

if any(isnan(xvar))
    xvar(isnan(xvar)) = 0;
    beta0(isnan(xvar)) = -inf; beta1(isnan(xvar)) = 0;
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
