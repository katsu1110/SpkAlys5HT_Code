function [ ppsz_b, ppsz_d, fig_h] = ppsz_mu( trials_base, trials_drug, drugname, p_flag)
% drug is red, base is blue

fig_h = -1;

% loop over the four windows
for w = unique([trials_base.window])
    
    traj.all{w} = [vertcat(trials_base([trials_base.window] == w).puptraj); ...
        vertcat(trials_drug([trials_drug.window] == w).puptraj)];
    
    traj.mn(w,:) = nanmean( traj.all{w} ) ;
    traj.sd(w,:) = nanstd( traj.all{w} ) ;
    traj.n_base{w} = 1 : sum( [trials_base.window] == w );
    traj.n_drug{w} = traj.n_base{w}(end) + 1 : ...
        traj.n_base{w}(end) + sum( [trials_drug.window] == w );
    
    traj.alldev{w} = bsxfun(@minus, traj.all{w}, traj.mn(w,:));
    
    ppsz_mn = nanmean( bsxfun(@minus, traj.all{w}, traj.mn(w,:)) , 2 );
    ppsz_base{w} = sort(ppsz_mn( traj.n_base{w} )) ;
    ppsz_drug{w} = sort(ppsz_mn( traj.n_drug{w} )) ;
    
    clear ppsz_mn;
    
end



if p_flag
    fig_h1 = plotpstraj(traj, drugname);
    fig_h2 = figure('Name', 'Before and After Exclusion', 'Position',  [554 546 1064 420]);
end

% equlize distributions
for w=unique([trials_drug.window])
    
    [h,p] = ttest2( ppsz_base{w}, ppsz_drug{w} );
    
    if p_flag
        tit = sprintf(['w%1.0f - ttest %1.2f - base(blue,n %1.0f, n_{excl} %1.0f)\n ' drugname '(red,n %1.0f, n_{excl} %1.0f)'], ...
            w, p, length(ppsz_base{w}),  length(traj.n_base{w}) - length(ppsz_base{w}), ...
            length(ppsz_drug{w}), length(traj.n_drug{w}) - length(ppsz_drug{w}));
        ppszhist(ppsz_base{w}, ppsz_drug{w}, fig_h2, w, tit);
        
    end
    
    % decreas population until they are not significantly different no more
    while h && ~isempty(ppsz_base{w}) && ~isempty(ppsz_drug{w})
        
        if mean(ppsz_base{w}) > mean(ppsz_drug{w})
            ppsz_base{w} = ppsz_base{w}( 1:end-1 );
            ppsz_drug{w} = ppsz_drug{w}( 2:end );
        else
            ppsz_drug{w} = ppsz_drug{w}( 1:end-1 );
            ppsz_base{w} = ppsz_base{w}( 2:end );
        end
        
        [h,p] = ttest2( ppsz_base{w}, ppsz_drug{w} );

        if p_flag
            tit = sprintf(['w%1.0f - ttest %1.2f - base(blue,n %1.0f, n_{excl} %1.0f)\n ' drugname '(red,n %1.0f, n_{excl} %1.0f)'], ...
                w, p, length(ppsz_base{w}),  length(traj.n_base{w}) - length(ppsz_base{w}), ...
                length(ppsz_drug{w}), length(traj.n_drug{w}) - length(ppsz_drug{w}));
            ppszhist(ppsz_base{w}, ppsz_drug{w}, fig_h2, w+4, tit)
        end
        
    end
    ppsz_b.(strcat('w' , num2str(w))) =mean(ppsz_base{w});
    ppsz_d.(strcat('w' , num2str(w))) =mean(ppsz_drug{w});
    
    
end

if p_flag
    fig_h(1) = fig_h1;
    fig_h(2) = fig_h2;
end

end


% plots histograms in a certain subplot w
function ppszhist(p1, p2, h, w, tit)

h(w) = subplot(2,4,w);
cla
cax = h(w);

% base
hist( cax, p2 );
h1 = findobj(cax,'Type','patch');
set(h1,'FaceColor', 'r', ...
    'EdgeColor', 'w', ...
    'facealpha', 0.7);
hold on;

% drug
hist( cax, p1 );
h1 = findobj(cax,'Type','patch');
set(h1, 'facealpha', 0.7);

title(tit, 'FontSize', 8, 'FontWeight', 'bold');

max_y = max(get(cax, 'Ylim'));
ylim([0 max_y + 1]);
max_y = max_y+1;

plot(cax, mean(p1), max_y, 'vk', 'MarkerFaceColor', 'k');
text(mean(p1), max_y-max_y/10, sprintf('%1.2f', mean(p1)), 'Parent', cax,'FontWeight', 'bold', 'FontSize', 8);

plot(cax ,mean(p2), max_y, 'vk', 'MarkerFaceColor', 'k');
text(mean(p2), max_y-max_y/5, sprintf('%1.2f', mean(p2)), 'Parent', cax, 'FontWeight', 'bold', 'FontSize', 8 );


end


% plot deviation from mean
function h = plotpstraj(traj, drugname)

h = figure('Name', 'Before Exclusion');

if size(traj.alldev{1},2) == 1000
    off = [0 2];
    inc = 2/999;
else 
    off = [0 0.5 1 1.5 2];
    inc = 0.5/224;
end
plot(0,0,'b'); hold on;
plot(0,0,'r');hold on;


for w = 1:length(traj.all)
    x = off(w) : inc : off(w+1);
    
    n1 = length(traj.n_base{w});
    n2 = length(traj.n_drug{w});
    
    plot(repmat(x,n1,1)' ,traj.alldev{w}(traj.n_base{w},:)', 'b'); hold on; % base
    plot(repmat(x,n2,1)',traj.alldev{w}(traj.n_drug{w},:)', 'r'); hold on; % drug
end

xlim([0 2]);
ylim_ = get(gca, 'ylim');


plot([0.5 0.5], ylim_, 'k--');
plot([1 1], ylim_, 'k--');
plot([1.5 1.5], ylim_, 'k--');
plot([0 2], [0 0], 'k', 'LineWidth', 2);

xlabel('approx time [ms]');
ylabel('ppsz deviation from mean');
legend('base', drugname);
set(gca, 'XTick', off);

end