function [mn_,se_,cl,extra_mn_,extra_se_,vals1,vals2]= plotTC_copy(ex,varargin)

% function plotTC(ex)
%
% plots a tuning curve of the currently run experiment

% history
% 08/05/14  hn: wrote it
% 01/02/15  hn: included options to plot up to three clusters (cl 0 to 2)
% 03/03/15  hn: for adaptation do not extract mean response for adapter.

disp('in TC')

j=1;
plot_hold_on = 0;
lineStyle = '-';
while j<nargin
    str = varargin{j};
    switch str
        case {'linestyle','LineStyle','lineStyle','Linestyle'}
            lineStyle = varargin{j+1};
            j=j+1;
        case 'PlotHoldOn'
         plot_hold_on = 1;
    end
    j=j+1;
end
plot_hold_on
set(gcf,'position',[1182         416         735         570])

if (~isfield(ex.Trials(1),'oRate') && ~isfield(ex.Trials(end),'oRate')) && ...
    (~isfield(ex.Trials(1),'Rate') && ~isfield(ex.Trials(end),'Rate'))
    disp('no spikes found')
    return;
end

% don't plot responses for bar stimulus
if strcmpi(ex.stim.type,'bar')
    return
end

itr = (find((abs([ex.Trials.Reward])>0)));
if length(itr)<2
    disp('too few trials')
    return
end 

% if we have fixation trials take those
if ~isempty(itr)
    Trials = ex.Trials(itr);
else Trials = ex.Trials;  % for testing purposes when we don't have fixation trials
end

% do we have sorted Spikes?  If so, plot these
if isfield(Trials,'Spikes')
    for n = 1:length(Trials)
        spks = Trials(n).Spikes;
        sts = Trials(n).Start - Trials(n).TrialStart;
        if ex.stim.vals.adaptation
            sts = sts(find(sts>ex.stim.vals.adaptationDur));
        end
        if length(sts)>=1
            Trials(n).oRate.cl1 = length(find(spks>=sts(1) & spks<=sts(end)+0.1));
        else
            Trials(n).oRate.cl1 = [];
        end
        Trials(n).oRate.cl0 = [];
        Trials(n).oRate.cl2 = [];
    end
end

% for backward compatibility: (we now read out spikes for diff. clusters)
if ~isfield(Trials(1).oRate,'cl1') && ~isfield(Trials(end).oRate,'cl1')
    disp('in backward compatibility')
    for n = 1:length(Trials)
        Trials(n).oRate.cl1 = Trials(n).oRate;
        Trials(n).oRate.cl0 = [];
        Trials(n).oRate.cl2 = [];
    end
end

% first get trials for blank stimulus
extra_mn = [];
extra_se = [];
blTrials = [];
if ex.exp.include_blank
    bl_tr = find([Trials.st] ==0);
    st_tr = find([Trials.st] ==1);
    blTrials = Trials(bl_tr);
    Trials = Trials(st_tr);
end

% tuning parameters and values
par1=[]; par2=[];
if isfield(ex,'exp') && isfield(ex.exp,'e1')
    par1 = ex.exp.e1.type;
end
if isfield(ex,'exp') && isfield(ex.exp,'e2')
    par2 = ex.exp.e2.type;
end

vals1 = []; vals2 = [];
if ~isempty(par1)
    vals1 = ex.exp.e1.vals;
end
if ~isempty(par2)
    vals2 = ex.exp.e2.vals;
end

        
% get the sorted trials for the remaining stimuli 
legendstr = {};
tr = {};
for n1 = 1:length(vals1)
    if length(vals2)>0
        for n2 = 1:length(vals2)
            idx = find(eval(['[Trials.' par1 '] == vals1(n1) & [Trials.' ...
                par2 '] == vals2(n2)']));
            tr{n2,n1} = idx;
            N(n2,n1) = length(idx);
            legendstr{n2} = [par2 '=' num2str(vals2(n2))];
        end
    else
        idx = find(eval(['[Trials.' par1 '] == vals1(n1)'])); 
        tr{n1} = idx;
        N(n1) = length(idx);
    end
end

% get the the spikes for the individual clusters
cl = {['cl1'],['cl2'],['cl0']};
[mn_,se_,cl,extra_mn_,extra_se_] = getRate4Clusters(Trials,blTrials,tr,cl);


% now plot the data
switch length(mn_)
    case 0
        return
    case 1
        subplot_pos{1} = [0.1    0.25    0.87    0.7];
    case 2
        % 2 subplots:
        subplot_pos{1} = [.1    0.60    0.87    0.35];
        subplot_pos{2} = [0.1    0.25    0.87    0.34];
    case 3
        % 3 subplots:
        subplot_pos{1} = [0.1    0.7    0.87   0.25];
        subplot_pos{2} = [ 0.1    0.47    0.87    0.22];
        subplot_pos{3} = [ 0.1    0.25    0.87    0.21];
end

if plot_hold_on
    for n = 1:length(subplot_pos)
        y_pos(n) = subplot_pos{n}(2);
    end
    c = get(gcf,'children');
    for n=1:length(c)
        h = get(c(n));
        % we delete existing legends
        if isfield(h,'Location')
            delete(c(n));
        end
    end
    c = get(gcf,'children');
    for n = 1:length(c);
        pos = get(c(n),'position');
        y_pos(n) = pos(2);
    end
    for n = 1:length(subplot_pos)
        [min_d idx] = min(abs(y_pos-subplot_pos{n}(2)));
        subplot_axis(n) = c(idx);
    end
end
    
        
cols = colormap(hsv(size(mn_{1},1)));

pre_ylim = [];
for n = 1:length(mn_)
    
    
    if ~plot_hold_on;
        subplot('position',subplot_pos{n})
        hold off
    else
        set(gcf,'currentaxes',subplot_axis(n));
        pre_ylim = get(gca,'ylim');
        hold on;
    end
    mn = mn_{n};
    se = se_{n};
    for n1 = 1:size(mn,1)
        if n==1;
            ms = 8;
            lw = 2;
        end
        errorbar(vals1,mn(n1,:),se(n1,:),'s','lineStyle',lineStyle,'color',cols(n1,:),...
            'markersize',ms,'markerfacecolor',cols(n1,:),'linewidth',lw);
        hold on;
        for n2 = 1:size(mn,2);
            text(1.03*vals1(n2),1.05*mn(n1,n2),num2str(N(n1,n2)),'fontsize',18);
        end
    end
    extra_mn = extra_mn_{n};
    extra_se = extra_se_{n};
    if ~isempty(extra_mn) && ~isempty(mn) && length(vals1)>=2 
        legendstr{length(legendstr)+1} = 'blank';
        errorbar(2*vals1(end)-vals1(end-1),extra_mn,extra_se,'ko','markersize',ms);
        hold on;
    end
    % format plot
    if sum(isnan(mn))>0
        ylim = [0 100];
    else
        ylim = [0 max(max(mn))+max([0 0.2*max(max(mn))+1])];
        if ~isempty(pre_ylim)
            ylim = [ylim(1) max([pre_ylim(2) ylim(2)])];
        end
    end
    set(gca,'ylim',ylim, 'fontsize',15,'fontweight','bold')
    if strcmpi(par1,'sf') || strcmpi(par1,'tf') || strcmpi(par1,'sz')|| strcmpi(par1,'co')
        set(gca,'xscale','log')
    end
    ylabel ('spike rate (ips)');
    if n == length(mn_)
        xlabel(par1)
    else set(gca,'xticklabel','');
    end
    if ex.stim.vals.adaptation
        plot(ones(1,2)*ex.stim.vals.adaptationOr,get(gca,'ylim'),'--k','linewidth',1)
    end
    
    if ~isempty(legendstr) 
        legend(legendstr,'Location','NorthEastOutside')
    end
    ntr = length(find((abs([Trials.Reward])>0)));
    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    if n == 1;
       title(['# of trials: ' num2str(ntr) ])
    end
    xp = xlim(2) - 0.1*diff(xlim);
    yp = ylim(1)+ 0.1*diff(ylim);


    t(n) = text(xp,yp,['cluster ' cl{n}(3)]);
end


% -------------------------------------------------------------------
function [mn,se,cl_new,extra_mn,extra_se] = getRate4Clusters(Trials,blTrials,tr,cl);
% check for which clusters we have spikes
for nc = 1:length(cl)
    if isfield(Trials(1).oRate,cl{nc})
        for n = 1:length(Trials)
                eval([cl{nc} '(n)= mean(Trials(n).oRate.' cl{nc} ');']);
        end
    else eval([cl{nc} '=NaN;'])
    end
end       
    
for nc = 1:length(cl)
    eval([cl{nc} '=' cl{nc} '(~isnan(' cl{nc} '));']);
    if isempty(eval(cl{nc}))
        cl{nc} = {};
    elseif eval(['sum(' cl{nc} ')==0'])
        cl{nc} = {};
    end
end
% get means and SEs for clusters that have spikes
cnt = 0;
cl_new = {};
mn = {};
se = {};
extra_mn = {};
extra_se = {};
for nc = 1:length(cl)
    if ~isempty(cl{nc})
        cnt = cnt+1;
        cl_new{cnt} = cl{nc};
        iRate  = [];
        for n1 = 1:size(tr,2)
            for n2 = 1:size(tr,1)
                itr = tr{n2,n1};
                for n3 = 1:length(itr);
                    iRate(n3) = eval(['Trials(itr(n3)).oRate.' cl{nc} ';']);
                end
                mn{cnt}(n2,n1) = mean(iRate); 
                se{cnt}(n2,n1) = std(iRate)/sqrt(length(itr)); 
            end
        end
        iRate  = [];
        extra_mn{cnt} = NaN;
        extra_mn{cnt} = NaN;
        for n1 = 1:length(blTrials)
            iRate(n1) = eval(['blTrials(n1).oRate.' cl{nc} ';']);
        end
        extra_mn{cnt}= mean(iRate); 
        extra_se{cnt}= std(iRate)/sqrt(length(iRate)); 

    end
end

