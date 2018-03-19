function [ex, argout] = evalSingleRC( exinfo, fname, varargin )
% [ex, argout] = evalSingleRC( exinfo, fname )
% 
% 
% Counterpart to evalSingleDG but specified for flashed gratings,
% performing subspace reverse correlation analysis, latency estimation and
% fitting the tuning curve.
% 
% 
% 
% 
% 
% @CL

%% define variables and parse input
argout = [];
j=1; 
lat_flag = true;
p_flag = false;

while j<= length(varargin)
   switch varargin{j}
       case 'exinfo'
           j=j+1;
       case 'lat_flag'
           lat_flag = varargin{j+1};
           j=j+1;
       case 'plot'
           p_flag = true;
   end
   j=j+1;
end

%% load the ex file 

% ignore the lfp for now
ex = loadCluster( fname, 'loadlfp', false ); 

% add pupil info
[ex_orig] = znormex(ex, exinfo, 1);

for ps = 1:3
    ex = ex_orig;
        switch ps                
                case 1
                        suffix = '_1';
                        ex.Trials = ex.Trials(ex.pupil_split(1).idx);
                case 2
                        suffix = '_2';
                        ex.Trials = ex.Trials(ex.pupil_split(2).idx);
                case 3
                        suffix = '';
        end
        
    % Mean and Variance
    [ res, spkstats, fitparam ] = RCsubspace(ex, varargin{:} );

    if lat_flag && p_flag && ps==3
        % name and save the figure of the ML response latency estimate
        if contains(fname, '5HT') && contains(fname, 'NaCl')
            set(gcf, 'Name', [exinfo.figname '_' exinfo.drugname]);
            savefig(gcf, [exinfo.fig_sdfs(1:end-4) '_mlfit_drug.fig'])
            close(gcf)
        else
            set(gcf, 'Name', [exinfo.figname '_base']);
            savefig(gcf, [exinfo.fig_sdfs(1:end-4) '_mlfit_base.fig'])
            close(gcf)
        end
    end

    % res = getLat2D(res, exinfo);% only if there are two stimulus dimensions evaluated

    nrep = res.sdfs.n;
    if ~isempty(res.netSpikesPerFrameBlank)
        nrep = [nrep; res.sdfs.extras{1}.n];
    end

    %% tc height diff
    minspk = min(spkstats.mn);
    maxspk = max(spkstats.mn);
    tcdiff = (maxspk - minspk) ./ mean([maxspk, minspk]);


    %% Anova
    idx = find(spkstats.(exinfo.param1)<1000);
    p_anova = nan(length(idx),length(idx));
    for i = 1:length(idx)
        for j = i:length(idx)

    %         mu = spkstats.mn(j);
    %         ct = prctile(squeeze(spkstats.bootstrap(i,1,:)), [5 95]);
    %         p_anova(i, j) = mu>ct(1) && mu<ct(2);

            A = squeeze(spkstats.bootstrap(i,1,:));
            B = squeeze(spkstats.bootstrap(j,1,:));

            p_anova(i, j) = min([sum(A<B)/1000, sum(B>A)/1000]);
        end
    end

    % p_anova = max(max(p_anova))

    %% fixation precision & microsaccade
    fixspan = fixationSpan(ex, 0);
    amp = [];
    peakv = [];
    duration = [];
    angle = [];
    for i = 1:length(fixspan.trial)
        amp = [amp fixspan.trial(i).microsaccade.amp];
        peakv = [peakv fixspan.trial(i).microsaccade.peakv];
        duration = [duration fixspan.trial(i).microsaccade.duration];
        angle = [angle fixspan.trial(i).microsaccade.angle];
    end
    argout = [argout {['fixationspan' suffix],  [fixspan.accuracy, fixspan.theta, fixspan.SI]', ...
            ['variance_eye' suffix], [fixspan.variance]',...
            ['recentering_win' suffix],  [fixspan.recentering.window_offsetX; fixspan.recentering.window_offsetY]', ...
            ['recentering_eye' suffix], [fixspan.recentering.eye_offsetX; fixspan.recentering.eye_offsetY]', ...
            ['microsac_amplitude' suffix], amp, ...
            ['microsac_peakv' suffix], peakv, ...
            ['microsac_duration' suffix], duration, ...
            ['microsac_angle' suffix], angle, ...
            ['microsac_counts' suffix], arrayfun(@(x) x.microsaccade.counts, fixspan.trial)}];

    %% assign output arguments

    argout =  [argout {['fitparam' suffix], fitparam, ...
        ['rateMN' suffix], spkstats.mn, ['ratePAR' suffix], spkstats.or, ['rateSME' suffix], spkstats.sem, ...
        ['tcdiff' suffix], tcdiff, ['sdfs' suffix], res.sdfs,...
        ['times' suffix], res.times, ...
        ['resdur' suffix], res.dur, ['resvars' suffix], res.vars2,...
        ['lat' suffix], res.latFP, ['lat2Hmax' suffix], res.lat,...
        ['p_anova' suffix], p_anova, ['nrep' suffix], nrep}];
end

end



function res = getLat2D(res, exinfo)

psth = res.sdfs.psth;
n_or = size(psth, 1);

[~, co_idx] = sort(res.sdfs.y(1,:));
n_co = size(psth, 2);

h = findobj('Tag', 'latency');
if isempty(h)
    figure('Tag', 'latency');
    col = 'r';
else
    set(h, 'Tag', 'done latency');
    figure(h);
    col = 'b';
end

for or = 1:n_or
    kk = 1;
    for co = co_idx
        
       [lat2hmax_orco(or,kk), latfp_orxco(or, kk), pPoisson(or, kk)] =...
           CL_newLatency_helper(res.sdfs.s{or,co}, -20:0.1:160);
       
        subplot(n_or, n_co, kk+ (or-1)*n_co); hold on;
        plot(0:0.1:160, res.sdfs.s{or,co}(201:end), 'Color', col); hold on;
        
        plot([latfp_orxco(or,kk) latfp_orxco(or,kk)], get(gca, 'YLim'), 'k' )
        xlim([-10,160]);
        kk = kk+1;
    end
end

res.latfp_xco = latfp_orxco;
res.latfp_pPoisson_xco = pPoisson;
res.lat2hmax_xco = lat2hmax_orco;



h = findobj('Tag', 'done latency');
if ~isempty(h)
    set(h, 'Name', exinfo.figname);
    savefig(h, exinfo.fig_raster);
    close(h);
end
end



function [lat2hmax, latfp, pPoisson] = CL_newLatency_helper(sdf, times)

lat2hmax     = -1;
dur     = -1;
latfp   = -1;
pPoisson = 0;
noise   = mean(sdf(200:400)); 

if max(sdf)>mean(noise)*3.5
    
    sdf2 = sdf-noise;     % normalizer for baseline variability
    
    % first time of half max of response
    idx = find( sdf2 >= (max(sdf2)/2), 1, 'first');  
    lat2hmax = times(idx)/10;

    try
    [latfp, ~, pPoisson] = friedmanpriebe(round(sdf(200:end)), ...
        'minTheta', 250);
    catch
        latfp = -10;
    end
    latfp = latfp/10;
end

end

