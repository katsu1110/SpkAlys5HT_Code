function dat = evalMU(fctX, fctY, exinfo)

if contains(fctX, 'time')
    dat = assignFctXY(fctX, fctY, exinfo);
elseif  contains(fctY, 'time')
    dat = assignFctXY(fctY, fctX, exinfo);
else
    [dat.x, ~, dat.xlab, ~] = assignFct(fctX, exinfo);
    % error and info are only in case the average varaince is plotted as a
    % function of binned mean firing rates
    [dat.y, dat.err, dat.ylab, dat.info] = assignFct(fctY, exinfo);
end

if size(dat.x) == size(dat.y)
    dat.info = sprintf([dat.info ' results %3f %3f'], corr(dat.x, dat.y));
end



dat.exinfo = exinfo(~isnan(dat.x));
if size(dat.y,1)>1
    dat.y  = dat.y(~isnan(dat.x),:);
else
    dat.y  = dat.y(~isnan(dat.x));
end
dat.x  = dat.x(~isnan(dat.x));

end


function [val, err, lab, addinfo] = assignFct(fctname, exinfo)
% 


lab = fctname;  % axis label
addinfo = [];   % additional information
err = [];       % error for errorbars

ff_min_nrep = 8;%4; % minimum #of stimulus repeats for inclusion to fano factor analysis
ff_min_spks = 0;%2.5; % minimum #of spikes in each stimulus condition for inclusion to fano factor analysis


alpha_p_recovery = 0.05;

% switch between potential function names and call the retrieve the
% information from the linked field in exinfo
switch fctname
    
    
    
    case 'meanfr diff base-drug (w recovery)'
        val = assignFct('meanfr base (w recovery)' , exinfo) ./ ...
            assignFct('meanfr drug (w recovery)' , exinfo);
        lab = [fctname ' (base - drug)'];
    
    case 'meanfr base (w recovery)'
        for i = 1:length(exinfo)
            if exinfo(i).recov_p>alpha_p_recovery && ~isempty(exinfo(i).recov_fname)
                
                stim_base = exinfo(i).ratepar;
                stim_drug = exinfo(i).ratepar_drug;
                i_both = ismember(stim_base, stim_drug);

                val(i) = mean(exinfo(i).ratemn(i_both));
            else
                val(i) = nan;
            end
        end
         
        
    case 'meanfr recovery'
        for i = 1:length(exinfo)
            if exinfo(i).recov_p>alpha_p_recovery && ~isempty(exinfo(i).recov_fname)
                
                stim_base = exinfo(i).ratepar;
                stim_drug = exinfo(i).ratepar_drug;
                i_both = ismember(stim_base, stim_drug);
                
        
                ex = loadCluster(exinfo(i).recov_fname{end}, ...
                    'ocul', exinfo(i).ocul, ...
                    'loadlfp', false);
                [~, spkrate] = znormex(ex, exinfo(i));
                
                
                stim_rec = [spkrate.(exinfo(i).param1)];
                i_all = ismember(stim_rec, intersect(stim_base, stim_drug));
                
                val(i) = mean([spkrate(i_all).mn]);
            else
                val(i) = nan;
            end
        end
        
        
    case 'meanfr drug (w recovery)'
        
        for i = 1:length(exinfo)
            
            if exinfo(i).recov_p>alpha_p_recovery && ~isempty(exinfo(i).recov_fname)
                stim_base = exinfo(i).ratepar;
                stim_drug = exinfo(i).ratepar_drug;
                i_both = ismember(stim_drug, stim_base);
                val(i) = mean( exinfo(i).ratemn_drug(i_both) );
                
            else
                val(i) = nan;
            end
        end
    
    case 'p modulation ttest'
        
        val = cellfun(@(x) x(1), {exinfo.pmodulation});
        
    case 'p modulation wilcoxon'
        
        val = cellfun(@(x) x(2), {exinfo.pmodulation});
        
    case 'p recovery ttest'
        
        val = cellfun(@(x) x(1), {exinfo.ret2base});
    
    case 'p recovery wilcoxon'
        
        val = cellfun(@(x) x(2), {exinfo.ret2base});
        
    case 'circular variance'
        val = [exinfo.circvar];
        
    case 'direction selectivity'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).dirsel;
        end
    
    case 'direction selectivity drug'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).dirsel_drug;
        end
            
    case 'direction selectivity diff'
        val = assignFct('direction selectivity drug' , exinfo) ./ ...
            assignFct('direction selectivity' , exinfo);
        lab = [fctname ' (drug/base)'];
        
    case 'fixation accuracy base'
        
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fixspan_base.accuracy;
        end
        
    case 'fixation accuracy drug'
        
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fixspan_drug.accuracy;
        end
        
    case 'fixation accuracy diff'
        
        val = assignFct('fixation accuracy base' , exinfo) - ...
            assignFct('fixation accuracy drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'delta n to 1st 5HT app'
        % number of drug experiments between 1st application of drug and this
        val = [exinfo.dn_id];
        val(val==inf) = -1;
        
    case 'delta t to 1st 5HT app'
        val = [exinfo.dt_cum_id];
        val(val==inf) = -1;
    
        lab = [fctname ' (sec)'];

    case 't of 5HT app before'
        val = [exinfo.dt_cum_id];
        val(val==inf) = -1;
    
        lab = [fctname ' (sec)'];

        
    case 'isolation quality base'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).spkqual_base;
        end
        
    case 'isolation quality drug'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).spkqual_drug;
        end
        
    case 'isolation quality max'
        for kk = 1:length(exinfo)
            val(kk) = max([exinfo(kk).spkqual_base, exinfo(kk).spkqual_drug]);
        end
        
        
    case 'eccentricity'
        val = [exinfo.ecc];
        
    case 'RF ecc corr'
        for kk = 1:length(exinfo)
            
            if isempty(exinfo(kk).RFw_corr)
                val(kk) = -100;
            else
                val(kk) = exinfo(kk).RFw_corr;
            end
        end
        
    case 'RF width x'
        for kk = 1:length(exinfo)
            
            if isempty(exinfo(kk).RFwx)
                val(kk) = -100;
            else
                val(kk) = exinfo(kk).RFwx;
            end
        end
    case 'RF width y'
        for kk = 1:length(exinfo)
            if isempty(exinfo(kk).RFwy)
                val(kk) = -100;
            else
                val(kk) = exinfo(kk).RFwy;
            end
        end
        
    case 'RF width mean'
        
        for kk = 1:length(exinfo)
            if exinfo(kk).isc2
                val(kk) = nan;
            else
                val(kk) = exinfo(kk).RFw;
            end
        end
        
        
    case 'pref size base'
        
        for kk = 1:length(exinfo)
            if strcmp(exinfo(kk).param1, 'sz')
                val(kk) = exinfo(kk).fitparam.mu;
            else
                val(kk) = -100;
            end
        end
        
        
    case 'pref size drug'
        
        for kk = 1:length(exinfo)
            if strcmp(exinfo(kk).param1, 'sz')
                val(kk) = exinfo(kk).fitparam_drug.mu;
            else
                val(kk) = -100;
            end
        end
        
    case 'pref size diff'
        val = assignFct('pref size base' , exinfo) - ...
            assignFct('pref size drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
        
    case 'nonparam area ratio'
        val = [exinfo.nonparam_ratio];
        
        
    case 'nonparam area ratio 2nd half'
        for i = 1:length(exinfo)
            val(i) = mean(exinfo(i).ff_drug.classic_2ndhalf.spkcnt_mn) /...
                mean(exinfo(i).ff.classic_2ndhalf.spkcnt_mn);
        end
        
    case 'MI'
        val = [exinfo.MI];
        
    case 'pf stim raw base'
        for kk =1:length(exinfo)
            [~, bin_i] = max( exinfo(kk).ratemn );
            val(kk) = exinfo(kk).ratepar(bin_i);
        end
    case 'pf stim raw drug'
        for kk =1:length(exinfo)
            [~, bin_i] = max( exinfo(kk).ratemn_drug );
            val(kk) = exinfo(kk).ratepar_drug(bin_i);
        end
    case 'pf stim raw diff'
        val = assignFct('pf stim raw base' , exinfo) - ...
            assignFct('pf stim raw drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'SI base'
        for kk =1:length(exinfo)
            if isfield(exinfo(kk).fitparam, 'SI')
                val(kk) = exinfo(kk).fitparam.SI;
            else
                val(kk) = -100;
            end
        end
        
        
    case 'SI drug'
        for kk =1:length(exinfo)
            if isfield(exinfo(kk).fitparam_drug, 'SI')
                val(kk) = exinfo(kk).fitparam_drug.SI;
            else
                val(kk) = -100;
            end
        end
        
        
    case 'SI diff'
        val = assignFct('SI base' , exinfo) - ...
            assignFct('SI drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'volt'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).volt;
        end
        
    case 'gain 2D'
        val = [exinfo.gain_2D];
        
        
    case 'off 2D'
        val = [exinfo.off_2D];
        
    case 'gain co'
        val = [exinfo.gain_co];
        
    case 'off co'
        val = [exinfo.off_co];
        
    case 'dose'
        val = [exinfo.dose];
        
    case  'r2 ag'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.sub.r2_ag;
        end
        val(val<0) = 0;
        
    case 'r2 cg'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.sub.r2_cg;
        end
        val(val<0) = 0;
        
    case 'r2 rg'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.sub.r2_rg;
        end
        val(val<0) = 0;
        
        
    case 'r2 rg cg diff'
        val = assignFct('r2 rg' , exinfo) - ...
            assignFct('r2 cg' , exinfo);
        
        
    case 'r2 ag cg diff'
        val = assignFct('r2 ag' , exinfo) - ...
            assignFct('r2 cg' , exinfo);
        
    case 'a ag'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.a_ag;
        end
        
    case 'a cg'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.a_cg;
        end
        
    case 'a rg'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.a_rg;
        end
        
    case 'phase selectivity'
        for kk = 1:length(exinfo)
            if isnan(exinfo(kk).tf_f1f0)
                val(kk) = nan; %exinfo(kk).phasesel;
            else
                try
                val(kk) = exinfo(kk).tf_f1f0;
                catch
                   disp(''); 
                end
            end
%             val(kk) = exinfo(kk).phasesel;
        end
        
    case 'phase selectivity drug'
        for kk = 1:length(exinfo)
            if isnan(exinfo(kk).tf_f1f0)
                val(kk) = nan;% exinfo(kk).phasesel_drug;
            else
                val(kk) = exinfo(kk).tf_f1f0;
            end
%                 val(kk) = exinfo(kk).phasesel_drug;
        end
        
        
    case 'phase selectivity diff'
        val = assignFct('phase selectivity', exinfo) - ...
            assignFct('phase selectivity drug', exinfo);
        lab = [lab ' (base-drug)'];
        
        
    case 'rmax base'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam.rmax;
        end
        
    case 'rmax drug'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam_drug.rmax;
        end
        
    case 'rmax diff'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam.rmax - exinfo(kk).fitparam_drug.rmax;
        end
        lab = [lab ' (base-drug)'];
    
    case 'c50 base'
        val = cellfun(@(x) x.c50, {exinfo.fitparam})*100;
%         val(val>1.1) = 1.1;
    case 'c50 drug'
        val = cellfun(@(x) x.c50, {exinfo.fitparam_drug})*100;
%         val(val>1.1) = 1.1;
        
    case 'c50 diff'
        val = assignFct('c50 base', exinfo) - ...
            assignFct('c50 drug', exinfo);
        lab = [lab ' (base-drug)'];
        
        
    case 'co fit m base'
        val = cellfun(@(x) x.m, {exinfo.fitparam});
    case 'co fit m drug'
        val = cellfun(@(x) x.m, {exinfo.fitparam_drug});
        
    case 'co fit n base'
        val = cellfun(@(x) x.n, {exinfo.fitparam});
        
    case 'co fit n drug'
        val = cellfun(@(x) x.n, {exinfo.fitparam_drug});
        
    case 'co fit n diff'
        val = assignFct('co fit n base', exinfo) - ...
            assignFct('co fit n drug', exinfo);
        lab = [lab ' (base-drug)'];
        
    case 'co fit m'
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).fitparam.m;
        end
        
    case 'blank diff'
        val = assignFct('blank base', exinfo) - ...
            assignFct('blank drug', exinfo);
        lab = [lab ' (base-drug)'];
                
    case 'blank base'
        for kk = 1:length(exinfo)
            if any(exinfo(kk).ratepar > 1000)
                val(kk) = exinfo(kk).ratemn(exinfo(kk).ratepar > 1000);
            else
                val(kk) = nan; 
            end
        end
%         val(val==0) = eps;
        
    case 'blank drug'
        for kk = 1:length(exinfo)
            if any(exinfo(kk).ratepar_drug > 1000)
                val(kk) = exinfo(kk).ratemn_drug(exinfo(kk).ratepar_drug > 1000);
            else
                val(kk) = nan; 
            end            
        end
%         val(val==0) = eps;

        
    case 'BRI ACF base'
        for kk = 1:length(exinfo)
            if ~isempty(exinfo(kk).bridx) && exinfo(kk).bridx(1)<1000
                val(kk) = exinfo(kk).bridx(1);
            else
                val(kk) = -1;
            end
        end
    case  'BRI ACF drug'
        for kk = 1:length(exinfo)
            if ~isempty(exinfo(kk).bridx) && exinfo(kk).bridx(2)<1000
                val(kk) = exinfo(kk).bridx(2);
            else
                val(kk) = -1;
            end
        end
        
    case  'BRI ACF diff'
        val = assignFct('BRI ACF base', exinfo) - ...
            assignFct('BRI ACF drug', exinfo);
        
        lab = [lab ' (base-drug)'];
        
    case 'BRI ISI base'
        for kk = 1:length(exinfo)
            if ~isnan(exinfo(kk).isi_frct)
                val(kk) = exinfo(kk).isi_frct(1);
            else
                val(kk) = -1;
            end
        end
    case 'BRI ISI drug'
        for kk = 1:length(exinfo)
            if ~isnan(exinfo(kk).isi_frct)
                val(kk) = exinfo(kk).isi_frct(2);
            else
                val(kk) = -1;
            end
        end
        
    case  'BRI ISI diff'
        val = assignFct('BRI ISI base', exinfo) - ...
            assignFct('BRI ISI drug', exinfo);
        
        lab = [lab ' (base-drug)'];
        
        
    case 'SMI'
        val = (cellfun(@max, {exinfo.ratemn}) - cellfun(@max, {exinfo.ratemn_drug}))  ./ ...
            (cellfun(@max, {exinfo.ratemn}) + cellfun(@max, {exinfo.ratemn_drug}));
        
    case 'electrode depth'
        
        val = [exinfo.ed];
        
    case 'smallest response'
        
        for kk=1:length(exinfo)
            mn = min(exinfo(kk).ratemn(exinfo(kk).ratepar<1000));
            val(kk) = mn(1);
        end
                
        
    case 'smallest response drug'
        
        for kk=1:length(exinfo)
            mn = min(exinfo(kk).ratemn_drug(exinfo(kk).ratepar_drug<1000));
            val(kk) = mn(1);
        end
        
    case 'smallest response overall'
                
        val = min([assignFct('smallest response', exinfo); ...
            assignFct('smallest response drug', exinfo)]);
       
    case 'response duration'
        for kk = 1:length(exinfo)
            val(kk) =  exinfo(kk).dur;
        end
    
        
    case 'response duration drug'
        for kk = 1:length(exinfo)
            val(kk) =  exinfo(kk).dur_drug;
        end
    
        
    case 'response duration diff'
        
        val = assignFct('response duration drug', exinfo) ...
            -assignFct('response duration', exinfo);
        lab = [fctname ' (base - drug)'];
        
        
    case 'latnoise'
        for kk = 1:length(exinfo)
            if exinfo(kk).isRC
                val(kk) = exinfo(kk).noise;
            else
                val(kk) = -1;
            end
        end
        
    case 'latnoise drug'
        for kk = 1:length(exinfo)
            if exinfo(kk).isRC
                val(kk) = exinfo(kk).noise_drug;
            else
                val(kk) = -1;
            end
        end
    
    case 'latnoise diff'
        
        val = assignFct('latnoise', exinfo) - ...
            assignFct('latnoise drug', exinfo);
        
        lab = [fctname '(base-drug)'];
        
    case 'pupil size base w1'
        val = assignPPSZ(exinfo, false, 'w1');
        
    case 'pupil size base w2'
        val = assignPPSZ(exinfo, false, 'w2');
    
    case 'pupil size base w3'
        val = assignPPSZ(exinfo, false, 'w3');
   
    case 'pupil size base w4'
        val = assignPPSZ(exinfo, false, 'w4');
    
    case 'pupil size drug w1'
        val = assignPPSZ(exinfo, 1, 'w1');
    
    case 'pupil size drug w2'
        val = assignPPSZ(exinfo, 1, 'w2');
        
    case 'pupil size drug w3'
        val = assignPPSZ(exinfo, 1, 'w3');
        
    case 'pupil size drug w4'
        val = assignPPSZ(exinfo, 1, 'w4');
        
    case 'latency base'
        for kk = 1:length(exinfo)
            if size(exinfo(kk).lat,1)>1
                if isempty(exinfo(kk).pfi)
                    val(kk) = exinfo(kk).lat(2,end);
                else
                    try
                        val(kk) = exinfo(kk).lat(2,exinfo(kk).pfi);
                    catch
                        c
                    end
                end
            else
                val(kk) = exinfo(kk).lat;
            end
        end
        val(isnan(val)) = -10;
    
    case 'latency drug'
        for kk = 1:length(exinfo)
            if size(exinfo(kk).lat,1)>1
                if isempty(exinfo(kk).pfi)
                    val(kk) = exinfo(kk).lat_drug(2,end);
                else
                    val(kk) = exinfo(kk).lat_drug(2,exinfo(kk).pfi_drug);
                end
            else
                val(kk) = exinfo(kk).lat_drug;
            end
        end
        val(isnan(val)) = -10;
    
    case 'latency hmax base'
        val = [exinfo.lat2Hmax];
    
    case 'latency hmax drug'
        val = [exinfo.lat2Hmax_drug];
    
    case 'latency hmax diff'
        val =  assignFct('latency hmax drug', exinfo) - ...
            assignFct('latency hmax base', exinfo);
        lab = 'latency hmax diff (drug-base)';
        
    case 'latency diff'
        val = assignFct('latency drug', exinfo) - ...
            assignFct('latency base', exinfo);
        lab = 'latency fp diff (drug-base)';
        
    case 'latency base corrected'
        val = [exinfo.lat_c];
        
    case 'latency drug corrected'
        val = [exinfo.lat_drug_c];
    
    case 'latency diff corrected'
        val = assignFct('latency drug corrected', exinfo) - ...
            assignFct('latency base corrected', exinfo);
        lab = 'latency diff corrected (drug-base)';
        
    case 'correction factor'
        val = [exinfo.reg_slope];
    
        
    case 'predicted latency'
        val = [exinfo.reg_off] + [exinfo.reg_slope] .* assignFct('latnoise', exinfo);
    
        
    case 'predicted latency drug'
        val = [exinfo.reg_off] + [exinfo.reg_slope].*assignFct('latnoise drug', exinfo);
    
        
    case 'predicted latency diff'
        val = assignFct('predicted latency drug', exinfo) -...
            assignFct('predicted latency', exinfo);
        lab = [fctname ' (drug - base)'];
        
    case 'predicted - true latency'
        val = assignFct('predicted latency', exinfo) - [exinfo.lat2Hmax];
    
        
    case 'predicted - true latency drug'
        val = assignFct('predicted latency drug', exinfo) - [exinfo.lat2Hmax_drug];
    
        
    case 'noise correlation corrected'
        val = [exinfo.rsc]-0.0062698*[exinfo.c0geomn];
        
    case 'noise correlation drug corrected'
        val = [exinfo.rsc_drug]-0.0062698*[exinfo.c0geomn_drug];
        
    case 'noise correlation diff corrected'
        val = assignFct('noise correlation corrected', exinfo) - ...
            assignFct('noise correlation drug corrected', exinfo);
        lab = 'noise correlation corrected diff (base-drug)';
        
        
    case 'noise correlation'
        val = [exinfo.rsc];
        lab = 'noise correlation';
        
    case 'noise correlation drug'
        val = [exinfo.rsc_drug];
        lab = 'noise correlation drug';
        
    case 'noise correlation diff'
        val = assignFct('noise correlation', exinfo) - ...
            assignFct('noise correlation drug', exinfo);
        lab = 'noise correlation diff (base-drug)';
        
        
    case  'signal correlation'
        val = [exinfo.rsig];
        lab = 'signal correlation';
        
    case 'signal correlation drug'
        val = [exinfo.rsig_drug];
        lab = 'signal correlation drug';
        
    case 'signal correlation diff'
        val = [exinfo.rsig] - [exinfo.rsig_drug];
        lab = 'signal correlation diff (base-drug)';
        
    case 'noise correlation 2nd half'
        val = [exinfo.rsc_2nd];
        
    case 'noise correlation drug 2nd half'
        val = [exinfo.rsc_2nd_drug];
        
    case 'noise correlation diff 2nd half'
        val = assignFct('noise correlation 2nd half', exinfo) - ...
            assignFct('noise correlation drug 2nd half', exinfo);
        
        
    case 'noise correlation 2nd half corrected'
        val = [exinfo.rsc_2nd] - 0.0053495*[exinfo.c0geomn_2nd];
        lab = [fctname ' (estimate = 0.0053495)'];
        
    case 'noise correlation drug 2nd half corrected'
        val = [exinfo.rsc_2nd_drug] - 0.0053495*[exinfo.c0geomn_2nd_drug];
        lab = [fctname ' (estimate = 0.0053495)'];
        
    case 'noise correlation diff 2nd half corrected'
        val = assignFct('noise correlation 2nd half corrected', exinfo) - ...
            assignFct('noise correlation drug 2nd half corrected', exinfo);
        
    case 'signal correlation 2nd half'
        val = [exinfo.rsig_2nd];
        
    case 'signal correlation drug 2nd half'
        val = [exinfo.rsig_2nd_drug];
    
    case 'signal correlation diff 2nd half'
        val = [exinfo.rsig_2nd] - [exinfo.rsig_2nd_drug];
        lab = 'signal correlation 2nd half diff (base-drug)';
        
    case 'mean spike rate base'
        
        for i =1:length(exinfo)
           stim_base = exinfo(i).ratepar;
           stim_drug = exinfo(i).ratepar_drug;
           i_both = ismember(stim_base, intersect(stim_drug,stim_base));
           val(i) = mean( exinfo(i).ratemn(i_both) );
        end
        lab = 'mean spike rate base';
        
    case 'mean spike rate drug'
       for i =1:length(exinfo)
           stim_base = exinfo(i).ratepar;
           stim_drug = exinfo(i).ratepar_drug;
           i_both = ismember(stim_drug, intersect(stim_drug,stim_base));
           val(i) = mean( exinfo(i).ratemn_drug(i_both) );
        end
        lab = 'mean spike rate drug';
        
    case 'mean spike rate diff'
       
        val = assignFct('mean spike rate base', exinfo) -...
            assignFct('mean spike rate drug', exinfo);
        lab = 'mean spike rate diff (base-drug)';
        
    case 'mean spike rate variance base'
%         val = cellfun(@mean, {exinfo.ratevars},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike rate variance';
        
    case 'mean spike rate variance drug'
%         val = cellfun(@mean, {exinfo.ratevars_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike rate variance drug';
        
    case 'mean spike rate variance diff'
%         val =  cellfun(@mean, {exinfo.ratevars},...
%             'UniformOutput', 0);
%         val2 = cellfun(@mean, {exinfo.ratevars_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
%         val2 = cell2mat(val2);
%         val = val - val2;
        lab = 'mean spike rate diff';
        
    case 'mean spike count base'
%         val = cellfun(@mean, {exinfo.spkCount_mn},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike count';
        
    case 'mean spike count drug'
%         val = cellfun(@mean, {exinfo.spkCount_mn_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike count drug';
        
    case 'mean spike count diff'
%         
%         val =  cellfun(@mean, {exinfo.spkCount_mn},...
%             'UniformOutput', 0);
%         val2 = cellfun(@mean, {exinfo.spkCount_mn_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
%         val2 = cell2mat(val2);
%         val = val - val2;
        lab = 'mean spike count diff';
        
    case 'mean spike count variance base'
%         val = cellfun(@mean, {exinfo.spkCount_var},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike count variance';
        
    case 'mean spike count variance drug'
%         val = cellfun(@mean, {exinfo.spkCount_var_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
        lab = 'mean spike count variance drug';
        
    case 'mean spike count variance diff'
%         val =  cellfun(@mean, {exinfo.spkCount_var},...
%             'UniformOutput', 0);
%         val2 = cellfun(@mean, {exinfo.spkCount_var_drug},...
%             'UniformOutput', 0);
%         val = cell2mat(val);
%         val2 = cell2mat(val2);
%         val = val - val2;
        lab = 'mean spike count variance diff';
        
    case 'r2'
        val = [exinfo.r2reg];
        lab = 'regression r2';
        
    case 'gain change'

        val = [exinfo.gslope];
        lab = 'gain change';
        
    case 'additive change'
        
        val = [exinfo.yoff];
        lab = 'additive change';
        
    case 'additive change (rel)'
        val = [exinfo.yoff_rel];
        
    case 'additive change (rel & w\o zeros)'
        
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).yoff_rel_wo_null;
        end
        
    case 'gain change (rel & w\o zeros)'
        
        for kk = 1:length(exinfo)
            val(kk) = exinfo(kk).gslope_rel_wo_null;
        end
        
                
    case 'ff 2nd half blank base'
                
        for kk = 1:length(exinfo)
            blank_stim = exinfo(kk).ff.classic_2ndhalf.stim >1000;
            
            if any(blank_stim)
                val(kk) = exinfo(kk).ff.classic_2ndhalf.ff(blank_stim);
            else
                val(kk) = nan;
            end
        end
        
        
    case 'ff 2nd half blank drug'
        
        for kk = 1:length(exinfo)
            blank_stim = exinfo(kk).ff_drug.classic_2ndhalf.stim >1000;
            if any(blank_stim)
                val(kk) = exinfo(kk).ff_drug.classic_2ndhalf.ff(blank_stim);
            else
                val(kk) = nan;
            end
        end
        
        
    case 'ff 2nd half blank diff'
        
        val = assignFct('ff 2nd half blank base', exinfo) -...
            assignFct('ff 2nd half blank drug', exinfo);
        
        lab = 'ff 2nd half blank diff (base-drug) ';
        
        
    case 'ff 2nd half peak resp base'
        
        for kk = 1:length(exinfo)
            [~,peak_stim] = max(exinfo(kk).ff.classic_2ndhalf.spkcnt_mn);
            val(kk) = exinfo(kk).ff.classic_2ndhalf.ff(peak_stim);
        end
        
        
    case 'ff 2nd half peak resp drug'
        
        for kk = 1:length(exinfo)
            [~,peak_stim] = max(exinfo(kk).ff.classic_2ndhalf.spkcnt_mn);
            val(kk) = exinfo(kk).ff_drug.classic_2ndhalf.ff(peak_stim);
        end
        
    case 'ff 2nd half peak resp diff'
    
        val = assignFct('ff 2nd half peak resp base', exinfo) -...
            assignFct('ff 2nd half peak resp drug', exinfo);
        
        lab = 'ff 2nd half peak resp diff (base-drug) ';
        
    case 'ff 2nd half base'
        for kk = 1:length(exinfo)
            val(kk) = getMeanFF(exinfo(kk).ff.classic_2ndhalf,...
                ff_min_nrep, ff_min_spks);
        end
        lab = 'mean FF 2nd half base';
        
    case 'ff 2nd half drug'
        
        for kk = 1:length(exinfo)
            val(kk) = getMeanFF(exinfo(kk).ff_drug.classic_2ndhalf,...
                ff_min_nrep, ff_min_spks); 
        end
        lab = 'mean FF 2nd half drug';

        
    case 'ff 2nd half diff'
        val = assignFct('fano factor 2nd half base', exinfo) -...
            assignFct('ff 2nd half drug', exinfo);
        lab = 'ff diff (base-drug) ';
        
        
    case 'ff 2nd half peak resp diff - blank diff'
        val = assignFct('ff 2nd half peak resp diff', exinfo) -...
            assignFct('ff 2nd half blank diff', exinfo);
        
        lab = 'ff diff (peak-blank) ';
        
        
    case 'ff 20+ base'
        for kk = 1:length(exinfo)
            val(kk) = getMeanFF(exinfo(kk).ff.classic_20plus,...
                ff_min_nrep, ff_min_spks);
        end
        %         val(val==0) = nan;
        lab = 'mean FF 20+ nd half base';
        
    case 'ff 20+ drug'
        for kk = 1:length(exinfo)
            val(kk) = getMeanFF(exinfo(kk).ff_drug.classic_20plus,...
                ff_min_nrep, ff_min_spks);
        end
        lab = 'mean FF 20+ half drug';
        
        
    case 'ff 20+ diff'
        val = assignFct('ff 20+ base', exinfo) -...
            assignFct('ff 20+ drug', exinfo);
        
        lab = 'ff 20+ (base-drug) ';
        
    case 'ff 20+ blank base'
        for kk = 1:length(exinfo)
            blank_stim = exinfo(kk).ff.classic_20plus.stim >1000;
            if any(blank_stim)
                val(kk) = exinfo(kk).ff.classic_20plus.ff(blank_stim);
            else
                val(kk) = nan;
            end
        end
        
    case 'ff 20+ blank drug'
        for kk = 1:length(exinfo)
            blank_stim = exinfo(kk).ff_drug.classic_20plus.stim >1000;
            if any(blank_stim)
                val(kk) = exinfo(kk).ff_drug.classic_20plus.ff(blank_stim);
            else
                val(kk) = nan;
            end
        end
        
    case 'ff 20+ blank diff'
        val = assignFct('ff 20+ blank base', exinfo) -...
            assignFct('ff 20+ blank drug', exinfo);
        lab = 'ff 20+ blank diff (base-drug) ';
            
    case 'ff 20+ peak resp base'
        for kk = 1:length(exinfo)
            [~, maxi] = max(exinfo(kk).ff_drug.classic_20plus.spkcnt_mn);
            if any(maxi)
                val(kk) = exinfo(kk).ff.classic_20plus.ff(maxi);
            else
                val(kk) = nan;
            end
        end
        
    case 'ff 20+ peak resp drug'
        for kk = 1:length(exinfo)
            [~, maxi] = max(exinfo(kk).ff_drug.classic_20plus.spkcnt_mn);
            if any(maxi)
                val(kk) = exinfo(kk).ff_drug.classic_20plus.ff(maxi);
            else
                val(kk) = nan;
            end
        end
        
    case 'ff 20+ peak resp diff'
        val = assignFct('ff 20+ peak resp base', exinfo) -...
            assignFct('ff 20+ peak resp drug', exinfo);
        lab = 'ff 20+ peak resp diff (base-drug) ';
        
    case 'ff 20+ peak resp diff - blank diff'
        val = assignFct('ff 20+ peak resp diff', exinfo) -...
            assignFct('ff 20+ blank diff', exinfo);
        lab = 'ff 20+ diff (peak-blank) ';
       
    case 'ff base'
        for kk = 1:length(exinfo)
            exinfo(kk).id
            val(kk) = getMeanFF(exinfo(kk).ff.classic,...
                ff_min_nrep, ff_min_spks);
        end
                
    case 'ff drug'
        for kk = 1:length(exinfo)
            val(kk) = getMeanFF(exinfo(kk).ff_drug.classic,...
                ff_min_nrep, ff_min_spks);
        end
        
    case 'ff diff'
         val = assignFct('ff base', exinfo) -...
            assignFct('ff drug', exinfo);
        lab = 'ff diff (base-drug) ';
        
    case 'ff fit base'
        for kk = 1:length(exinfo)
            val(kk) = [exinfo(kk).ff.fit];
        end
        
    case 'ff fit drug'
        for kk = 1:length(exinfo)
            val(kk) = [exinfo(kk).ff_drug.classic];
        end
        
    case 'ff fit diff'
        for kk = 1:length(exinfo)
            val(kk) = [exinfo(kk).ff.fit] - [exinfo(kk).ff_drug.fit];
        end
        
    case 'wave width'
        for kk = 1:length(exinfo)
            if exinfo(kk).wdt(1)>1
                val(kk) = exinfo(kk).wdt(1);
            else
                val(kk) = nan;
            end
        end
        
    case 'pupil size'
        
    case 'spike field coherence'
        
    case 'tc ampl base'
        val = [exinfo.tcdiff];
        lab = 'tuning curve difference base';
        
    case 'tc ampl drug'
        val = [exinfo.tcdiff_drug];
        lab = 'tuning curve difference drug';
        
    case 'tc ampl diff'
        val = [exinfo.tcdiff] - [exinfo.tcdiff_drug];
        lab = 'tuning curve difference diff';
    case 'fano factor mitchel bins'
        val = [0.05 0.1 0.1778 0.31 0.5623 1 1.7783 3.10 5.6234 10 20 40];
%         val = exp( binrng(1)+mean(diff(binrng))/2 : mean(diff(binrng)) : binrng(9)+0.2+mean(diff(binrng))/2);
        lab = 'mean spike count (per 100ms)';
        
    case 'fano factor mitchel variance'
        
        funname = 'mitchel_2ndhalf';
%         funname = 'mitchel';
        
        exinfo = exinfo([exinfo.is5HT]);
        
        binrng = [0.05 0.1 0.1778 0.31 0.5623 1 1.7783 3.10 5.6234 10 20 40];
        err = zeros( length(binrng), 2 );
        val = zeros( length(binrng), 2 );
        
        %%% baseline
        val_base_var    = cellfun(@(x) x.(funname).var, {exinfo.ff}, 'UniformOutput', 0);
        val_base_var    = horzcat(val_base_var{:});
        val_base_mn     = cellfun(@(x) x.(funname).mn, {exinfo.ff}, 'UniformOutput', 0);
        val_base_mn     = horzcat(val_base_mn{:});
        
        [~, bin_i] = histc(val_base_mn, binrng); % bin according to mean values
        for kk = 1:max(bin_i) % average variance across all data with similar mean values
            % there are also units with 0 spike counts. They are ignored
            % for the fano factor analysis here, therefore the bin starts
            % with 1 instead of 0
            base{kk} = val_base_var(bin_i==kk);
            val(kk, 1) = nanmean( base{kk} );
            n_base(kk) = sum( bin_i==kk );
            err(kk, 1) = 2 * ( nanstd( base{kk} ) ./ sqrt(n_base(kk)) ); %sem
        end
        
        
        %%% drug
        val_drug_var    = cellfun(@(x) x.(funname).var, {exinfo.ff_drug}, 'UniformOutput', 0);
        val_drug_var    = horzcat(val_drug_var{:});
        val_drug_mn     = cellfun(@(x) x.(funname).mn, {exinfo.ff_drug}, 'UniformOutput', 0);
        val_drug_mn     = horzcat(val_drug_mn{:});

        
        [~, bin_i] = histc(val_drug_mn, binrng);% bin according to mean values
        for kk = 1:max(bin_i)  % average variance across all data with similar mean values
            
            drug{kk} = val_drug_var(bin_i==kk);
            val(kk, 2) = nanmean( drug{kk} );
            n_drug(kk) = sum( bin_i==kk );
            err(kk, 2) = 2 * ( nanstd( drug{kk} ) ./ sqrt(n_drug(kk)) );
        end
        
        val(isnan(val)) = 0;
        err(isnan(err)) = 0;
        
        % compare the distribution
        for kk = 1:max(bin_i)
            if all(isnan(base{kk})) || all(isnan(drug{kk}))
                p(kk) = nan;
            else
                p(kk) = ranksum(base{kk}, drug{kk});
            end
        end
        
        lab = 'spike count variance';
        addinfo = ['p-values: ' sprintf('%1.2f  \t', p) ...
            '\n n Base = ' sprintf('%1.0f   \t', n_base) ...
            '\n n 5HT = ' sprintf('%1.0f   \t', n_drug)];
        
    case 'gauss fit mu base'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.mu;
            end
        end
    
        
    case 'gauss fit mu drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.mu;
            end
        end
    
        
    case 'gauss fit mu diff'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.mu - exinfo(kk).fitparam_drug.mu;
            end
        end
        val(val>90) = 180-val(val>90);
        lab = [fctname '(base-drug)'];
        
    case 'gauss fit sig base'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.sig;
            end
        end
    
        
    case 'gauss fit sig drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.sig;
            end
        end
    
        
    case 'gauss fit sig diff'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = [exinfo(kk).fitparam.sig] - [exinfo(kk).fitparam_drug.sig];
            end
        end
        lab = [fctname '(base-drug)'];
        
        
        
         case 'gauss fit a base'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.a;
            end
        end
    
        
    case 'gauss fit a drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.a;
            end
        end
    
        
    case 'gauss fit a diff'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = [exinfo(kk).fitparam.a] - [exinfo(kk).fitparam_drug.a];
            end
        end
        lab = [fctname '(base-drug)'];
        
        
    case 'gauss fit b base'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.b;
            end
        end
    
        
    case 'gauss fit b drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.b;
            end
        end
    
        
    case 'gauss fit b diff'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = [exinfo(kk).fitparam.sig] - [exinfo(kk).fitparam_drug.b];
            end
        end
        lab = [fctname '(base-drug)'];
        
    case 'gauss ratio fit width center'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.wc;
            end
        end
    
        
    case 'gauss ratio fit width center drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.wc;
            end
        end
    
        
    case 'gauss ratio fit width surround'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.ws;
            end
        end
    
        
    case 'gauss ratio fit width surround drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.ws;
            end
        end
    
        
    case 'gauss ratio fit gain surround'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.ks;
            end
        end
    
        
    case 'gauss ratio fit gain surround drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.ks;
            end
        end
    
        
    case 'gauss ratio fit gain center'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam )
                val(kk)  = exinfo(kk).fitparam.kc;
            end
        end
    
        
    case 'gauss ratio fit gain center drug'
        for kk =1:length(exinfo)
            if ~isempty ( exinfo(kk).fitparam_drug )
                val(kk)  = exinfo(kk).fitparam_drug.kc;
            end
        end
    
        
    case 'geometric mean c1-c0' 
        for kk =1:length(exinfo)
                val(kk) = exinfo(kk).c0geomn_2nd; 
%              val(bin) = exinfo(bin).c0geomn; 
        end

        lab = [fctname ' 2nd half'];
        
end
end



function dat = assignFctXY(fcttime, fctother, datinfo)


[dat.y, dat.ylab] = assignFct(fctother, datinfo);

dat.xlab = fcttime;

switch fcttime
    
    case 'experiment time'
        
        for neuron = unique([datinfo.id])
            neuron_idx = find([datinfo.id] == neuron);
            
            [b, idx] = sort(unique([datinfo(neuron_idx).date]));
            
            for i = 1:length(b)
                dat.x(neuron_idx([datinfo(neuron_idx).date] == b(i))) = idx(i);
            end
        end
        
    case 'trial time'
        dat.x = 1:length(datinfo);
end

end



function val = assignPPSZ(exinfo, drug_flag, w_name)

if drug_flag
    fname = 'ppsz_mu_drug';
else
    fname = 'ppsz_mu';
end

for i = 1:length(exinfo)
    
    if isfield(exinfo(i).(fname), w_name)
        val(i) = exinfo(i).(fname).(w_name);
    end
    
end

end


function val = getMeanFF(ff_struct, minrep, minspk)


idx = ff_struct.stimrep >= minrep & ...
    ff_struct.spkcnt_mn > minspk ;
try
    % val = nanmean( ff_struct.stimrep(idx).*ff_struct.ff(idx) )/ ...
    %     sum(ff_struct.stimrep(idx));
    val = nanmean(ff_struct.ff(idx));
catch
   disp('') 
end
end

