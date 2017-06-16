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



% swith between potential function names and call the retrieve the
% information from the linked field in exinfo
switch fctname
        
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
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).spkqual_base;
        end
        
    case 'isolation quality drug'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).spkqual_drug;
        end
        
    case 'isolation quality max'
        for bin = 1:length(exinfo)
            val(bin) = max([exinfo(bin).spkqual_base, exinfo(bin).spkqual_drug]);
        end
        
        
    case 'blank/or ratio base'
        for bin = 1:length(exinfo)
            
            if isempty(exinfo(bin).respratio)
                val(bin) = -1;
            else
                val(bin) = exinfo(bin).respratio;
            end
        end
        
    case 'blank/or ratio drug'
        for bin = 1:length(exinfo)
            
            if isempty(exinfo(bin).respratio_drug)
                val(bin) = -1;
            else
                val(bin) = exinfo(bin).respratio_drug;
            end
        end
        
    case 'blank/or ratio diff'
        
        val = assignFct('blank/or ratio base' , exinfo) - ...
            assignFct('blank/or ratio drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'eccentricity'
        val = [exinfo.ecc];
        
    case 'RF ecc corr'
        for bin = 1:length(exinfo)
            
            if isempty(exinfo(bin).RFw_corr)
                val(bin) = -100;
            else
                val(bin) = exinfo(bin).RFw_corr;
            end
        end
        
    case 'RF width x'
        for bin = 1:length(exinfo)
            
            if isempty(exinfo(bin).RFwx)
                val(bin) = -100;
            else
                val(bin) = exinfo(bin).RFwx;
            end
        end
    case 'RF width y'
        for bin = 1:length(exinfo)
            if isempty(exinfo(bin).RFwy)
                val(bin) = -100;
            else
                val(bin) = exinfo(bin).RFwy;
            end
        end
        
    case 'RF width mean'
        
        for bin = 1:length(exinfo)
            if isempty(exinfo(bin).RFw)
                val(bin) = -100;
            else
                val(bin) = exinfo(bin).RFw;
            end
        end
        
        
    case 'pref size base'
        
        for bin = 1:length(exinfo)
            if strcmp(exinfo(bin).param1, 'sz')
                val(bin) = exinfo(bin).fitparam.mu;
            else
                val(bin) = -100;
            end
        end
        
        
    case 'pref size drug'
        
        for bin = 1:length(exinfo)
            if strcmp(exinfo(bin).param1, 'sz')
                val(bin) = exinfo(bin).fitparam_drug.mu;
            else
                val(bin) = -100;
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
        for bin =1:length(exinfo)
            [~, bin_i] = max( exinfo(bin).ratemn );
            val(bin) = exinfo(bin).ratepar(bin_i);
        end
    case 'pf stim raw drug'
        for bin =1:length(exinfo)
            [~, bin_i] = max( exinfo(bin).ratemn_drug );
            val(bin) = exinfo(bin).ratepar_drug(bin_i);
        end
    case 'pf stim raw diff'
        val = assignFct('pf stim raw base' , exinfo) - ...
            assignFct('pf stim raw drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'SI base'
        for bin =1:length(exinfo)
            if isfield(exinfo(bin).fitparam, 'SI')
                val(bin) = exinfo(bin).fitparam.SI;
            else
                val(bin) = -100;
            end
        end
        
        
    case 'SI drug'
        for bin =1:length(exinfo)
            if isfield(exinfo(bin).fitparam_drug, 'SI')
                val(bin) = exinfo(bin).fitparam_drug.SI;
            else
                val(bin) = -100;
            end
        end
        
        
    case 'SI diff'
        val = assignFct('SI base' , exinfo) - ...
            assignFct('SI drug' , exinfo);
        lab = [fctname ' (base-drug)'];
        
    case 'volt'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).volt;
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
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.sub.r2_ag;
        end
        val(val<0) = 0;
        
    case 'r2 cg'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.sub.r2_cg;
        end
        val(val<0) = 0;
        
    case 'r2 rg'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.sub.r2_rg;
        end
        val(val<0) = 0;
        
        
    case 'r2 rg cg diff'
        val = assignFct('r2 rg' , exinfo) - ...
            assignFct('r2 cg' , exinfo);
        
        
    case 'r2 ag cg diff'
        val = assignFct('r2 ag' , exinfo) - ...
            assignFct('r2 cg' , exinfo);
        
    case 'a ag'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.a_ag;
        end
        
    case 'a cg'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.a_cg;
        end
        
    case 'a rg'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.a_rg;
        end
        
    case 'phase selectivity'
        for bin = 1:length(exinfo)
%             if isnan(exinfo(bin).tf_f1f0)
%                 val(bin) = exinfo(bin).phasesel;
%             else
%                 val(bin) = exinfo(bin).tf_f1f0;
%             end
            val(bin) = exinfo(bin).phasesel;
        end
        
    case 'phase selectivity drug'
        for bin = 1:length(exinfo)
%             if isnan(exinfo(bin).tf_f1f0)
%                 val(bin) = exinfo(bin).phasesel_drug;
%             else
%                 val(bin) = exinfo(bin).tf_f1f0;
%             end
                val(bin) = exinfo(bin).phasesel_drug;
        end
        
        
    case 'phase selectivity diff'
        val = assignFct('phase selectivity', exinfo) - ...
            assignFct('phase selectivity drug', exinfo);
        lab = [lab ' (base-drug)'];
        
        
    case 'rmax base'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam.rmax;
        end
        
    case 'rmax drug'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam_drug.rmax;
        end
        
    case 'rmax diff'
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam.rmax - exinfo(bin).fitparam_drug.rmax;
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
        for bin = 1:length(exinfo)
            val(bin) = exinfo(bin).fitparam.m;
        end
    case 'blank diff'
        val = assignFct('blank base', exinfo) - ...
            assignFct('blank drug', exinfo);
        lab = [lab ' (base-drug)'];
        
        
        
        for kk = 1:length(exinfo)
            
            xb = exinfo(kk).ratemn(exinfo(kk).ratepar > 1000);
            xd = exinfo(kk).ratemn_drug(exinfo(kk).ratepar_drug > 1000);
            
            if ~isempty(xb) && ~isempty(xd) 
                if xb>= .5 && xd >0
                    val(kk) = log(xd/xb);
                else
                    val(kk) = nan; 
                end
            else
                val(kk) = nan;
            end
            
        end
        
        
    case 'blank base'
        for kk = 1:length(exinfo)
            if any(exinfo(kk).ratepar > 1000)
                val(kk) = exinfo(kk).ratemn(exinfo(kk).ratepar > 1000);
            else
                val(kk) = nan; 
            end
        end
        
    case 'blank drug'
        for kk = 1:length(exinfo)
            if any(exinfo(kk).ratepar_drug > 1000)
                val(kk) = exinfo(kk).ratemn_drug(exinfo(kk).ratepar_drug > 1000);
            else
                val(kk) = nan; 
            end            
        end
        
    case 'BRI ACF base'
        for bin = 1:length(exinfo)
            if ~isempty(exinfo(bin).bridx) && exinfo(bin).bridx(1)<1000
                val(bin) = exinfo(bin).bridx(1);
            else
                val(bin) = -1;
            end
        end
    case  'BRI ACF drug'
        for bin = 1:length(exinfo)
            if ~isempty(exinfo(bin).bridx) && exinfo(bin).bridx(2)<1000
                val(bin) = exinfo(bin).bridx(2);
            else
                val(bin) = -1;
            end
        end
        
    case  'BRI ACF diff'
        val = assignFct('BRI ACF base', exinfo) - ...
            assignFct('BRI ACF drug', exinfo);
        
        lab = [lab ' (base-drug)'];
        
    case 'BRI ISI base'
        for bin = 1:length(exinfo)
            if ~isnan(exinfo(bin).isi_frct)
                val(bin) = exinfo(bin).isi_frct(1);
            else
                val(bin) = -1;
            end
        end
    case 'BRI ISI drug'
        for bin = 1:length(exinfo)
            if ~isnan(exinfo(bin).isi_frct)
                val(bin) = exinfo(bin).isi_frct(2);
            else
                val(bin) = -1;
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
        
        for bin=1:length(exinfo)
            
            mn = min(exinfo(bin).ratemn);
            val(bin) = mean(mn(1));
            
        end
        
        lab = fctname;
        
        
    case 'smallest response drug'
        
        for bin=1:length(exinfo)
            mn = min(exinfo(bin).ratemn_drug);
            val(bin) = mean(mn(1));
        end
        
        lab = fctname;
        
    case 'smallest response overall'
        
        
        val = min([assignFct('smallest response', exinfo); ...
            assignFct('smallest response drug', exinfo)]);
        
        lab = fctname;
        
    case 'response duration'
        for bin = 1:length(exinfo)
            val(bin) =  exinfo(bin).dur;
        end
        lab = fctname;
        
    case 'response duration drug'
        for bin = 1:length(exinfo)
            val(bin) =  exinfo(bin).dur_drug;
        end
        lab = fctname;
        
    case 'response duration diff'
        
        val = assignFct('response duration drug', exinfo) ...
            -assignFct('response duration', exinfo);
        lab = [fctname ' (base - drug)'];
        
        
    case 'latnoise'
        for bin = 1:length(exinfo)
            if exinfo(bin).isRC
                val(bin) = exinfo(bin).noise;
            else
                val(bin) = -1;
            end
        end
        lab = fctname;
        
    case 'latnoise drug'
        for bin = 1:length(exinfo)
            if exinfo(bin).isRC
                val(bin) = exinfo(bin).noise_drug;
            else
                val(bin) = -1;
            end
        end
        lab = fctname;
        
    case 'latnoise diff'
        
        val = assignFct('latnoise', exinfo) - ...
            assignFct('latnoise drug', exinfo);
        
        lab = [fctname '(base-drug)'];
        
    case 'pupil size base w1'
        val = assignPPSZ(exinfo, false, 'w1');
        lab = fctname;
        
    case 'pupil size base w2'
        val = assignPPSZ(exinfo, false, 'w2');
        lab = fctname;
        
    case 'pupil size base w3'
        val = assignPPSZ(exinfo, false, 'w3');
        lab = fctname;
        
    case 'pupil size base w4'
        val = assignPPSZ(exinfo, false, 'w4');
        lab = fctname;
        
        
    case 'pupil size drug w1'
        val = assignPPSZ(exinfo, 1, 'w1');
        lab = fctname;
        
    case 'pupil size drug w2'
        val = assignPPSZ(exinfo, 1, 'w2');
        lab = fctname;
        
    case 'pupil size drug w3'
        val = assignPPSZ(exinfo, 1, 'w3');
        lab = fctname;
        
    case 'pupil size drug w4'
        val = assignPPSZ(exinfo, 1, 'w4');
        lab = fctname;
        
    case 'latency base'
        for bin = 1:length(exinfo)
            if size(exinfo(bin).lat,1)>1
                if isempty(exinfo(bin).pfi)
                    val(bin) = exinfo(bin).lat(2,end);
                else
                    try
                        val(bin) = exinfo(bin).lat(2,exinfo(bin).pfi);
                    catch
                        c
                    end
                end
            else
                val(bin) = exinfo(bin).lat;
            end
        end
        val(isnan(val)) = -10;
        lab = fctname;
        
    case 'latency drug'
        for bin = 1:length(exinfo)
            if size(exinfo(bin).lat,1)>1
                if isempty(exinfo(bin).pfi)
                    val(bin) = exinfo(bin).lat_drug(2,end);
                else
                    val(bin) = exinfo(bin).lat_drug(2,exinfo(bin).pfi_drug);
                end
            else
                val(bin) = exinfo(bin).lat_drug;
            end
        end
        val(isnan(val)) = -10;
        lab = fctname;
        
    case 'latency hmax base'
        val = [exinfo.lat2Hmax];
        lab = fctname;
        
    case 'latency hmax drug'
        val = [exinfo.lat2Hmax_drug];
        lab = fctname;
        
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
        lab = fctname;
        
    case 'latency drug corrected'
        val = [exinfo.lat_drug_c];
        lab = fctname;
        
    case 'latency diff corrected'
        val = assignFct('latency drug corrected', exinfo) - ...
            assignFct('latency base corrected', exinfo);
        lab = 'latency diff corrected (drug-base)';
        
    case 'correction factor'
        val = [exinfo.reg_slope];
        lab = fctname;
        
    case 'predicted latency'
        val = [exinfo.reg_off] + [exinfo.reg_slope] .* assignFct('latnoise', exinfo);
        lab = fctname;
        
    case 'predicted latency drug'
        val = [exinfo.reg_off] + [exinfo.reg_slope].*assignFct('latnoise drug', exinfo);
        lab = fctname;
        
    case 'predicted latency diff'
        val = assignFct('predicted latency drug', exinfo) -...
            assignFct('predicted latency', exinfo);
        lab = [fctname ' (drug - base)'];
        
    case 'predicted - true latency'
        val = assignFct('predicted latency', exinfo) - [exinfo.lat2Hmax];
        lab = fctname;
        
    case 'predicted - true latency drug'
        val = assignFct('predicted latency drug', exinfo) - [exinfo.lat2Hmax_drug];
        lab = fctname;
        
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
            w = exinfo(i).nrep;
            val(i) = sum( exinfo(i).ratemn .* w  )/sum(w);
        end
        lab = 'weighted mean spike rate base';
        
    case 'mean spike rate drug'
       for i =1:length(exinfo)
            w = exinfo(i).nrep_drug;
            val(i) = sum( exinfo(i).ratemn_drug .* w  ) / sum(w);
        end
        lab = 'weighted mean spike rate drug';
        
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
%         for i = 1:length(exinfo)
%            if min(exinfo(i).ratemn)>0
%                 val(i) = exinfo(i).gslope;
%            else
%                val(i) = nan;
%            end
%         end
        val = [exinfo.gslope];
        lab = 'gain change';
        
    case 'additive change'
        val = [exinfo.yoff];
        lab = 'additive change';
        
    case 'additive change (rel)'
        val = [exinfo.yoff_rel];
        
    case 'gain change (rel)'
        val = [exinfo.gslope_rel];
        
    case 'fano factor 2nd half base'
        for bin = 1:length(exinfo)
            base = exinfo(bin).ff.classic_2ndhalf;
            drug = exinfo(bin).ff_drug.classic_2ndhalf;
            
            trials = getPartialTrials(exinfo(bin).trials_c1);
            trials_drug = getPartialTrials(exinfo(bin).trials_c1_drug);
            par1 = unique(trials.param);
            par2 = unique(trials_drug.param);
            
            s2 = ismember(par1, par2);
            s1 = ismember(par2, par1);
            idx = drug.spkcnt_mn(s1) < 2.5 | base.spkcnt_mn(s2) < 2.5 | ...
                drug.stimrep(s1) < 4 | base.stimrep(s2) < 4 ;
            
            
            base.ff(idx) = nan;
            
            idx2 = isnan(base.ff);
            base.stimrep(idx2) = nan;
            
            val(bin) = nansum(base.ff .*base.stimrep)...
                / nansum(base.stimrep);
        end
        val(val==0) = nan;
        lab = 'weighted fano factor 2nd half base';
        
    case 'fano factor 2nd half drug'
        
        for bin = 1:length(exinfo)
            
            base = exinfo(bin).ff.classic_2ndhalf;
            drug = exinfo(bin).ff_drug.classic_2ndhalf;
            
            trials = getPartialTrials(exinfo(bin).trials_c1);
            trials_drug = getPartialTrials(exinfo(bin).trials_c1_drug);
            par1 = unique(trials.param);
            par2 = unique(trials_drug.param);
            
            s2 = ismember(par1, par2);
            s1 = ismember(par2, par1);
            idx = drug.spkcnt_mn(s1) < 2.5 | base.spkcnt_mn(s2) < 2.5 | ...
                drug.stimrep(s1) <4 | base.stimrep(s2) <4 ;
            
            drug.ff(idx) = nan;
            
            idx2 = isnan(drug.ff);
            drug.stimrep(idx2) = nan;
            
            val(bin) = nansum(drug.ff .* drug.stimrep )...
                / nansum(drug.stimrep);
            
        end
        val(val==0) = nan;
        lab = 'weighted fano factor 2nd half drug';

        
    case 'fano factor 2nd half diff'
        val = assignFct('fano factor 2nd half base', exinfo) -...
            assignFct('fano factor 2nd half drug', exinfo);
        
        lab = 'weighted fano factor diff (base-drug) ';
        
    case 'fano factor base'
        for bin = 1:length(exinfo)
            exinfo(bin).id
            val(bin) = nansum(exinfo(bin).ff.classic.ff .*...
                exinfo(bin).ff.classic.stimrep)...
                / nansum(exinfo(bin).ff.classic.stimrep);
        end
        val(val==0) = nan;

        lab = 'weighted fano factor base';
                
    case 'fano factor drug'
        for bin = 1:length(exinfo)
            val(bin) = nansum(exinfo(bin).ff_drug.classic.ff .*...
                exinfo(bin).ff_drug.classic.stimrep )...
                / nansum(exinfo(bin).ff_drug.classic.stimrep);
        end
        val(val==0) = nan;

        lab = 'weighted fano factor drug';
        
    case 'fano factor diff'
         val = assignFct('fano factor base', exinfo) -...
            assignFct('fano factor drug', exinfo);
        lab = 'weighted fano factor diff (base-drug) ';
        
    case 'fano factor fit base'
        for bin = 1:length(exinfo)
            val(bin) = [exinfo(bin).ff.fit];
        end;
        
    case 'fano factor fit drug'
        for bin = 1:length(exinfo)
            val(bin) = [exinfo(bin).ff_drug.classic];
        end;
        
    case 'fano factor fit diff'
        for bin = 1:length(exinfo)
            val(bin) = [exinfo(bin).ff.fit] - [exinfo(bin).ff_drug.fit];
        end
        
    case 'wave width'
        for bin = 1:length(exinfo)
            if exinfo(bin).wdt(1)>1
                val(bin) = exinfo(bin).wdt(1);
            else
                val(bin) = nan;
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
        for bin = 1:max(bin_i) % average variance across all data with similar mean values
            % there are also units with 0 spike counts. They are ignored
            % for the fano factor analysis here, therefore the bin starts
            % with 1 instead of 0
            base{bin} = val_base_var(bin_i==bin);
            val(bin, 1) = nanmean( base{bin} );
            n_base(bin) = sum( bin_i==bin );
            err(bin, 1) = 2 * ( nanstd( base{bin} ) ./ sqrt(n_base(bin)) ); %sem
        end
        
        
        %%% drug
        val_drug_var    = cellfun(@(x) x.(funname).var, {exinfo.ff_drug}, 'UniformOutput', 0);
        val_drug_var    = horzcat(val_drug_var{:});
        val_drug_mn     = cellfun(@(x) x.(funname).mn, {exinfo.ff_drug}, 'UniformOutput', 0);
        val_drug_mn     = horzcat(val_drug_mn{:});

        
        [~, bin_i] = histc(val_drug_mn, binrng);% bin according to mean values
        for bin = 1:max(bin_i)  % average variance across all data with similar mean values
            
            drug{bin} = val_drug_var(bin_i==bin);
            val(bin, 2) = nanmean( drug{bin} );
            n_drug(bin) = sum( bin_i==bin );
            err(bin, 2) = 2 * ( nanstd( drug{bin} ) ./ sqrt(n_drug(bin)) );
        end
        
        val(isnan(val)) = 0;
        err(isnan(err)) = 0;
        
        % compare the distribution
        for bin = 1:max(bin_i)
            if all(isnan(base{bin})) || all(isnan(drug{bin}))
                p(bin) = nan;
            else
                p(bin) = ranksum(base{bin}, drug{bin});
            end
        end
        
        lab = 'spike count variance';
        addinfo = ['p-values: ' sprintf('%1.2f  \t', p) sprintf('\n')...
            'n Base = ' sprintf('%1.0f   \t', n_base) sprintf('\n')...
            '  n 5HT = ' sprintf('%1.0f   \t', n_drug)];
        
    case 'gauss fit mu base'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.mu;
            end
        end
        lab = fctname;
        
    case 'gauss fit mu drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.mu;
            end
        end
        lab = fctname;
        
    case 'gauss fit mu diff'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.mu - exinfo(bin).fitparam_drug.mu;
            end
        end
        val(val>90) = 180-val(val>90);
        lab = [fctname '(base-drug)'];
        
    case 'gauss fit sig base'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.sig;
            end
        end
        lab = fctname;
        
    case 'gauss fit sig drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.sig;
            end
        end
        lab = fctname;
        
    case 'gauss fit sig diff'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = [exinfo(bin).fitparam.sig] - [exinfo(bin).fitparam_drug.sig];
            end
        end
        lab = [fctname '(base-drug)'];
        
        
        
         case 'gauss fit a base'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.a;
            end
        end
        lab = fctname;
        
    case 'gauss fit a drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.a;
            end
        end
        lab = fctname;
        
    case 'gauss fit a diff'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = [exinfo(bin).fitparam.a] - [exinfo(bin).fitparam_drug.a];
            end
        end
        lab = [fctname '(base-drug)'];
        
        
    case 'gauss fit b base'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.b;
            end
        end
        lab = fctname;
        
    case 'gauss fit b drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.b;
            end
        end
        lab = fctname;
        
    case 'gauss fit b diff'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = [exinfo(bin).fitparam.sig] - [exinfo(bin).fitparam_drug.b];
            end
        end
        lab = [fctname '(base-drug)'];
        
    case 'gauss ratio fit width center'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.wc;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit width center drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.wc;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit width surround'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.ws;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit width surround drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.ws;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit gain surround'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.ks;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit gain surround drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.ks;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit gain center'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam )
                val(bin)  = exinfo(bin).fitparam.kc;
            end
        end
        lab = fctname;
        
    case 'gauss ratio fit gain center drug'
        for bin =1:length(exinfo)
            if ~isempty ( exinfo(bin).fitparam_drug )
                val(bin)  = exinfo(bin).fitparam_drug.kc;
            end
        end
        lab = fctname;
        
    case 'geometric mean c1-c0' 
        for bin =1:length(exinfo)
                val(bin) = exinfo(bin).c0geomn_2nd; 
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


