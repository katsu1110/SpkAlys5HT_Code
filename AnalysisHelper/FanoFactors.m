function [ ff, ff_fit, mitchel, church ] = FanoFactors( ex, stim, Mn, Vars, nrep, param1 )
% all kind of fano factor calculations are performed
% - casual mn/vars
% - power law fit
% - mitchell et al 
% - churchland et al 


%------------------------------------------------------ Fano Factor
ff.ff = ( Vars ./ Mn );
ff.stim = stim; %stimulus values
ff.spkcnt_mn = Mn;
ff.spkcnt_var  = Vars;
ff.stimrep = nrep;

[FitPa, ~,~,~,~] = FitPower( Mn , Vars); %HNs function
ff_fit = FitPa.exponent;

%------------------------------------------------------- Mitchel FF
mitchel = [];
% [mitchel.mn, mitchel.var] = Mitcheletal(ex.Trials, param1); 

%------------------------------------------------------- Churchl FF
church = [];
% church = Churchlandetal( ex.Trials );
    
end

