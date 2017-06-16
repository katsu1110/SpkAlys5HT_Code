function exinfo = addLM2Struct( exinfo )


load('linreg.mat');

lat_std_t = 200:400;
k = 1;

for i = 1:length( exinfo )
    
    
    exinfo(i).reg_slope = 0;
    exinfo(i).reg_off  = 0;
    exinfo(i).noise = 0;
    exinfo(i).noise_drug = 0;
    
    exinfo(i).lat_c = 0;
    exinfo(i).lat_drug_c = 0;
    
    exinfo(i).predicted_change = 0;
    exinfo(i).fig_latjackknife = '';
        
    
    if exinfo(i).isRC
        
        
        if exinfo(i).ismango
            idx = exinfo(i).id == [simlm.id] & strcmp({simlm.monkey}, 'ma');
        else
            idx = exinfo(i).id == [simlm.id]+0.5 & strcmp({simlm.monkey}, 'ka');
        end
        
        if sum(idx) == 1
            exinfo(i).reg_slope = simlm(idx).std_slope;
            exinfo(i).reg_off  = simlm(idx).std_off;
            exinfo(i).noise = mean(sqrt(exinfo(i).resvars(lat_std_t)));
            exinfo(i).noise_drug = mean(sqrt(exinfo(i).resvars_drug(lat_std_t)));
            
            exinfo(i).lat_c = exinfo(i).lat2Hmax - exinfo(i).noise*exinfo(i).reg_slope;
            exinfo(i).lat_drug_c = exinfo(i).lat2Hmax_drug - exinfo(i).noise_drug*exinfo(i).reg_slope;
            
            exinfo(i).predicted_lat_change = (exinfo(i).reg_off + exinfo(i).noise*exinfo(i).reg_slope) - ...
                (exinfo(i).reg_off + exinfo(i).noise_drug*exinfo(i).reg_slope);
            
            exinfo(i).fig_latjackknife = simlm(idx).figname;
        else
            
            disp(exinfo(i).figname);
            
        end
    end
    
end

