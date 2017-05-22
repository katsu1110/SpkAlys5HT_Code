function stats = getFixWindStats( exinfo )
% get min max and median window


ID = unique([exinfo.id]);


for kk = 1:length(ID)
    
    temp = exinfo(find([exinfo.id]==ID(kk), 1, 'first'));
    load(temp.fname);
   [ fwh(kk), fww(kk) ] = getFixWind( ex );

end

% compute the min, median and max of each 
stats.H = [min(fwh), median(fwh), max(fwh)];
stats.W = [min(fww), median(fww), max(fww)];
stats.H_raw = fwh;
stats.W_raw = fww;

end