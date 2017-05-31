function s = getStringID(id, monkey)
% returns the corresponding string prefix to identify matching entries in
% the sorting table

id = floor(id);
if nargin == 1 || id < 215 || (id>224 && id<226) || (id>233 && id~=248)|| strcmp(monkey, 'ka')
    
    if id <10
        s = '000';
    elseif id <100
        s = '00';
    elseif id < 1000
        s = '0';
    else
        s = '';
    end
    
    s = [s num2str(id)];
    
else
    
    s = num2str(id);
    
end
end