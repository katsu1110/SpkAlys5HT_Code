function A = readSpikingTable( fname )
% reads the spike sorting csv file named fname and returns a struct
% containing the partial mokey file name and the isolation quality
%
%
% @CL 22.11.2016

if ~isempty(strfind(fname, 'Kaki'))
    A = importfile(fname);
else
    A = importfile2(fname);
end
    
for kk = 1:length(A)
    
    q= A(kk).isolationquality;
    if isempty(strfind(q,'c2')) && isempty(strfind(q,'c1'))
        A(kk).isoQc1 = str2double(q);
        A(kk).isoQc2 = 6;

            
    elseif isempty(strfind(q,'c2')) && ~isempty(strfind(q,'c1'))
        
        [~,si] = regexp(q,'=', 'match');
        
        A(kk).isoQc1 = str2double(q(si(1)+1:end));
        A(kk).isoQc2 = 6;

    else
        
        [~,si] = regexp(q,'=', 'match');
        [~,ei] = regexp(q,';', 'match');
        
        A(kk).isoQc1 = str2double(q(si(1)+1:ei-1));
        
        if isempty(strfind(q,'c2'))
            A(kk).isoQc2 = 6;
        else
            A(kk).isoQc2 = str2double(q(si(2)+1:end));
           
        end
    end

    if isnan( A(kk).isoQc2)
        A(kk).isoQc2 = 6;
    end      
    if isnan( A(kk).isoQc1) 
        A(kk).isoQc1 = 6;
    end      
    
end







