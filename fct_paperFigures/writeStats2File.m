function stats = writeStats2File(s, fname)
% writes the staticts s attached to a plot from mugui into a text file with
% the given filename fname and reads out the main values relevant for the
% plot


% write the string to file
fileID = fopen([fname(1:end-3) 'txt'],  'wt');
fprintf(fileID, '%s\n', s);
fclose(fileID);

%%% read statistics into an array
% n 5HT
[~,idxn]  = regexp(s, 'N = '); 
i5ht = idxn(1); l = 3;
while contains( s( i5ht:i5ht+l ), ')') && l>0
    l = l-1;
end
stats.n5HT = str2double( s( i5ht:i5ht+l ) );

% n NaCl
inacl = idxn(4); l = 3;
while contains( s( inacl:inacl+l ), ')') && l>0
    l = l-1;
end
stats.nNaCl = str2double( s( inacl:inacl+l ) );



% median values
[~, idxmed]= regexp(s, '(5, 50, 95');

% median in 5HT
stats.med5HTx = getMedFromString(s, idxmed(1));
stats.med5HTy = getMedFromString(s, idxmed(2));

% median NaCl
% ... in x
stats.medNaClx = getMedFromString(s, idxmed(4));
stats.medNaCly = getMedFromString(s, idxmed(5));


% paired signrank p value
p = getSampleTest(s);
stats.p_xy5HT = p.xy5HT;
stats.p_xyNaCl = p.xyNaCl;
stats.p_x5HTvsxNaCl = p.x5HTvsNaCl;
stats.p_y5HTvsyNaCl = p.y5HTvsNaCl;
stats.p_xy5HTvsxyNaCl = p.xy5HTvsxyNaCl;


% correlation rho and p value
[stats.spearmanRho5HT, stats.spearmanP5HT, stats.spearmanRhoNaCl, stats.spearmanPNaCl] = ...
    getCorrStats(s, 'Spearman: rho=');
[stats.pearsonRho5HT, stats.pearsonP5HT, stats.pearsonRhoNaCl, stats.pearsonPNaCl] = ...
    getCorrStats(s, 'Pearson: rho=');
end




function [rho5HT, p5HT, rhoNaCl, pNaCl] = getCorrStats(s, fname)
% extract correlation rho and p value
l = 5;
l2 = 8;

[~,idxp]  = regexp(s, fname);  idxp = idxp+1;
rho5HT = str2double( s( idxp(1):idxp(1)+l ) );
rhoNaCl = str2double( s( idxp(2):idxp(2)+l ) );

if rho5HT>=0
    p5HT = str2double( s( idxp(1)+l+3:idxp(1)+l+l2+4 ) );
else
    p5HT = str2double( s( idxp(1)+l+4:idxp(1)+l+l2+4 ) );
end

if rhoNaCl>=0
    pNaCl = str2double( s( idxp(2)+l+3:idxp(2)+l2+l+4 ) );
else
    pNaCl = str2double( s( idxp(2)+l+4:idxp(2)+l2+l+4 ) );
end

end



function p = getSampleTest(s)


alpha = 0.05; % alpha level
l = 8; 

% test indeces
%%% indeces for paired tests
[~,idxPairedSR]  = regexp(s, 'paired signrank p=');  idxPairedSR = idxPairedSR+1; % non-parametric
[~,idxPairedTT]  = regexp(s, 'paired t-test p=');  idxPairedTT = idxPairedTT+1; % parametric 

p.xy5HT = str2double( s( idxPairedSR(1):idxPairedSR(1)+l ) );
p.xyNaCl = str2double( s( idxPairedSR(2):idxPairedSR(2)+ l) );


%%% indeces for two sample tests
[~,idxp2sTT]  = regexp(s, '2-sample ttest p='); idxp2sTT = idxp2sTT+1; % parametric
[~,idxp2sW]  = regexp(s, 'wilcoxon p='); idxp2sW = idxp2sW+1;  % non-parametric

p.x5HTvsNaCl = str2double( s( idxp2sW(1):idxp2sW(1)+l ) );
p.y5HTvsNaCl = str2double( s( idxp2sW(2):idxp2sW(2)+l ) );
p.xy5HTvsxyNaCl = str2double( s( idxp2sW(3):idxp2sW(3)+l ) );


% entries containing the test p-Value for normal distribution
% if p>alpha, than the sample is assumed normally distributed
% [~,idxp]  = regexp(s, 'test for normal dist p=');  idxp = idxp+1;
% kstest_p_5HTx = str2double( s( idxp(1):idxp(1)+l ) );
% kstest_p_5HTy = str2double( s( idxp(2):idxp(2)+l ) );
% kstest_p_5HTyx = str2double( s( idxp(3):idxp(3)+l ) );
% kstest_p_NaClx = str2double( s( idxp(4):idxp(3)+l ) );
% kstest_p_NaCly = str2double( s( idxp(5):idxp(4)+l ) );
% kstest_p_NaClxy = str2double( s( idxp(6):idxp(4)+l ) );

%%% test for normal distribution concerning both samples (5HT and NaCl)
% if kstest_p_5HTx>alpha && kstest_p_NaClx > alpha
%     p.x5HTvsNaCl = str2double( s( idxp2sTT(1):idxp2sTT(1)+ l) ); % assuming normal distribution
% else
%     p.x5HTvsNaCl = -str2double( s( idxp2sW(1):idxp2sW(1)+l ) );
% end

% if kstest_p_5HTy>alpha && kstest_p_NaCly > alpha
%     p.y5HTvsNaCl = str2double( s( idxp2sTT(2):idxp2sTT(2)+l ) );% assuming normal distribution
% else
%     p.y5HTvsNaCl = -str2double( s( idxp2sW(2):idxp2sW(2)+l ) );
% end 

% if kstest_p_5HTyx>alpha && kstest_p_NaClxy
%     p.xy5HTvsNaCl = str2double( s( idxp2sTT(3):idxp2sTT(3)+l ) );% assuming normal distribution
% else
%     p.xy5HTvsNaCl = -str2double( s( idxp2sW(3):idxp2sW(3)+l ) );
% end 




%%% test for normal distribution concerning paired samples (5HT or NaCl)
% [~,idxp]  = regexp(s, 'normal dist X-Y h=');  idxp = idxp+1;
% kstest_h_5HTxy = str2double( s( idxp(1):idxp(1)+1 ) );
% kstest_h_NaClxy = str2double( s( idxp(2):idxp(2)+1 ) );


% if h is 1, the hypotheses that claims normal distribution is rejected
% i.e., if h is 0, the normal distribution is assumed
% if ~kstest_h_5HTxy
%     p.xy5HT = str2double( s( idxPairedTT(1):idxPairedTT(1)+l) ); % assuming normal distribution
% else
%     p.xy5HT = -str2double( s( idxPairedSR(1):idxPairedSR(1)+l ) );
% end


% if ~kstest_h_NaClxy
%     p.xyNaCl = str2double( s( idxPairedTT(2):idxPairedTT(2)+l ) );% assuming normal distribution
% else
%     p.xyNaCl = -str2double( s( idxPairedSR(2):idxPairedSR(2)+ l) );
% end



end






function med = getMedFromString(s, idxmed)

s2 = s(idxmed:idxmed+100);
idx = strfind(s2, '/');
med = str2double( s2( idx(1)+1: idx(2)-1 ) );

end



