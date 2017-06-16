function [ idx ] = cmpFiles( expInfo, varargin )
% loads a list with file names and compares whhich are represented in the
% expInfo data struct


j = 1;
while j<=length(varargin)
    switch varargin{j}
        case 'sfn_all'
            exnamefile = textscan(fopen('Z:\Corinna\filenames\allSFN_poster.txt', 'r'), '%s');
            exnamefile = exnamefile{1}(1:2:end);
            hn_files = cellfun(@(x, i0, i1) x(i0:i1), exnamefile, strfind(exnamefile, 'mango'),...
                strfind(exnamefile, '.mat'), 'UniformOutput', false);

        case 'sfn_or'
            exnamefile = textscan(fopen('Z:\Corinna\filenames\ORISFN_poster.txt', 'r'), '%s');
            exnamefile = exnamefile{1}(1:2:end);
            hn_files = cellfun(@(x, i0, i1) x(i0:i1), exnamefile, strfind(exnamefile, 'mango'),...
                strfind(exnamefile, '.mat'), 'UniformOutput', false);

        case 'rc_all'
            load('hn_rc.mat', 'hn');
            hn_files = strrep(hn, 'J:\data\', '');
    end
    j = j +1;
end




cl_files = cellfun(@(x, i0, i1) x(i0:i1-1), {expInfo.fname}, strfind({expInfo.fname}, 'mango'),...
    strfind({expInfo.fname}, '.mat'), 'UniformOutput', false)';


idx = ismember(cl_files, hn_files); 

end

