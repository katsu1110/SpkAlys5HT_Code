function h = setPaperPlotsProps(varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



figwidth = 17.4;                     %# figure width

j = 1;
while j<=length(varargin)
    switch varargin{j}
        case 'wdt'
            figwidth = varargin{j+1};
            j = j+1;
    end
    j = j+1;
end



global uni_k uni_r xarg_sc xarg_reg xarg_tc xarg_hist fsz_add;

% figure settings
inchf = 2.54;              
figwidth_in = figwidth/inchf;        %# left/right margins from page borders
h = figure('Units', 'inches', 'Position', [5 3.5 figwidth_in figwidth_in]);


% uni_k = [50 65 75]/255; % anthrazit
% uni_r = [165 30 55]/255; % uni red
% uni_k = [65 90 140]/255; % dark blue
% uni_r = [228,26,28]/255; % red
uni_r = 'r';
uni_k = 'k';


alpha = 0.7;
fsz = 8;
sz = 10;

xarg_sc = {'fsz', fsz, 'alpha', alpha, 'sz', sz, 'unity'};
xarg_reg = {'fsz', fsz, 'alpha', alpha, 'sz', sz, 'cross'};
xarg_tc = {'fsz', fsz, 'sz', 15, 'square', 'alpha', alpha};
xarg_hist = {'fsz', fsz, 'sz', 2, 'txtsz', fsz-2, 'alpha', alpha+0.1};
fsz_add = 0; % increase the font size with this value to emphasize labels and titles



end

