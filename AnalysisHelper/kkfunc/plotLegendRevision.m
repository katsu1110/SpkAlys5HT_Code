function h = plotLegendRevision(data, varargin)

switch nargin
    case 0
        data = 'summary';
end

switch data
    case 'nc4'
        numbers = [76,21;37,10;100,30;35,22];
    case 'nc6'
        numbers = [76,21;37,10;97,30;35,22];
    case 'nc8'
        numbers = [75,20;20,2;57,20;34,20];
    case 'ff4'
        numbers = [76, 21; 37, 10; 100, 30; 35, 22];
    case 'ff6'
        numbers = [76, 21; 36, 10; 97, 30; 35, 22];
    case 'ff8'
        numbers = [73, 20; 10, 1; 50, 10; 31, 18];
    case 'ff_peak-blank'
        numbers = [39,14; 21,5;63,17;20,15];
    case 'eye'
        numbers = [68, 20; 25, 5; 90, 25; 23, 15];
    case 'summary'
        numbers = [61,19; 1, 0; 42, 10; 11, 9;...
            63, 19; 2, 0; 45, 10; 15, 10;...
            33,13; 12,2; 54,14; 11, 9];
end

close all
sz = 20;
a = 0.5;
h = figure;
for i = 1:4
    if ~strcmp(data, 'summary')
         switch i
            case 1
                param = 'or';
                val = numbers(1,:);
            case 2
                param = 'sf';
                val = numbers(2,:);
            case 3
                param = 'co';
                val = numbers(3,:);
            case 4
                param = 'sz';
                val = numbers(4,:);
         end
    else
         switch i
            case 1
                param = 'or';
            case 2
                param = 'sf';
            case 3
                param = 'co';
            case 4
                param = 'sz';
         end
        val = numbers;
    end
    col1 = getCol4Stim(1, param);
    scatter(i-0.1, 1,  sz, 'o', 'markerfacecolor', col1, 'markeredgecolor', col1, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    col2 = getCol4Stim(0, param);
    scatter(i+0.1, 1,  sz, 'o',  'markerfacecolor', col2, 'markeredgecolor', col2, 'markerfacealpha', a, 'markeredgealpha', 1)
    hold on;
    text(i-0.1,1.4,upper(param),'FontSize',5)
    if ~strcmp(data, 'summary')
        text(i-0.18, 0.6, [num2str(val(1)) '/' num2str(val(2))],'FontSize',5)
    else
        text(i-0.18, 0.6, [num2str(val(i,1)) '/' num2str(val(i,2))],'FontSize',5)
        text(i-0.18, 0.3, [num2str(val(i+4,1)) '/' num2str(val(i+4,2))],'FontSize',5)
        text(i-0.18, 0, [num2str(val(i+8,1)) '/' num2str(val(i+8,2))],'FontSize',5)
    end
end
text(-0.25, 1, 'Symbol','FontSize',5)
% text(-0.5, 0.6, 'n','FontSize',5)
text(4.5, 0.6, 'n','FontSize',5)
text(4.5, 1, '5HT / NaCl','FontSize',5)
if strcmp(data, 'summary')
   text(-0.25, 0.6, '(Fig. A, C)', 'FontSize',5)
   text(-0.22, 0.3, '(B, D)', 'FontSize',5)
   text(-0.22, 0, '(E)', 'FontSize',5)
end
% xlim([-0.5 5])
% ylim([0.5 1.5])
whitebg('w')
% set(0,'defaultfigurecolor',[1 1 1])
set(gcf,'position', [866   817   374   161])
axis off


% switch data
%     case 'nc4'
%         numbers = [76,21;37,10;100,30;35,22];
%     case 'nc6'
%         numbers = [76,21;37,10;97,30;35,22];
%     case 'nc8'
%         numbers = [75,20;20,2;57,20;34,20];
%     case 'ff4'
%         numbers = [76, 21; 37, 10; 100, 30; 35, 22];
%     case 'ff6'
%         numbers = [76, 21; 36, 10; 97, 30; 35, 22];
%     case 'ff8'
%         numbers = [73, 20; 10, 1; 50, 10; 31, 18];
%     case 'ff_peak-blank'
%         numbers = [39,14; 21,5;63,17;20,15];
% end
% 
% close all
% sz = 20;
% a = 0.5;
% h = figure;
% for i = 1:4
%      switch i
%         case 1
%             param = 'or';
%             val = numbers(1,:);
%         case 2
%             param = 'sf';
%             val = numbers(2,:);
%         case 3
%             param = 'co';
%             val = numbers(3,:);
%         case 4
%             param = 'sz';
%             val = numbers(4,:);
%     end
%     col1 = getCol4Stim(1, param);
%     scatter(i-0.1, 1,  sz, 'o', 'markerfacecolor', col1, 'markeredgecolor', col1, 'markerfacealpha', a, 'markeredgealpha', 1)
%     hold on;
%     col2 = getCol4Stim(0, param);
%     scatter(i+0.1, 1,  sz, 'o',  'markerfacecolor', col2, 'markeredgecolor', col2, 'markerfacealpha', a, 'markeredgealpha', 1)
%     hold on;
%     text(i-0.1,1.4,upper(param),'FontSize',5)
%     text(i-0.18, 0.6, [num2str(val(1)) '/' num2str(val(2))],'FontSize',5)
% end
% text(-0.3, 1, 'Symbol','FontSize',5)
% text(0.2, 0.6, 'n','FontSize',5)
% text(4.5, 1, '5HT / NaCl','FontSize',5)
% xlim([-0.5 5])
% ylim([0.5 1.5])
% whitebg('w')
% % set(0,'defaultfigurecolor',[1 1 1])
% set(gcf,'position', [616   809   560   103])
% axis off

switch data
    case 'eye'
        savedir = 'Z:\Corinna\SerotoninPaper\Figures\Figure08_EyeData\raw_figs\legend_';
    otherwise
        savedir = 'Z:\Corinna\SerotoninPaper\Figures\Figure09_NoiseCorr_and_FF\raw_figs\legend_';
end
savefig(h, [savedir data '.fig'])
