function plotGLMser(exinfo, param1, varargin)

switch nargin
        case 1
                param1 = 0;
end

% select stmulus parameter
inc = [];
if ~isequal(param1, 0)
        for i = 1:length(exinfo)
                if isequal(exinfo(i).param1, param1)
                        inc = [inc i];
                end
        end
end

% vectorize
b_pf_drug = [exinfo(inc).b_pf_drug];
b_pf_pupil = [exinfo(inc).b_pf_pupil];
b_pf_intr = [exinfo(inc).b_pf_intr];
b_upf_drug = [exinfo(inc).b_upf_drug];
b_upf_pupil = [exinfo(inc).b_upf_pupil];
b_upf_intr = [exinfo(inc).b_upf_intr];




close all
figure;

subplot(2,1,1)
errorbar([1 2 3], [mean(b_pf_drug) mean(b_pf_pupil) mean(b_pf_intr)], ...
        [std(b_pf_drug) std(b_pf_pupil) std(b_pf_intr)],'o')
hold on;

subplot(2,1,2)
errorbar([1 2 3], [mean(b_upf_drug) mean(b_upf_pupil) mean(b_upf_intr)], ...
        [std(b_upf_drug) std(b_upf_pupil) std(b_upf_intr)],'o')

% for i = 1:length(b_pf_drug)
%         subplot(2,1,1)
%         plot([1 2 3], [b_pf_drug(i) b_pf_pupil(i) b_pf_intr(i)],'-o','linewidth',0.5) 
%         hold on;
%         
%         subplot(2,2,2)
%         plot([1 2 3], [b_upf_drug(i) b_upf_pupil(i) b_upf_intr(i)],'o','linewidth',0.5) 
%         hold on;
% end