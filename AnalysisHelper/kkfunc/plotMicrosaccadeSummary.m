function [para] = plotMicrosaccadeSummary(dat)
%% summarize microsaccade data given by Corinna's exinfo analysis
% INPUT: dat = get(gcf, 'UserData') after plotting a scatter in 'mugui'
%
%
% written by Katsuhisa (12.07.17)
% +++++++++++++++++++++++++++++++++++++++++++++++++++

% initialization
for m = 1:2
    for i = 1:4
        switch i
            case 1
                a = 1;
                b = 1;
            case 2
                a = 1;
                b = 2;
            case 3
                a = 2;
                b = 1;
            case 4
                a = 2;
                b = 2;
        end
        para.animal(m).ser(a).exp(b).amp = [];
        para.animal(m).ser(a).exp(b).peakv = [];
        para.animal(m).ser(a).exp(b).duration = [];
        para.animal(m).ser(a).exp(b).angle = [];
    end
end

% store data
for i = 1:length(dat.expInfo)
    % is the drug 5HT?
    if dat.is5HT(i)==1
        d = 1;
    else
        d = 2;
    end
    
    % is the animal mango?
    if dat.expInfo(i).ismango==1
        m = 1;
    else
        m = 2;
    end
    
    % baseline
    para.animal(m).ser(d).exp(1).amp = [para.animal(m).ser(d).exp(1).amp, dat.expInfo(i).microsac_amplitude];
    para.animal(m).ser(d).exp(1).peakv = [para.animal(m).ser(d).exp(1).peakv, dat.expInfo(i).microsac_peakv];
    para.animal(m).ser(d).exp(1).duration = [para.animal(m).ser(d).exp(1).duration, dat.expInfo(i).microsac_duration];
    para.animal(m).ser(d).exp(1).angle = [para.animal(m).ser(d).exp(1).angle, dat.expInfo(i).microsac_angle];
    
    % drug
    para.animal(m).ser(d).exp(2).amp = [para.animal(m).ser(d).exp(2).amp, dat.expInfo(i).microsac_amplitude_drug];
    para.animal(m).ser(d).exp(2).peakv = [para.animal(m).ser(d).exp(2).peakv, dat.expInfo(i).microsac_peakv_drug];
    para.animal(m).ser(d).exp(2).duration = [para.animal(m).ser(d).exp(2).duration, dat.expInfo(i).microsac_duration_drug];
    para.animal(m).ser(d).exp(2).angle = [para.animal(m).ser(d).exp(2).angle, dat.expInfo(i).microsac_angle_drug];
end

% visualize
figure;
for i = 1:4
    switch i
        case 1
            name = 'amp';
            xlab = 'log amplitude';
        case 2
            name = 'peakv';
            xlab = 'log peak v';
        case 3
            name = 'duration';
            xlab = 'log duration';
        case 4
            name = 'angle';
            xlab = 'angle';
    end
    for k = 1:4
        switch k
            case 1
                ylab = 'Mango (vs 5HT)';
                a = 1;
                b = 1;
            case 2                
                ylab = 'Mango (vs NaCl)';
                a = 1;
                b = 2;
            case 3                
                ylab = 'Kaki (vs 5HT)';
                a = 2;
                b = 1;
            case 4                
                ylab = 'Kaki (vs NaCl)';
                a = 2;
                b = 2;
        end
        subplot(4,4,i+4*(k-1))
        if i <= 3
            x1 = log(para.animal(a).ser(b).exp(1).(name));
            x2 = log(para.animal(a).ser(b).exp(2).(name));
%             x1 = sqrt(para.animal(a).ser(b).exp(1).(name));
%             x2 = sqrt(para.animal(a).ser(b).exp(2).(name));
        else
            x1 = para.animal(a).ser(b).exp(1).(name);
            x2 = para.animal(a).ser(b).exp(2).(name);
        end
        histogram(x1)
        hold on;
        histogram(x2)
        med1 = median(x1);
        med2 = median(x2);
        hold on;
        yy = get(gca, 'YLim');
        plot([med1 med1], yy, '-b')
        hold on;
        plot([med2 med2], yy, '-r')
        set(gca, 'box', 'off')
        set(gca, 'TickDir', 'out')       
        
        % stats
        [~,p1] = ttest2(x1,x2);
        p2 = ranksum(x1,x2);
        title(['ppar = ' num2str(p1) ', pnon = ' num2str(p2)])
                      
        if i == 1
            ylabel(ylab)  
        end
        if k == 4
            xlabel(xlab)
        end       
    end    
    
end