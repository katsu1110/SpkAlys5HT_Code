function [ marker ] = markerAssignment( param1, monkey )
% assigns the marker according to the monkey
% circles indicate ma
% squares indicate ka


if strcmp(monkey, 'ma')
    marker = 'o';
    
elseif strcmp(monkey, 'ka')
    marker = 's';
    
end

