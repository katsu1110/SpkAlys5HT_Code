function A = getPartialTrials(A, n)
% divides the struct A into n parts and returns the elements of the last part
% n is an optinal argument. the default value is 2.
%
% @CL 


if nargin==1
    n = 2;
end

j = ceil(length(A)/n);
A = A(j:end);

