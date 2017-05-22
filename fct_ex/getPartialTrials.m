function A = getPartialTrials(A)
% divides the struct A into n parts and returns the elements of the last part


n = 2;

j = ceil(length(A)/n);
A = A(j:end);

