function [c,ceq] = mycon(X, M)
c = zeros(length(M)-1,1);
for i=2:length(M);
c(i-1) = M{i}(X);
end
ceq = [];