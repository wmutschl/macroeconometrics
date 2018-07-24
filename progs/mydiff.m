function y=mydiff(data,m)

%input: data --> vector or matrix
%       m --> lag at which to difference
%output: y --> data in first differences

y=data(m+1:end,:)-data(1:end-m,:);