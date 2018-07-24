%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of Random Numbers with Laplace distribution %  
%             with MATLAB Implementation                 %
%                                                        %
% Author: M.Sc. Eng. Hristo Zhivomirov          05/01/15 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = randl(m, n)

% function: x  = randl(m, n)
% m - number of matrix rows
% n - number of matrix columns
% x - matrix with Laplacian distributed numbers with mu = 0 and sigma = 1

% generation of a numbers with Uniform distribution
u1 = rand(m, n);
u2 = rand(m, n);

% generation of a numbers with Laplace distribution
x = log(u1./u2);

end
