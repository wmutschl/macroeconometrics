function [Bh,e,xtx,xty] = ols(y,phi)
% ols: estimate a system of equations: [Bh,e,xtx,xty] = ols(y,phi)
%               Y(T*nvar) = XB + u, X: T*k, B: k*nvar.
%    where Bh: the estimated B; column: nvar; row: number of r.h.s. variables.
%          e:  estimated residual e = y -xBh,  T*nvar
%          xtx:  X'X: k-by-k
%          xty:  X'Y: k-by-nvar
%          phi:  X; T-by-k; column: number of r.h.s. variables (including
%                                                          deterministic terms)
%          y:    Y: T-by-nvar
%
%          See also sye.m and syed.m
% Copyright (C) 1997-2012 Tao Zha
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% ** setup of orders and lengths **
[u d v]=svd(phi,0); %trial
%xtx = phi'*phi;      % X'X, k*k (ncoe*ncoe)
vd=v.*(ones(size(v,2),1)*diag(d)'); %trial
dinv = 1./diag(d);    % inv(diag(d))
vdinv=v.*(ones(size(v,2),1)*dinv'); %trial
xtx=vd*vd';
xtxinv = vdinv*vdinv';
%xty = phi'*y;        % X'Y
uy = u'*y; %trial
xty = vd*uy; %trial
%Bh = xtx\xty;        %inv(X'X)*(X'Y), k*m (ncoe*nvar).
Bh = xtxinv*xty;
%e = y - phi*Bh;      % from Y = XB + U, e: (T-lags)*nvar
e = y - u*uy;