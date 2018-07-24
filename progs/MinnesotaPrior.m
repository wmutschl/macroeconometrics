function SigmaPi = MinnesotaPrior(OverallTightness,CrossEqTightness,HarmonicLagDecay,OmegaVec,k)
% MinnesotaPrior: Constructs the covariance matrix for the prior when the Minnesota prior
%                 is used.
%
% USAGE:
%
%       SigmaPi = MinnesotaPrior(OverallTightness,CrossEqTightness,HarmonicLagDecay,OmegaVec,p,k)
%
% REQUIRED INPUT:  OverallTightness (scalar), with the overall tightness hyperparameter.
%
%                  CrossEqTightness (scslar), with the cross-equation tightness hyperpartameter.
%
%                  HarmonicLagDecay (scalar), with the harmonic lag decay hyperparameter.
%
%                  OmegaVec (vector) with the variances of the residuals.
%
%                  k (integer), the lag order of the VAR
%
% REQUIRED OUTPUT: SigmaPi (matrix), with the resulting prior variances.
%
% NOTE: The parameters on lags are assumed to be organized as follows:
%
%               PiMatrix = [Pi(1) ... Pi(k)]
%
% The output is then the prior covariance matrix of vec(PiMatrix), where vec is the column stacking
% operator. Let Pi(i,j;k) ber element (i,j) in Pi(k). The prior variance of this element is now
%
%          var(Pi(i,j;l)) = OverallTightness/(l^HarmonicLagDecay),           if i=j
%
%                           OverallTightness*CrossEqTightness*OmegaVec(i)/OmegaVec(j)*(l^HarmonicLagDecay)
%                                                                            otherwise.
%
%
%                       Written by: Anders Warne
%                                   New Area Wide Model Project
%                                   DG-R/MPR
%                                   European Central Bank (ECB)
%                                   Email: anders.warne@ecb.europa.eu
%                                   Copyright © 2006-2015 European Central Bank.
%
%                       First version: January 12, 2007.
%                        This version: January 12, 2015.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LICENSE INFORMATION:
%
%      YADA is free software: you can redistribute it and/or modify
%      it under the terms of the GNU General Public License as published by
%      the Free Software Foundation, either version 3 of the License, or
%      (at your option) any later version.
%
%      This program is distributed in the hope that it will be useful,
%      but WITHOUT ANY WARRANTY; without even the implied warranty of
%      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%      GNU General Public License for more details.
%
%      You should have received a copy of the GNU General Public License
%      along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%      YADA is released under the GNU General Public License, Version 3,
%      29 June 2007 <http://www.gnu.org/licenses/>. The current release of
%      the program was last modified by the ECB on the "This version" date
%      above.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGELOG:
%
% * 16-01-2007: Fixed a bug when dealing with the scale factor OmegaVec.
%
% * 08-02-2007: Updated the documentation.
%
% * 13-11-2007: Updated the documentation.
%
% * 23-05-2008: Updated the documentation.
%
% * 19-12-2008: Updated the documentation.
%
% * 02-05-2009: Updated the documentation.
%
% * 05-01-2010: Updated the documentation.
%
% * 05-01-2011: Updated the documentation.
%
% * 02-01-2012: Updated the documentation.
%
% * 02-01-2013: Updated the documentation.
%
% * 24-01-2014: Updated the documentation.
%
% * 12-01-2015: Updated the documentation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

p = length(OmegaVec);
if size(OmegaVec,1)<p;
   OmegaVec = OmegaVec';
end;
%
% add the CrossEqTightness parameter to 
%
SigmaPii = ones(p,p)*CrossEqTightness;
for i=1:p;
   SigmaPii(i,i) = 1;
end;
%
% add the OverallTightness and the scale factor
%
SigmaPii = OverallTightness*SigmaPii.*(OmegaVec*(ones(p,1)./OmegaVec)');
%
% now it remains to deal with harmonic lag decay
%
SigmaPi = zeros(p,p*k);
SigmaPi(:,1:p) = SigmaPii;
for j=2:k;
   %
   % add
   %
   SigmaPi(:,1+(p*(j-1)):p*j) = (1/(j^HarmonicLagDecay))*SigmaPii;
end;
%
% now we take the vectorization into account
%
SigmaPi = diag(vec(SigmaPi));

%
% end of MinnesotaPrior.m
%
