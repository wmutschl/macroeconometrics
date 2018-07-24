function f=f_SMOD(B0inv,SIGMAUHAT,LRMat,J,Acomp,Rshort,Rshorth,Rlong,selSR,selSRh,selLR,idxSRh)
% =======================================================================
% Evaluates the system of nonlinear equations vech(SIGMAUHAT) = vech(B0inv*B0inv')
% subject to the short-run or long-run restrictions.
% =======================================================================
% f=f_SMOD(B0inv,SIGMAUHAT,LRMat,Rshort,Rlong,selSR,selLR)
% -----------------------------------------------------------------------
% INPUTS
%   - B0inv     : candidate for short-run impact matrix. [nvars x nvars]
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
%	- LRMat     : total long-run impact matrix. [nvars x nvars] 
%                 - for VAR model A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)
%                 - for VECM model Upsylon = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))
%   - Rshort    : Matrix of short-run restrictions on impact matrix B_0^{-1}. [nvars x nvars]
%   - Rlong     : Matrix of long-run restrictions. [nvars x nvars]
%   - selSR     : Index for short-run restrictions
%   - selLR     : Index for long-run restrictions
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value, see below
% -----------------------------------------------------------------------
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu
tmp = transpose(SIGMAUHAT-B0inv*B0inv');
vechBBt_SIGU = tmp(triu(true(size(tmp))));

JAhJtB0inv = [];
for h = idxSRh
    rshorth = Rshorth{h};
    tmp = (J*Acomp^h*J')*B0inv;
    JAhJtB0inv = [JAhJtB0inv; tmp(selSRh{h}) - rshorth(selSRh{h})];
end

LRMAT = LRMat*B0inv;
f=[vechBBt_SIGU;
   B0inv(selSR) - Rshort(selSR);
   JAhJtB0inv;
   LRMAT(selLR) - Rlong(selLR);
   ];
    