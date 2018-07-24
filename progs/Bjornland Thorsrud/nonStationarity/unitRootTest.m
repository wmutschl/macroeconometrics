function result=unitRootTest(y,p,nlag,ex,exin)
% PIURPOSE: Test for unit roots using either Dickey Fuller (df) or
% Augmented Dickey Fuller (adf) test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% y = a time-series vector (t x 1). Should be in levels
%
% p = order of time polynomial in the null-hypothesis
%                 p = -1, no deterministic part
%                 p =  0, for constant term
%                 p =  1, for constant plus time-trend
%                 p >  1, for higher order polynomial
%
% nlags = # of lagged changes of y included     
%
% varargin = options input. String followed by argument
%   'vname','' = string. If supplied a table with test statistics will be
%   produced
%
% Output:
%
% result = structure with the following fields: 
%   .meth   = 'adf'
%   .beta   = estimate of the (autoregressive) parameters
%   .adf    = ADF t-statistic
%   .crit   = (6 x 1) vector of critical values
%           [1% 5% 10% 90% 95% 99%] quantiles    
%   .nlag   = nlag   
%   .p      = same as input
%   .nobs   = number of observations used in regression
%
% Usage: 
%
% result=unitRootTest(y,p,nlag)
%
% Notes: References: Said and Dickey (1984) 'Testing for Unit Roots in
% Autoregressive Moving Average Models of Unknown Order', 
% Biometrika, Volume 71, pp. 599-607.
%
% Written by...
% James P. LeSage, Dept of Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com
%
% The model estimated and tested is:
% Dy(t)= a + b y(t-1) + D_1(y(t-1)) + ... + D_n(y(t-n))+ e(t)
%
% The unit root test is carried out under the null hypothesis of a unit 
% root against the alternative hypothesis of no unit root. The important
% parameter is b=ø-1, where ø is the parameter in; 
% y(t)= a + ø y(t-1) + d_1(y(t-1)) + ... + d_n(y(t-n))+ e(t)
% i.e. a unit root implies b=ø-1 = 0 
%
% Once a value for the test statistic is computed it can be compared to the
% relevant critical value for the Dickey–Fuller Test. If the test statistic
% is less (this test is non symmetrical so we do not consider an absolute 
% value) than (a larger negative) the critical value, then the null
% hypothesis of a unit root is rejected and no unit root is present.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.vname=exin;
options.method='adf';
options.nvarxc=1;

% Check input
if p<-1
    error('unitRootTest:input','p less than -1 in unitRootTest');
end;
[nobs,n]=size(y);
if n>1
    error('unitRootTest:input','cannot handle a matrix -- only vectors');    
end;
if nlag<1
    error('unitRootTest:input','nlag must be larger than or equal to 1');        
end;
if (nobs - 2*nlag)+1<1
    error('unitRootTest:input','nlags too large in adf, negative dof');        
end;

func=str2func(options.method);
result=feval(func,y,p,nlag);

% Add output
result.nlag=nlag;
result.p=p;
result.(options.method)=(result.beta(1,1))/result.bstd(1,1);
result.crit=ztcrit(result.nobs,p);        
result.critc=rztcrit(result.nobs,options.nvarxc,p);        
result.meth=options.method;

if ~isempty(options.vname)
    print(result,options.vname)
    result.vname=options.vname;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result=adf(x,p,nlag)

y=tdiff(x,1);
ych=latMlag(y,nlag);
yl=latMlag(x,1);

y2=y(2+nlag:end,1);

ych=ych(2+nlag:end,:);
yl=yl(2+nlag:end,:);
x2=[yl ych];

if (p > -1)
    pt=ptrend(p,size(x2,1));
    x2 = [x2 pt];
else
    pt=[];
end;

[nobs,nvarx]=size(x2);

beta=(x2\y2)';
resid=y2-x2*beta';
so=(resid'*resid)/(nobs-nvarx);
invx=diag(inv(x2'*x2));
bstd=sqrt(so*invx);

% output
result.beta=beta;
result.bstd=bstd;
result.nobs=nobs;
result.pt=pt;
result.nvarx=nvarx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print(result,vname)

if isempty(vname)
    Vname{1}='variable1';
else
    Vname=[];
    [~,nsize]=size(vname);
    nmax = min(nsize,16); % truncate vnames to 16-characters
    Vname{1}=vname(1,1:nmax);
end; % end of nflag issue

fid=1;
fprintf(fid,'\nAugmented DF test for unit root variable: %20s \n',Vname{1});  
in.cnames = char('ADF t-statistic','# of lags','constant','time polynomial','AR(1) estimate');
in.fmt = char('16.6f','%5d','%5d','%5d','%16.6f');

if result.p<0
    c=0;
    tp=0;
else
    c=1;
    if result.p>0
        tp=result.p;
    else
        tp=0;
    end;            
end
tmp = [result.adf result.nlag c tp result.beta(1,1)+1];
mprint(tmp,in);
 
in2.cnames = char('1% Crit Value','5% Crit Value','10% Crit Value');
in2.rnames = char(' ','adf','adfc');
in2.fmt = '%16.3f';
in2.fid = fid;
mprint([result.crit(1:3)';result.critc(1:3)'],in2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
