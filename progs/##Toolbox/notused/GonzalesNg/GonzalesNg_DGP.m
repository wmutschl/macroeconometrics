clear all;
bigt = 200; init = 0;
capt=init*bigt;
nobs=bigt+capt;
nvar=3;
sigma=eye(2);
   dx =   randn(nobs,nvar-1);
   e  =   randn(nobs,1) ;              %/* Generate random walks           */
   u  = randn(nobs,1);
   u =  recserar(u,0,1);
   x  = recserar(dx,zeros(1,nvar-1),ones(1,nvar-1)) ;
      %/* Accumulate the innovations      */
     %/* x is a non-cointegrated system  */
    %/* of 2 I(1) variables             */
   %/* cointegrating vector: (1,-1,-1) */
   y  =  x(:,1) + 2*x(:,2) + e ;      %/* dgp */
   data=[y x];
   ddata=diff(data,1);

data=trimr(data,capt,0);
ddata=trimr(ddata,capt,0);

[h,pValue,stat,cValue,mles] = jcitest(data,'model','H2','lags',2,'test','trace','display','full')
%let alpha[3,1]=-1,1,1.;
%gama=zeros(3,1);
%gama[1]=1.0;
%nr1=1;
%retp(data,ddata);