function [m_se_mA0,m_se_mPhiOne,m_tv_mA0,m_tv_mPhiOne] = GetBootstrapStdError_svec_var(mmA0,mmPhiOne,mA0,mPhiOne,K,mB_Res,mC1_Res)
% /*
% **	Computes bootstrap standard errors and t-values for
% **	contemp. and long run impact matrix
% **	
% **	Input: 	mmA0:  		K^2 x M, each column holds vectorized bootstrapped contemp. impact matrix 
% **						(M: # of bootstrap replications)			
% **
% **			mmPhiOne:  	K^2 x M, each column holds vectorized bootstrapped long run impact matrix 
% **			
% **			mA0:   		K x K, point estimate of contemp. impact matrix	
% **
% **			mPhiOne:	K x K, of point estimate of long run impact matrix
% **
% **			mB_Res:		K x K matrix, zeros indicate zero restrictions in contemp. impact matrix
% **
% **			mC1_Res:	K x K matrix, zeros indicate zero restrictions in long run impact matrix
% **
% **
% **	Output: m_se_mA0:		K x K, bootstrap std errors of contemp. impact matrix	
% **
% **			m_se_mPhiOne 	K x K, bootstrap std errors of long run impact matrix	
% **
% **			m_tv_mA0:		K x K, bootstap t-values of contemp. impact matrix
% **			
% **		 m_tv_mPhiOne:		K x K, bootstrap t-values of long run impact matrix	
% **			
% **	
% **
% */

 
m_se_mA0     = reshape(sqrt(mean(((mmA0-vec(mA0))^2)')),K,K)';			%/* compute bootstrap standard errors for A0 	*/
m_se_mPhiOne = reshape(sqrt(mean(((mmPhiOne-vec(mPhiOne))^2)')),K,K)';	%/* compute bootstrap standard errors for C1 	*/

m_tv_mA0 	   = mA0./m_se_mA0 ; 	  	 
m_tv_mPhiOne  = mPhiOne./m_se_mPhiOne;
 
vBsel =  find(vec(mB_res)==0);
vC1sel = find(vec(mC1_res)==0);

v_se_mPhiOne = vec(m_se_mPhiOne);
v_tv_mPhiOne = vec(m_tv_mPhiOne);

v_se_mA0 = vec(m_se_mA0);
v_tv_mA0 = vec(m_tv_mA0);
     
if isempty(vC1sel);
    m_se_mPhiOne = reshape(v_se_mPhiOne,K, K)';
    m_tv_mPhiOne = reshape(v_tv_mPhiOne,K, K)';   
else  
    v_se_mPhiOne(vC1sel) = zeros(size(vC1sel,1),1);
    v_tv_mPhiOne(vC1sel) = zeros(size(vC1sel,1),1);
    m_se_mPhiOne = reshape(v_se_mPhiOne,K, K)';
    m_tv_mPhiOne = reshape(v_tv_mPhiOne,K, K)';   
end
	 

if ismiss(vBsel);        
    m_se_mA0 = reshape(v_se_mA0,K, K)';
    m_tv_mA0 = reshape(v_tv_mA0,K, K)';
else	 
    v_se_mA0(vBsel) = zeros(size(vBsel,1),1);
    v_tv_mA0(vBsel) = zeros(size(vBsel,1),1);
    m_se_mA0 = reshape(v_se_mA0,K, K)';
    m_tv_mA0 = reshape(v_tv_mA0,K, K)';
end