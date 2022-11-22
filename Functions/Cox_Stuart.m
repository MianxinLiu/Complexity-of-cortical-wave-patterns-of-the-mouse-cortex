%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COX-STUART TEST   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[H,p_value]=Cox_Stuart(V,alpha)
%%%%%%%%%%%%%%%%%
%%% Performs  Cox-Stuart test of the null hypothesis of trend
%%% absence in the vector V,  against the alternative of trend. 
%%% The result of the test is returned in H = 1 indicates
%%% a rejection of the null hypothesis at the alpha significance level. H = 0 indicates
%%% a failure to reject the null hypothesis at the alpha significance level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS
%V = time series [vector]
%alpha =  significance level of the test [scalar]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The significance level of a test is a threshold of probability a agreed
%to before the test is conducted. A typical value of alpha is 0.05. If the p-value of a test is less than alpha,
%the test rejects the null hypothesis. If the p-value is greater than alpha, there is insufficient evidence 
%to reject the null hypothesis. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS
%H = test result [1] Reject of Null Hypthesis [0] Insufficient evidence to reject the null hypothesis
%p_value = p-value of the test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% From Matlab Help %%%%%%%%%%%%%%%
%The p-value of a test is the probability, under the null hypothesis, of obtaining a value
%of the test statistic as extreme or more extreme than the value computed from
%the sample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% References 
%Cox, D. R., and A. Stuart (1955), Some quick tests for trend in location and
%dispersion, Biometrika, 42, 80–95.
%Conover, W. J. (1999), Practical Nonparametric Statistics, 3rd ed., John
%Wiley, New York.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simone Fatichi -- simonef@dicea.unifi.it
%   Copyright 2009
%   $Date: 2009/10/03 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=reshape(V,length(V),1); 
alpha = alpha/2; %
n=length(V);
if mod(n,2) == 0
    V1=V(1:n/2);
    V2=V((n/2+1):end);
else
    V1=V(1:(n-1)/2);
    V2=V(((n+1)/2+1):end);
end
S=sign(V1-V2);
xp=length(find(S==1));
xn=length(find(S==-1));
xz=length(find(S==0));
X=min(xp,xn);
N=xp+xn-xz;
if N<=0
    Z=0;
else
    p=0.5; q= 0.5;
    Z=(abs(X- N*p)-0.5)/sqrt(N*p*q);
end
p_value=2*(1-normcdf(abs(Z),0,1));
pz=norminv(1-alpha,0,1);
H=abs(Z)>pz; 
return