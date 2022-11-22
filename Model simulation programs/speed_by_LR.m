function [TheoryWaveVol,TheoryWaveVol_largeWave]=speed_by_LR(k,alpha,f)
% compute theoretical wave speed by linear regression

% 1. use the whole range
q=k/(2*pi);
Y=f';
A=[ones(length(q),1),q'];
sig=max(abs(alpha)); %
w=exp(-alpha.^2/(2*sig^2))';
A=w.*A;
Y=w.*Y;
b=(A'*A)\(A'*Y);
TheoryWaveVol=-b(2);

% 2. only use the range with f>0
id_positive=find(f>0);
q=q(id_positive);alpha=alpha(id_positive);f=f(id_positive);
Y=f';
A=[ones(length(q),1),q'];
sig=max(abs(alpha)); %
w=exp(-alpha.^2/(2*sig^2))';
A=w.*A;
Y=w.*Y;
b=(A'*A)\(A'*Y);
TheoryWaveVol_largeWave=-b(2);

