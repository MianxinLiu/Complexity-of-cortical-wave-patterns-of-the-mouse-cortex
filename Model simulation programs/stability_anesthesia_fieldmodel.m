function [k,alpha,f]=stability_anesthesia_fieldmodel(p,branch)
% get the dispersion relation of the steady states
% p is the anesthetic parameter
% branch = 'top' or 'bottom' for seeing the steady states in upper or lower
% branch in the bifurcation graph


% parameter setting

I_sc=0.018; % subcortical firing strength (/ms)
[Tau_e,Tau_i] = deal(40);			% E,I neuron time constant (ms)
[ge,gi]=deal(0.156,0.875);	% E,I synaptic gain strengths , ACh,GABA effect (mV*ms)
tau_di = 20; % IPSP decay time  (ms)
tau_de=5; % EPSP decay time (ms)
[Ve_rev, Vi_rev,Ve_rest, Vi_rest]=deal(0,-70,-62.5,-64);% reversal,resting potential for E,I  (mV)
[Ne_cc,Ne_local,Ni_local]=deal(200,85,120);%column-column E;intracolumn E,I;
[Qe_max,Qi_max,theta_e,theta_i,sigma_e,sigma_i]=deal(0.03,0.06,-58.5,-58.5,3,5);% sigmoid maximum firing rates (/ms),threshold (mV), 'width' (mV)

delta_lambda_i=0.08; % contribution factor
delta_lambda_e=0.01; % contribution factor
lambda_i=1+delta_lambda_i*p; % scale factor
lambda_e=1-delta_lambda_e*p; % scale factor
delta_D=3/7; % contribution factor
lambda_D=1-delta_D*p; % scale factor
Di=0.07*lambda_D;
De = Di*0.01;		% excitatory diffusion set at 1% of inhibitory value

v = 0.056;				% axonal conduction velocity (mm/ms)
r = 0.1;			% characteristic length scale for axons connectivity (mm)
gamma_e_cc=v/r;  % cortical-cortical EPSP decay rate (/ms)

% auxiliary parameters
Cse=pi/(sqrt(3) * sigma_e);Csi=pi/(sqrt(3) * sigma_i);
if strcmp(branch,'top')
    branch=3; % branch == 3 ==> top branch;
elseif strcmp(branch,'bottom')
    branch=1; % branch == 1 ==> bottom branch
end
[Phi_e_roots,Phi_i_roots,Ve_roots,Vi_roots,EIR] = findroots();
Phi_e=Phi_e_roots(branch);
Phi_i=Phi_i_roots(branch);
Ve=Ve_roots(branch);
Vi=Vi_roots(branch);
EIR=EIR(branch);

Cse=pi/(sqrt(3) * sigma_e);Csi=pi/(sqrt(3) * sigma_i);
Ee=exp(Cse*(theta_e-Ve));Ei=exp(Csi*(theta_i-Vi));
J=zeros(8);
J(1,1)=-(1+lambda_e*ge*Phi_e+lambda_i*gi*Phi_i)/Tau_e;
J(1,3)=(lambda_e*ge/Tau_e)*(Ve_rev-Ve);
J(1,5)=(lambda_i*gi/Tau_e)*(Vi_rev-Ve);
J(2,2)=-(1+lambda_e*ge*Phi_e+lambda_i*gi*Phi_i)/Tau_i;
J(2,3)=(lambda_e*ge/Tau_i)*(Ve_rev-Vi);
J(2,5)=(lambda_i*gi/Tau_i)*(Vi_rev-Vi);
J(3,4)=1;
J(4,1)=(1/tau_de)^2*Ne_local*Qe_max*Cse*Ee/(1+Ee)^2;
J(4,3)=-(1/tau_de)^2;J(4,4)=-2/tau_de;J(4,7)=(1/tau_de)^2*Ne_cc;
J(5,6)=1;
J(6,2)=(1/(tau_di*lambda_i)^2)*Ni_local*Qi_max*Csi*Ei/(1+Ei)^2;
J(6,5)=-1/(tau_di*lambda_i)^2;J(6,6)=-2/(tau_di*lambda_i);
J(7,8)=1;
J(8,1)=gamma_e_cc^2*Qe_max*Cse*Ee/(1+Ee)^2;
J(8,7)=-gamma_e_cc^2;
J(8,8)=-2*gamma_e_cc;

k=0:0.01:20;alpha=k;f=k;
for i=1:length(k)
    J_temp=J;
    J_temp(1,1)=J_temp(1,1)-De*k(i)^2/Tau_e;
    J_temp(2,2)=J_temp(2,2)-Di*k(i)^2/Tau_i;
    J_temp(8,7)=J_temp(8,7)-(v*k(i))^2;
    EigenValues=eig(J_temp);
    [~,main]=max(real(EigenValues));
    alpha(i)=real(EigenValues(main));f(i)=imag(EigenValues(main))/(2*pi);
end
alpha=alpha*1000;f=f*1000;

fig=figure();
set(gcf,'Position',[100,50,400,400]);
left_color=[0,0,1];right_color =[1,0,0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left;
plot(k,alpha,'-','LineWidth',2,'Color',[0,0,1]);hold on;
ylim([-20,3]);
xlabel('wave number k (/ 10 mm)'); % since the whole region length is Lx=10 mm
ylabel('Eig (/s)');
yyaxis right;
plot(k,f,'--','LineWidth',2,'Color',[1,0,0]);hold on;
ylim([0,4]);
ylabel('frequency (Hz)');
title(sprintf('p=%.1f',p));


    function [Phi_e,Phi_i,Ve,Vi,EIratio] = findroots()
        % Find steady states
        n_roots = 0;
        Nsearch = 8000; Nsearch_max = 8*128000; % if no roots found, keep increasing search density until Nsearch_max reached
        while ((n_roots == 0 || n_roots == 2) && Nsearch <= Nsearch_max) % NB: We treat zero roots or 2 roots as error condition (expect 1 or 3 roots)
            Qe_1 = linspace(0, Qe_max, Nsearch)';
            err_Qe = Qe_resid(Qe_1);
            chs_Qe = err_Qe(1:end-1) .* err_Qe(2:end); % Count the no. of Qe roots (look for sign-change in residuals)
            chs_Qe_index = find(chs_Qe < 0); % the left one's index of sign-changing
            n_roots = length(chs_Qe_index);
            if (n_roots == 0 || n_roots == 2)
                Nsearch = 2*Nsearch;
            end
        end
        if n_roots >= 1 % Form a bracketting interval for each root, each row is an interval
            brack_Qe = [Qe_1(chs_Qe_index)  Qe_1(chs_Qe_index + 1)]; % Qe_1 is column
        else
            error('failed to find a root');
        end
        [Qe,Qi,Ve,Vi,Phi_e,Phi_i,EIratio] = deal(ones(n_roots, 1)); % find each roots
        for j = 1: n_roots
            fun = @ (Qe_1) Qe_resid(Qe_1);
            Qe(j) = fzero(fun, brack_Qe(j,:));
            Ve(j) = invQsige(Qe(j));
            Qi(j) = EqA(Qe(j), Ve(j));
            Vi(j) = invQsigi(Qi(j));
            Phi_e(j)=(Ne_cc+Ne_local)*Qe(j)+I_sc;
            Phi_i(j)= Ni_local*Qi(j);
            EIratio(j)=abs(lambda_e*ge*Phi_e(j)*(Ve_rev-Ve(j)+Ve_rev-Vi(j))/...
                (lambda_i*gi*Phi_i(j)*(Vi_rev-Ve(j)+Vi_rev-Vi(j))));
        end
        
        
        %------------------------------------------------------------------------
        function err_Qe = Qe_resid(Qe_1)     % error residual in Qe
            Ve_1 = invQsige(Qe_1);
            Qi_2 = EqA(Qe_1, Ve_1);
            Vi_2 = invQsigi(Qi_2);
            Qe_3 = EqB(Qi_2, Vi_2);
            err_Qe = Qe_3 - Qe_1;
        end
        %------------------------------------------------------------------------
        function Qi_2 = EqA(Qe_1, Ve_1) % Given Qe_1 and Ve_1, use Eq.A to compute Qi_2.
            Qi_2 = (Ve_1 - Ve_rest  -(Ve_rev - Ve_1)*lambda_e*ge...
                .* ((Ne_cc + Ne_local) * Qe_1 + I_sc))./ ...
                ((Vi_rev - Ve_1)*lambda_i*gi* Ni_local);
            Qi_2(Qi_2 < 0 | Qi_2 > Qi_max) = NaN;
        end
        %------------------------------------------------------------------------
        function Qe_3 = EqB(Qi_2, Vi_2)   % Given Qi_2 and Vi_2, use Eq.B to compute Qe_3.
            Qe_3 = ((Vi_2 - Vi_rest - ((Vi_rev - Vi_2)*lambda_i*gi...
                * Ni_local .* Qi_2)) ./ ((Ve_rev - Vi_2)*lambda_e*ge) - ...
                I_sc) / (Ne_cc + Ne_local);
            Qe_3(Qe_3 < 0 | Qe_3 > Qe_max) = NaN;
        end
        %------------------------------------------------------------------------
        function  invsig = invQsige(val)  % Inverse of excitatory sigmoid function; output in millivolts
            invsig = NaN*ones(size(val));
            ok = find(val > 0 & val < Qe_max);
            invsig(ok) = theta_e - log(Qe_max./val(ok) - 1)/Cse;
        end
        %------------------------------------------------------------------------
        function  invsig = invQsigi(val) % Inverse of inhibitory sigmoid function
            invsig = NaN*ones(size(val));
            ok = find(val > 0 & val < Qi_max);
            invsig(ok) = theta_i - log(Qi_max./val(ok) - 1)/Csi;
        end
        
    end
end



