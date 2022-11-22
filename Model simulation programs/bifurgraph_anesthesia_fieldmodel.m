function [EIR_curve_up,EIR_curve_middle,EIR_curve_down]=bifurgraph_anesthesia_fieldmodel()
% plot the bifurcation curve for the steady-state without spatical couple


% parameter setting
I_sc=0.018; % subcortical firing strength (/ms)
[ge,gi]=deal(0.156,0.875);	% E,I synaptic gain strengths , ACh,GABA effect (mV*ms)
[Ve_rev, Vi_rev,Ve_rest, Vi_rest]=deal(0,-70,-62.5,-64);% reversal,resting potential for E,I  (mV)
[Ne_cc,Ne_local,Ni_local]=deal(200,85,120);%column-column E;intracolumn E,I;
[Qe_max,Qi_max,theta_e,theta_i,sigma_e,sigma_i]=deal(0.03,0.06,-58.5,-58.5,3,5);% sigmoid maximum firing rates (/ms),threshold (mV), 'width' (mV)
delta_lambda_i=0.08; % contribution factor
delta_lambda_e=0.01; % contribution factor
% auxiliary parameters
Cse=pi/(sqrt(3) * sigma_e);Csi=pi/(sqrt(3) * sigma_i);

p=-1:0.002:1;
rootnums=zeros(1,length(p));
equilibriumpoints=cell(length(p),1);
EI_ratiopoints=cell(length(p),1);
for k=1:length(p)
    Ldi=1+p(k)*delta_lambda_i; % scale factor
    Lde=1-p(k)*delta_lambda_e; % scale factor
    [~,~,Ve,~,EIratio] = findroots(Ldi,Lde);
    equilibriumpoints{k}=Ve';
    rootnums(k)=length(equilibriumpoints{k});
    EI_ratiopoints{k}=EIratio';
end
c=1;
while length(equilibriumpoints{c})==1
    c=c+1;
end
bp1=c;
while length(equilibriumpoints{c})~=1
    c=c+1;
end
bp2=c-1;
equilibrium_curve_up=zeros(1,bp2);
equilibrium_curve_down=zeros(1,length(p)-bp1+1);
equilibrium_curve_middle=zeros(1,bp2-bp1+1);
EIR_curve_up=zeros(1,bp2);
EIR_curve_down=zeros(1,length(p)-bp1+1);
EIR_curve_middle=zeros(1,bp2-bp1+1);

for k=1:bp2
    equilibrium_curve_up(k)=equilibriumpoints{k}(end);
    EIR_curve_up(k)=EI_ratiopoints{k}(end);
end
for k=1:length(p)-bp1+1
    equilibrium_curve_down(k)=equilibriumpoints{bp1+k-1}(1);
    EIR_curve_down(k)=EI_ratiopoints{bp1+k-1}(1);
end
for k=1:bp2-bp1+1
    equilibrium_curve_middle(k)=equilibriumpoints{bp1+k-1}(2);
    EIR_curve_middle(k)=EI_ratiopoints{bp1+k-1}(2);
end

fig = figure(); clf;
set(fig, 'Position', [100 100 400 400]);
set(fig, 'Name', 'Bifurcation Graph of Ve');

yyaxis left;
plot(p(1:bp2),equilibrium_curve_up,'k-', 'linewidth',2); hold on;
plot(p(bp1:bp2),equilibrium_curve_middle,'k--', 'linewidth',2); hold on;
plot(p(bp1:end),equilibrium_curve_down,'k-', 'linewidth',2); hold on;
ylabel('Ve (mV)', 'fontsize',12);
xlabel('p', 'fontsize',12);
xlim([p(1),p(end)]);
ylim([-64 -55]);
set(gca,'Ycolor','b');

yyaxis right;
plot(p(1:bp2),EIR_curve_up,'k-', 'linewidth',2); hold on;
plot(p(bp1:bp2),EIR_curve_middle,'g--', 'linewidth',2); hold on;
plot(p(bp1:end),EIR_curve_down,'c-', 'linewidth',2); hold on;
ylabel('E-I current ratio', 'fontsize',16);
ylim([0.7,1.4]);
set(gca,'Ycolor','r');



    function [Phi_e,Phi_i,Ve,Vi,EIratio] = findroots(lambda_i,lambda_e)
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


