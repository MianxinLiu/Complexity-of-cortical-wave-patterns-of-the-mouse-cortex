function Ve_rec=anesthesia_fieldmodel(p,mu,LRposition)
%% field model simulation 
% input parameters: mu is the long range strength
% p is the anesthetic degree
% structure LRposition has two subfields: 'linkout' and 'receive'
% both LRposition.linkout and LRposition.receive are N x 2 matrices 
% they are the outgoing and incoming position of the N long-range
% connections: the j-th projection from LRposition.linkout(j,:) to 
% LRposition.receive(j,:)



%% model parameters
N = 100;		%  the size of network is N*N
Lx=10;			% ength of the region(mm)
I_sc=0.018; % subcortical firing   strength (/ms)


[Tau_e,Tau_i] = deal(40);		% E,I neuron time constant (ms)
[ge,gi]=deal(0.156,0.875);	% E,I synaptic gain strengths , ACh,GABA effect (mV*ms)
tau_de = 5;		% EPSP decay time (ms)
tau_di = 20; % IPSP decay time  (ms)
[Ve_rev, Vi_rev,Ve_rest,Vi_rest]=deal(0,-70,-62.5,-64);% reversal,resting potential for E,I  (mV)
[Ne_cc,Ne_local,Ni_local]=deal(200,85,120);%column-column E;intracolumn E,I;
[Qe_max,Qi_max,theta_e,theta_i,sigma_e,sigma_i]=deal(0.03,0.06,-58.5,-58.5,3,5);% sigmoid maximum firing rates (/ms),threshold (mV), 'width' (mV)
v = 0.056;				% axonal conduction velocity (mm/ms)
r = 0.1;			% characteristic length scale for axons connectivity (mm)

a = 5.2;		% noise scale-factor a, if noise, then a=5.2
% a= 0;   % if no noise, a=0


delta_lambda_i=0.08; % contribution factor
delta_lambda_e=0.01; % contribution factor
lambda_i=1+delta_lambda_i*p; % scale factor
lambda_e=1-delta_lambda_e*p; % scale factor
delta_D=3/7; % contribution factor
lambda_D=1-delta_D*p; % scale factor
Di=0.07*lambda_D;    % I diffusion strength, already include the anesthesia effect
De = Di*0.01;		% E diffusion strength, already include the anesthesia effect

linkout_site=LRposition.linkout;
receive_site=LRposition.receive;

if size(linkout_site,1)~=size(receive_site,1)
    error('long-range position error');
end
is_it_longrange=~isempty(linkout_site); % 1 if include a long-range connection

if is_it_longrange
    LR_receive_position=receive_site(:,1)+(receive_site(:,2)-1)*N; % transform the idex
    LR_out_position=linkout_site(:,1)+(linkout_site(:,2)-1)*N; % transform the idex
    ran_id=randperm(length(LR_out_position));
    LR_out_position=LR_out_position(ran_id); % randomly permute
    M=size(LR_out_position,1); % total long-range links number
    distance_in_grid=sqrt((linkout_site(:,1)-receive_site(:,1)).^2+(linkout_site(:,2)-receive_site(:,2)).^2); % the grid distance between positions
end
%% technical parameters
rand_state = floor(sum(100*clock)); rng(rand_state);% random seed to initialize random number generator

branch=1; % branch == 1 ==> bottom branch

dx =Lx/N;			% spatial  step resolution (mm)
T_running = 3*10^3;		% length of simulation (ms)
dt = 0.4;			% simulation time step (ms) 
Nsteps = T_running/dt; % number of time-steps for simulation
dt_max1 = 0.25*Tau_i*dx^2/(Di+eps);
dt_max2 = dx/(sqrt(2)*v);% + eps for case D2 = 0
if dt > min(dt_max1,dt_max2) % Courant conditions: 1 for inhibitory diffusion, 2 for wave equation
    error('dt too large! , lager than  %f ',min(dt_max1,dt_max2));
end
t_discard = 0*10^3; % the initial discarding time (ms)

if is_it_longrange
    delay_step = round(Lx*distance_in_grid/(v*N*dt)); % steps of the delay ,dx * separation is the absolute separete distance, hypot=sqrt( ^2+ ^2)
end

Recordtime= 10; % record time step (ms)
Recordstep = Recordtime/dt;	% per screen update

%% auxiliary parameters
z=1:N;zm1=[N,1:N-1];zp1=[2:N,1];xm1y=[];xp1y=[];xym1=[];xyp1=[];
for j=1:N
    for i=1:N
        xm1y=[xm1y;zm1(i)+(z(j)-1)*N];
        xp1y=[xp1y;zp1(i)+(z(j)-1)*N];
        xym1=[xym1;z(i)+(zm1(j)-1)*N];
        xyp1=[xyp1;z(i)+(zp1(j)-1)*N];
    end
end % auxiliary matrix for Laplacian operator      
B_noise = a*sqrt(I_sc/dt); % noise-amplitude coefficients for subcortical flux (note 1/sqrt(dt) factor)
W1=(v/r*dt)^2;W2=2-2*(v/r*dt)-(v/r*dt)^2;W3=2*(v/r*dt)-1;W4=(v*dt/dx)^2; % wave equation coefs
Ke=1 - 2*dt/tau_de;Ki=1 - 2*dt/(tau_di*lambda_i);
ge2t=(1/tau_de)^2*dt;gi2t=dt/(tau_di*lambda_i)^2;
D1x=De/dx^2;D2x=Di/dx^2;
Cse=pi/(sqrt(3) * sigma_e);Csi=pi/(sqrt(3) * sigma_i);

%% set initial conditions

%  First find locate steady state(s)
[~, ~, Ve0, Vi0] = SO_ss_finder();
num_roots = length(Ve0);
if num_roots > 1
    [Ve0, Vi0] = deal( Ve0(branch), Vi0(branch)); % since they are ordered increasingly, 1,3 is the low,top branch
end % initialize the grids around the homogeneous steady-state values

[Ve,Vi]= deal(ones(N^2,1)*Ve0+rand(N^2,1),ones(N^2,1)*Vi0+rand(N^2,1));

% set special initial activation (if any)
% [Ve,Vi]= deal(ones(N^2,1)*Ve0,ones(N^2,1)*Vi0); % give initial deviation to exclude steady
% Ve_initial=reshape(Ve,[100,100]); % initially activate a region
% Ve_initial([66,67],[34,35])=-54;
% Ve_initial([66,67],[34,35])=-60;
% Ve=reshape(Ve_initial,[10000,1]);

Qe = Qe_max./(1 + exp((theta_e-Ve)*Cse));
Qi = Qi_max./(1 + exp((theta_i-Vi)*Csi));
[phi_e,phi_e_old] = deal(Qe);
Phi_e = Ne_cc*phi_e + Ne_local*Qe + I_sc;
Phi_i =                  Ni_local*Qi;
[F_e,F_i,phi_e_longrange]  = deal(zeros(N^2,1));



if is_it_longrange
    Qe_LR=zeros(M,Nsteps); % record the rate of outgoing sites
    % totally there are M sites/links
end

% some recording array
Ve_rec=zeros(N^2,T_running/Recordtime);
kr=1;



%% computing
for k = 1: Nsteps

    % 1. update wave equations
    TEMP=W1*Qe+W2*phi_e+W3*phi_e_old+W4*(phi_e(xm1y)+phi_e(xp1y)+phi_e(xym1)+phi_e(xyp1)-4*phi_e);
    phi_e_old = phi_e; 
    phi_e = TEMP;
    % 2. update the  synaptic flux equations (include pt fluxes and sc noise)
    if is_it_longrange
       % Qe_LR(k)=Qe(LR_out_position); 
        for s=1:M
            Qe_LR(s,k)=Qe(LR_out_position(s));
        end
        for s=1:M 
            if k > delay_step(s) %%% two-point connections-- only active for t > t_delay
                phi_e_longrange(LR_receive_position(s)) = mu * Qe_LR(s,k-delay_step(s));
            end
        end
    end
    %%% here F= dPhi/dt to avoid the second-order time derivative
    F_e=F_e*Ke+ge2t*(Ne_cc*phi_e+Ne_local*Qe+I_sc+phi_e_longrange+B_noise*randn(N^2,1)-Phi_e); 
    F_i=F_i*Ki+gi2t*(Ni_local*Qi-Phi_i);
    Phi_e=Phi_e+F_e*dt;
    Phi_i=Phi_i+F_i*dt;
   
    
    % 3. update the soma voltages
    Ve = Ve + dt/Tau_e*(Ve_rest - Ve + ...
        lambda_e*ge*( Ve_rev - Ve ).*Phi_e + lambda_i*gi*( Vi_rev - Ve ).*Phi_i + ...
        D1x*(Ve(xm1y)+Ve(xp1y)+Ve(xym1)+Ve(xyp1)-4*Ve));
    Vi = Vi + dt/Tau_i*(Vi_rest - Vi + ...
        lambda_e*ge*( Ve_rev - Vi  ).*Phi_e + lambda_i*gi*( Vi_rev - Vi ).*Phi_i + ...
        D2x*(Vi(xm1y)+Vi(xp1y)+Vi(xym1)+Vi(xyp1)-4*Vi));

    % 4. update the firing rates
    Qe = Qe_max./(1 + exp((theta_e-Ve)*Cse));
    Qi = Qi_max./(1 + exp((theta_i-Vi)*Csi));
    
    if mod(k, Recordstep) == 0 % record   
        Ve_rec(:,kr)=Ve;
        kr=kr+1;
    end
end


Ve_rec(:,1:t_discard/Recordtime)=[];
Ve_rec=reshape(Ve_rec,[100,100,T_running/Recordtime]);

    function [Qe_root, Qi_root, Ve_root, Vi_root] = SO_ss_finder()
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
        [Qe_root,Qi_root,Ve_root,Vi_root] = deal(ones(n_roots, 1)); % find each roots
        for j = 1: n_roots
            fun = @ (Qe_1) Qe_resid(Qe_1);
            Qe_root(j) = fzero(fun, brack_Qe(j,:));
            Ve_root(j) = invQsige(Qe_root(j));
            Qi_root(j) = EqA(Qe_root(j), Ve_root(j));
            Vi_root(j) = invQsigi(Qi_root(j));
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






