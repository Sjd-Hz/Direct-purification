clear,clc
close all;

H = [1;0]; %|0>
V = [0;1]; %|1>
zero=H;
Id2 = eye(2,2);
rx=[0,1;1,0]; %Pauli X
x2=kron(rx,rx);
si=1/sqrt(2)*(kron(H,V)+kron(V,H));
rsi=si*si';
rmax=rsi;
%%%define operators for entanglement distillation process:
CNOTed = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];

Mz0=[1,0;0,0];
CZ0=kron(Id2,kron(Id2,kron(Mz0,Mz0))); %Measuring the second pair in state 0

Mz1=[0,0;0,1];
CZ1=kron(Id2,kron(Id2,kron(Mz1,Mz1)));%Measuring the second pair in state 1

for r_idx = 1:100
    r = r_idx * 0.01; 

    %%%%%%%%%%%%%%%%%%%%%%%%% amplitude damping channel:
    e0 = [1 0; 0 sqrt(1-r)];
    e1 = [0 sqrt(r); 0 0];
    E0=kron(Id2,e0);
    E1=kron(Id2,e1);
    rer=(E0*rsi*E0'+E1*rsi*E1');
    fid_AD(r_idx)=trace(rer*rsi);
    % Measurement operators:
    th=atan(1/sqrt(1-r));
    M0=([cos(th),0;0,sin(th)]);
    M1=([sin(th),0;0,cos(th)]);
    %Two qubit measurement operators, only apply on the second qubit:
    M20=kron(Id2,M0);
    M21=kron(Id2,M1);

    %%%%%%%%%%%%%%%%%%%%%%%%% pre-distillation process
    rM0=(M20*rsi*M20');
    rM1=(x2*M21*rsi*M21'*x2');
    rm=rM0+rM1;   
    rhofn0=(E0*rm*E0'+E1*rm*E1');
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Entanglement distillation process
    rho_input=kron(rhofn0,rhofn0);

    %%%sacrificial measurement result: 11
    rsd=(PartialTrace(CZ1*CNOTed*rho_input*CNOTed'*CZ1',[3,4]));
    p1(r_idx)=trace(rsd);  
    rho_1=(rsd/trace(rsd));
    fid1(r_idx)=trace(rho_1*rsi);

    %%%sacrificial measurement result: 00
    rfd=(PartialTrace(CZ0*CNOTed*rho_input*CNOTed'*CZ0',[3,4]));
    p0(r_idx)=trace(rfd);
    rho_0=(rfd/ p0(r_idx));
    fid0(r_idx)=trace(rho_0*rmax);
    %%%In case of measuremet result 00, we only keep the result if fidelity is improved compared to the damped fideity:
    if fid0(r_idx)>=fid_AD(r_idx)
        fid0(r_idx)=fid0(r_idx);
        p0_det(r_idx)=p0(r_idx);
        p_tot(r_idx)=p0_det(r_idx)+p1(r_idx);
    else
        p_tot(r_idx)=p1(r_idx);
        fid0(r_idx)=NaN;
        p0_det(r_idx)=NaN;
    end
  
    F0=2*(1-r)/(2-r);
    p1(r_idx)=1/2*(F0^2);
    if r<0.45
        ptot(r_idx)=F0^2+(1-F0)^2;
    else
        ptot(r_idx)=1/2*(F0^2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Weak measurement and its reversal protection scheme (WMRPS):  
    fid_wmr = @(x) -(x(1) + x(2) + r - x(1)*r - 2*(1 - x(1))^(1/2)*(1 - x(2))^(1/2)*(1 - r)^(1/2) - 2)/((2*x(2)*r - 2)*(x(1) - 1) - 2*x(2) + 2);
    p_wmr = @(x) ((x(1) - 1)*(r - 1))/2 - x(2)/2 + (r*(x(1) - 1)*(x(2) - 1))/2 + 1/2;
    %%% x is the weak measurement strength    
    x1_values = 0:0.01:0.99; %pre-WM strength
    x2_values = 0:0.01:0.99; %post-WM strength
    optimal_val_tmp = -Inf; 
    optimal_params_tmp = zeros(2, 1);
        
    %%%Finding the maximum fidelity of WMRPS:
    for x1_idx = 1:length(x1_values)
            for x2_idx = 1:length(x2_values)
                x = [x1_values(x1_idx), x2_values(x2_idx)];
                value = fid_wmr(x);
                prob = p_wmr(x);            
                % Update the temporary optimal value and parameters if necessary
                if value > optimal_val_tmp 
                    optimal_val_tmp = value;
                    optimal_params_tmp = x;
                    optimal_p_tmp = prob;
                end
            end
    end    
    Fid_WMR(r_idx) = optimal_val_tmp;
    optimal_params(:, r_idx) = optimal_params_tmp;
    p_WMR(r_idx) = optimal_p_tmp;
end

figure(1)
tf=0.01:0.01:1;
L1=plot( tf,fid0,'b-.','LineWidth',2); hold on
L2=plot(tf,fid1,'r-','LineWidth',2); hold on
L3=plot(tf,fid_AD,'k--','LineWidth',2); hold on
L4= plot(tf, Fid_WMR, '-','color','[0.9290 0.6940 0.1250]', 'LineWidth', 2);hold on
axis tight
grid on
xlabel('r','FontSize',12);
ylabel('Fidelity','FontSize',12);
legend({'Fid_0','Fid_1','Fid_{AD}','Fid_{WMRPS}'},'Location','southwest');


figure(2)
L11=plot(tf,p0_det,'b-.','LineWidth',2); hold on
L22=plot(tf,p1,'r-','LineWidth',2); hold on
L33=plot(tf,p_tot,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2); hold on
L44 = plot(tf, p_WMR, '-','color','[0.9290 0.6940 0.1250]', 'LineWidth',2);

legend({'p_0','p_1','p_{tot}','p_{WMRPS}'},'Location','northeast');
axis tight
grid on
xlabel('r','FontSize',12);
ylabel('Probability','FontSize',12);

