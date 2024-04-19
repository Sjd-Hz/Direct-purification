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

r=0.45;
F_max=0.999;

%%%amplitude damping channel:
e0 = [1 0; 0 sqrt(1-r)];
e1 = [0 sqrt(r); 0 0];
E0=kron(Id2,e0);
E1=kron(Id2,e1);
rer=(E0*rsi*E0'+E1*rsi*E1');
% Measurement operators:
th=atan(1/sqrt(1-r));
M0=([cos(th),0;0,sin(th)]);
M1=([sin(th),0;0,cos(th)]);
%Two qubit measurement operators, only apply on the second qubit:
M20=kron(Id2,M0);
M21=kron(Id2,M1);
rM0=(M20*rsi*M20');
rM1=(x2*M21*rsi*M21'*x2');
rm=rM0+rM1;

rhofn0_det=(E0*rm*E0'+E1*rm*E1');

%%%%%%%%%%%  Entanglement distillation %%%%%%%%%%

F_0=trace(rhofn0_det*rsi);

kk=30; %set k round
fid_tot(1)=F_0;
fid_tot_det(1)=F_0;
%fid0_det(1)=F_0;
for ii=2:kk
   if(ii==2)
       rho_input{ii}=kron(rhofn0_det,rhofn0_det);
   else
       rho_input{ii}=kron(rho_0_det{ii-1},rho_0_det{ii-1});
   end
   %%%coincidence measurement result: 11
   rs=(PartialTrace(CZ1*CNOTed*rho_input{ii}*CNOTed'*CZ1',[3,4]));
   p1_det(ii)=trace(rs);
   rho_1_det=(rs/trace(rs));
   fid1_det(ii)=trace(rho_1_det*rsi);
   %%%coincidence measurement result: 00
   rf=(PartialTrace(CZ0*CNOTed*rho_input{ii}*CNOTed'*CZ0',[3,4]));
   p0_det(ii)=trace(rf);
   rho_0_det{ii}=(rf/ p0_det(ii));   
   fid0_det(ii)=trace(rho_0_det{ii}*rmax);
   
   %%% Finding the probability of '11' by considering previous operations
   %%% probabilities
   mult=1;
   kd=ii-1;
   for dd=2:ii-1
       kd=kd-1;
   mult=mult*(p0_det(dd)^(2^(kd)));
   end
   p1r_det(ii)=mult*p1_det(ii);
   p0r_det(ii)=mult*p0_det(ii);
   fid_tot_det(ii)=(fid0_det(ii)* p0_det(ii)+fid1_det(ii)* p1_det(ii))/( p0_det(ii)+ p1_det(ii)); %%%total fidelity at each round
   if fid0_det(ii)>F_max %The condition to stop the entnaglement distillation process
     m=ii;
     break;
   end
end
figure(1)
L1=plot( fid0_det,'b-->','LineWidth',2); hold on
L2=plot( fid1_det,'->','Color',[0.4660 0.6740 0.1880],'LineWidth',2); hold on
axis tight
xlim([2 Inf])

grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Fidelity','FontSize',12);
legend({'Fid_0^{k}','Fid_1^{k}'},'Location','best');


figure(2)
L11=plot(p1r_det,'b-->','LineWidth',2); hold on
L22=plot(p0r_det,'-.>','Color',[0.4660 0.6740 0.1880],'LineWidth',2); hold on
legend({'p_{1}^{k}','p_{0}^{k}'},'Location','best');
axis tight
xlim([2 Inf])
grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Net probability','FontSize',12);

