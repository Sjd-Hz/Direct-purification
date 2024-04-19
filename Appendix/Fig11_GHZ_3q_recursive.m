clear;
close all;
clc
%%%Identity operators
Id=[1,0;0,1];
Id2=kron(Id,Id);
Id3=kron(Id2,Id);
%%% Unilateral x rotation:
rx=[0,1;1,0];
rx2=kron(rx,Id);
rx3=kron(rx2,Id);
%%% X rotation on all qubits:
x2=kron(rx,rx);
x3=kron(x2,rx);
%%% Defining initial maximally GHZ state
al=1/sqrt(2);
be=1/sqrt(2);
h0=[1,0]';
h1=[0,1]';

z2=kron(h0,h0);
z3=kron(z2,h0);
v2=kron(h1,h1);
v3=kron(v2,h1);

si3=(al)*z3+(be)*v3;
rsi3=si3*si3';
rsx3=rx3*rsi3*rx3';
%%%%Distillation operations:
sw=[1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
%%% Swapping operations are employed to rearrange the target and sacrificial entangled states so that they are positioned next to each other. 
o1=kron(Id2,kron(sw,Id2));
o2=kron(Id,kron(sw,Id3));
o3=kron(Id3,kron(sw,Id));
Op=o3*o2*o1;

CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0]; 
Cnot3=kron(CNOT,kron(CNOT,CNOT));
%%%Sacrificial meausrement '0':
Mz0=[1,0;0,0];
Mz20=kron(Id,Mz0);
CZ0=kron(Mz20,kron(Mz20,Mz20));
%%%Sacrificial meausrement '1':
Mz1=[0,0;0,1];
Mz21=kron(Id,Mz1);
CZ1=kron(Mz21,kron(Mz21,Mz21));

%Amplitude damping channel:
r=0.25;

e0 = [1 0; 0 sqrt(1-r)];
e1 = [0 sqrt(r); 0 0];

e2=cell(1,4);
e2{1}=kron(e0,e0);
e2{2}=kron(e0,e1);
e2{3}=kron(e1,e0);
e2{4}=kron(e1,e1);

e3=cell(1,8);
k=1;
for j=1:4
e3{k}=kron(e2{j},e0);
e3{k+1}=kron(e2{j},e1);
k=k+2;
end

Rer=0;
 for ii=1:1:8
        R=e3{ii}*rsx3*e3{ii}';
        Rer=R+Rer;      
 end
fider=trace(Rer*rsi3);

%%% Weak measurement operations:    
N=3;
t=atan(1/(1-r)^((N/2)-1));
M1=[cos(t),0;0,sin(t)];
M0=[sin(t),0;0,cos(t)];
%%% Unilateral weak measurement to only apply on the first qubit     
M30=kron(M0,Id2);
M31=kron(M1,Id2);
    
rin0=M30*rsx3*M30';
rin1=x3*M31*rsx3*M31'*x3';
%%% Passing through ADC:
Rf0=0;
Rf1=0;
for bb=1:1:8
    R0=e3{bb}*rin0*e3{bb}';
    Rf0=R0+Rf0;
    R1=e3{bb}*rin1*e3{bb}';
    Rf1=R1+Rf1;
end
rff=(Rf0+Rf1);
rfx=rx3*rff*rx3';
fid0(1)=(trace(rfx*rsi3));
in0=kron(rfx,rfx);    
rf=Op*in0*Op';
k=20; %set k round

rho_input=cell(1,k);
f_max=0.99;
sum=0;
for ii=2:k
   
   if(ii==2)
       rho_input{ii}=rf;
   else
     rho_per=kron(rho_0{ii-1},rho_0{ii-1});
     rho_input{ii}=Op*rho_per*Op';
       
   end
   %%%sacrificial measurement result: 1
   rs=(PartialTrace(CZ1*Cnot3*rho_input{ii}*Cnot3'*CZ1',[2,4,6]));
   p1(ii)=trace(rs);
   rho_1{ii}=(rs/p1(ii));
   fid1{ii}=trace(rho_1{ii}*rsi3);
   %%%sacrificial measurement result: 0
   rf=(PartialTrace(CZ0*Cnot3*rho_input{ii}*Cnot3'*CZ0',[2,4,6]));
   p0(ii)=trace(rf);
   rho_0{ii}=(rf/p0(ii));
   fid0(ii)=trace(rho_0{ii}*rsi3);
   
   %%% Finding the probability of '1' and '0' by considering previous operations probabilities
   mult=1;
   kd=ii-1;
   for dd=2:ii-1
       kd=kd-1;
       mult=mult*(p0(dd)^(2^(kd)));
   end
   p1r(ii)=mult*p1(ii);
   p0r(ii)=mult*p0(ii);
   if fid0(ii)>f_max %The condition to stop the entnaglement distillation process
     m=ii;
     break;
   end
end

figure(1)
L1=plot( fid0,'m--o','LineWidth',2); hold on
grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Fidelity','FontSize',12);
legend({'Fid_{0}^{k}'});
title(['GHZ^3 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})

figure(2)
L11=plot(p1r,'m--o','LineWidth',2); hold on
L22=plot(p0r,'-.o','Color',[0.9290 0.6940 0.1250],'LineWidth',2); hold on
legend({'P_{1}^{k}','P_{0}^{k}'});
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})
xlim([2 Inf])
grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Net probability','FontSize',12);
title(['GHZ^3 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);