clear;
%close all;
clc

r=0.25;
f_max=0.99;

Id=[1,0;0,1];
Id2=kron(Id,Id);
Id3=kron(Id2,Id);
Id4=kron(Id3,Id);
Id5=kron(Id4,Id);
sw=[1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];

rx=[0,1;1,0];
%Unilateral Pauli-X:
rx2=kron(rx,Id);
rx3=kron(rx,kron(Id,Id));
rx4=kron(rx3,Id);

%Applying Pauli-X on all qubits
x2=kron(rx,rx);
x3=kron(x2,rx);
x4=kron(x3,rx);

h0=[1,0]';
h1=[0,1]';
z2=kron(h0,h0);
z3=kron(z2,h0);
z4=kron(z3,h0);
v2=kron(h1,h1);
v3=kron(v2,h1);
v4=kron(v3,h1);

si4=(1/sqrt(2))*(z4+v4);
rsi4=si4*si4';% The initial maximally entangled GHZ state
rsx4=rx4*rsi4*rx4';%Changing the structure of the denisty matrix

%Amplitude damping channel:
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

e4=cell(1,16);
k=1;
for j=1:8
e4{k}=kron(e3{j},e0);
e4{k+1}=kron(e3{j},e1);
k=k+2;
end

N=4;
t=atan(1/(1-r)^((N/2)-1));
M0=[cos(t),0;0,sin(t)];
M1=[sin(t),0;0,cos(t)];

M40=kron(Id,kron(Id,kron(Id,M0)));
M41=kron(Id,kron(Id,kron(Id,M1)));

rm=M40*rsx4*M40'+x4*M41*rsx4*M41'*x4';
  rff=0;
    for bb=1:1:16
        Rf=e4{bb}*rm*e4{bb}';
        rff=rff+Rf;    
    end
rfx=(rx4*rff*rx4');
%trace(rfx)
fid0(1)=trace(rfx*rsi4);
rin0=kron(rfx,rfx);

%%%define operators for entanglement distillation process:
%%Swap the qubits to bring sacrificial and target qubits after each other
o1=kron(Id3,kron(sw,Id3));
o2=kron(Id2,kron(sw,Id4));
o3=kron(Id,kron(sw,Id5));
o4=kron(Id4,kron(sw,Id2));
o5=kron(Id3,kron(sw,Id3));
o6=kron(Id5,kron(sw,Id));
Op=o6*o5*o4*o3*o2*o1;
rf=Op*rin0*Op';

CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0]; 
Cnot4=kron(CNOT,(kron(CNOT,kron(CNOT,CNOT))));
rcnot=Cnot4*rf*Cnot4;


Mz0=[1,0;0,0];
Mz20=kron(Id,Mz0);
CZ0=kron(Mz20,kron(Mz20,kron(Mz20,Mz20)));


Mz1=[0,0;0,1];
Mz21=kron(Id,Mz1);
CZ1=kron(Mz21,kron(Mz21,kron(Mz21,Mz21)));

k=20; %set k round
rho_input=cell(1,k);
rho_success=cell(1,k);
rho_fail=cell(1,k);
yield_n=cell(1,k);
tot_effi=cell(1,k);
 p_fail=cell(1,k);
sum=0;
for ii=2:k
   
   if(ii==2)
       rho_input{ii}=rf;
   else
       rho_per=kron(rho_0{ii-1},rho_0{ii-1});  
       rho_input{ii}=Op*rho_per*Op';       
   end
  rs=(PartialTrace(CZ1*Cnot4*rho_input{ii}*Cnot4'*CZ1',[2,4,6,8]));
   p1(ii)=trace(rs);
   rho_1{ii}=(rs/p1(ii));
    %disp( rho_1{ii})
   fid1{ii}=trace(rho_1{ii}*rsi4);
   %%%00
   rf=(PartialTrace(CZ0*Cnot4*rho_input{ii}*Cnot4'*CZ0',[2,4,6,8]));
   p0(ii)=trace(rf);
   rho_0{ii}=(rf/p0(ii));
  %disp( rho_0{ii})
   fid0(ii)=trace(rho_0{ii}*rsi4);
   
    %%% Finding the probability of '11' and '00' by considering previous operations probabilities
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
title(['GHZ^4 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);
legend({'Fid_{0}^{k}'});
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})

figure(2)
L11=plot(p1r,'m--o','LineWidth',2); hold on
L22=plot(p0r,'-.o','Color',[0.9290 0.6940 0.1250],'LineWidth',2); hold on
legend({'P_{1}^{k}','P_{0}^{k}'});
xlim([2 Inf])
grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Net probability','FontSize',12);
title(['GHZ^4 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})