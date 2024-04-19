clear;
close all;
clc
%%%Identity operators
Id=[1,0;0,1];
Id2=kron(Id,Id);
Id3=kron(Id2,Id);
Id4=kron(Id3,Id);
Id5=kron(Id4,Id);
Id6=kron(Id5,Id);
Id7=kron(Id6,Id);

%%% Unilateral x rotation:
rx=[0,1;1,0];
rx2=kron(rx,Id);
rx3=kron(rx2,Id);
rx4=kron(rx3,Id);
rx5=kron(rx4,Id);

%%% X rotation on all qubits:
x2=kron(rx,rx);
x3=kron(x2,rx);
x4=kron(x3,rx);
x5=kron(x4,rx);

%%% Defining initial maximally GHZ state
h0=[1,0]';
h1=[0,1]';

z2=kron(h0,h0);
z3=kron(z2,h0);
z4=kron(z3,h0);
z5=kron(z4,h0);

v2=kron(h1,h1);
v3=kron(v2,h1);
v4=kron(v3,h1);
v5=kron(v4,h1);
al=1/sqrt(2);
be=1/sqrt(2);
si5=(al)*z5+(be)*v5;
rsi5=si5*si5';
rsx5=rx5*rsi5*rx5';

%%%%Distillation operations:
sw=[1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
%%% Swapping operations are employed to rearrange the target and sacrificial entangled states so that they are positioned next to each other. 
o1=kron(Id4,kron(sw,Id4));
o2=kron(Id3,kron(sw,Id5));
o3=kron(Id2,kron(sw,Id6));
o4=kron(Id,kron(sw,Id7));
o5=kron(Id5,kron(sw,Id3));
o6=kron(Id4,kron(sw,Id4));
o7=kron(Id3,kron(sw,Id5));
o8=kron(Id6,kron(sw,Id2));
o9=kron(Id5,kron(sw,Id3));
o10=kron(Id7,kron(sw,Id));
Op=o10*o9*o8*o7*o6*o5*o4*o3*o2*o1;

CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0]; 
Cnot5=kron(CNOT,kron(CNOT,(kron(CNOT,kron(CNOT,CNOT)))));

%%%Sacrificial meausrement '0':
Mz0=[1,0;0,0];
Mz20=kron(Id,Mz0);
CZ0=kron(Mz20,kron(Mz20,kron(Mz20,kron(Mz20,Mz20))));

%%%Sacrificial meausrement '1':
Mz1=[0,0;0,1];
Mz21=kron(Id,Mz1);
CZ1=kron(Mz21,kron(Mz21,kron(Mz21,kron(Mz21,Mz21))));

jj=0;
for r=0:0.01:1
    jj=jj+1;
    %%% Amplitude damping:
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
    
    e5=cell(1,32);
    k=1;
    for j=1:16
    e5{k}=kron(e4{j},e0);
    e5{k+1}=kron(e4{j},e1);
    k=k+2;
    end

    %%% The damped state after passing through ADC:
    Rer=0;
     for ii=1:1:32
            R=e5{ii}*rsi5*e5{ii}';
            Rer=R+Rer;      
     end
    fider(jj)=trace(Rer*rsi5);
    %%% Weak measurement operations:
    N=5;
    t=atan(1/(1-r)^((N/2)-1));
    M1=[cos(t),0;0,sin(t)];
    M0=[sin(t),0;0,cos(t)];
    %%% Unilateral weak measurement to only apply on the first qubit
    M50=kron(M0,Id4);
    M51=kron(M1,Id4);
    
    rin0=M50*rsx5*M50';
    rin1=x5*M51*rsx5*M51'*x5';
    %%% Passing through ADC:
    Rf0=0;
    Rf1=0;
    for bb=1:1:32
        R0=e5{bb}*rin0*e5{bb}';
        Rf0=R0+Rf0;
        R1=e5{bb}*rin1*e5{bb}';
        Rf1=R1+Rf1;
    end
    rff=(Rf0+Rf1);         
    rfx=(rx5*rff*rx5');    
    fidb=trace(rfx*rsi5);
    %%% Distillation process:
    rin0=kron(rfx,rfx);    
    rf=Op*rin0*Op';
    rho_input=rf;
    %%%sacrificial measurement result: 1
    rs=(PartialTrace(CZ1*Cnot5*rho_input*Cnot5'*CZ1',[2,4,6,8,10]));
    p1(jj)=trace(rs);
    rho_1=(rs/p1(jj));       
    fid1(jj)=trace(rho_1*rsi5);
    %%%sacrificial measurement result: 0
    rf=(PartialTrace(CZ0*Cnot5*rho_input*Cnot5'*CZ0',[2,4,6,8,10]));
    p0=trace(rf);
    rho_0=(rf/p0);
    fid0(jj)=trace(rho_0*rsi5);
    %%%In case of measuremet result 0, we only keep the result if fidelity is improved compared to the damped fideity:
    if fid0(jj)>=fider(jj)
       fid0_5(jj)=fid0(jj);
       p0_5(jj)=p0;
       p_tot(jj)=p0_5(jj)+p1(jj);
    else
       p_tot(jj)=p1(jj);
       fid0_5(jj)=NaN;
       p0_5(jj)=NaN;
    end  
end
figure(1)
L1=plot( fid0_5,'b-','LineWidth',2); hold on
L2=plot(fid1,'m-','LineWidth',3); hold on
L3=plot(fider,'k-','LineWidth',3); hold on
axis tight
grid on
xlabel('r','FontSize',12);
ylabel('Fidelity','FontSize',12);
xticks([10 20 30 40 50 60 70 80 90 100])
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})

figure(2)
L11=plot(p0_5,'b-','LineWidth',2); hold on
L22=plot(p1,'m-','LineWidth',3); hold on
L33=plot(p_tot,'g-','LineWidth',3); hold on
axis tight
grid on
xlabel('r','FontSize',12);
ylabel('Probability','FontSize',12);
xticks([10 20 30 40 50 60 70 80 90 100])
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})

