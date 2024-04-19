clear;
%close all;
clc

r=0.25;
f_max=0.99;

Id=[1,0;0,1];

rx=[0,1;1,0];
rz=[1,0;0,-1];
rx2=kron(rx,Id);
rx3=kron(rx2,Id);
rx4=kron(rx3,Id);
rx5=kron(rx4,Id);

x2=kron(rx,rx);
x3=kron(x2,rx);
x4=kron(x3,rx);
x5=kron(x4,rx);

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

si5=1/sqrt(2)*(z5+v5);
rsi5=si5*si5';
rsx5=rx5*rsi5*rx5';

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

N=5;
t=atan(1/(1-r)^((N/2)-1));
M1=[cos(t),0;0,sin(t)];
M0=[sin(t),0;0,cos(t)];


M20=kron(M0,Id);
M21=kron(M1,Id);

M30=kron(M20,Id);
M31=kron(M21,Id);

M40=kron(M30,Id);
M41=kron(M31,Id);

M50=kron(M40,Id);
M51=kron(M41,Id);

rin0=M50*rsx5*M50';
rin1=x5*M51*rsx5*M51'*x5';

    Rf0=0;
    Rf1=0;
    for bb=1:1:32
        R0=e5{bb}*rin0*e5{bb}';
        Rf0=R0+Rf0;
        R1=e5{bb}*rin1*e5{bb}';
        Rf1=R1+Rf1;
    end
     %assume((0<r)&(r<1))
     rff=(Rf0+Rf1);
     
rfx=(rx5*rff*rx5');
%trace(rfx)
fid0(1)=trace(rfx*rsi5);
%rin=kron(rsi3,rsi3);
rin=kron(rfx,rfx);


Id2=kron(Id,Id);
Id3=kron(Id2,Id);
Id4=kron(Id3,Id);
Id5=kron(Id4,Id);
Id6=kron(Id5,Id);
Id7=kron(Id6,Id);
sw=[1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];

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
rf=Op*rin*Op';

CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0]; 
Cnot5=kron(CNOT,kron(CNOT,(kron(CNOT,kron(CNOT,CNOT)))));

Mz0=[1,0;0,0];
Mz20=kron(Id,Mz0);
CZ0=kron(Mz20,kron(Mz20,kron(Mz20,kron(Mz20,Mz20))));

Mz1=[0,0;0,1];
Mz21=kron(Id,Mz1);
CZ1=kron(Mz21,kron(Mz21,kron(Mz21,kron(Mz21,Mz21))));


sum=0;
F=0;
N=5;
F(1)=(2*(1-r)^(N-1))/((1-r)^(N-2)+1)
      
for ii=2:20
   
   if(ii==2)
       rho_input{ii}=rf;
   else
       rho_per=kron(rho_0{ii-1},rho_0{ii-1});    
       rho_input{ii}=Op*rho_per*Op';    
   end
  rs=(PartialTrace(CZ1*Cnot5*rho_input{ii}*Cnot5'*CZ1',[2,4,6,8,10]));
   p1(ii)=trace(rs);
   rho_1{ii}=(rs/p1(ii));
    %disp( rho_1{ii})
   fid1{ii}=trace(rho_1{ii}*rsi5);
   %%%00
   rf=(PartialTrace(CZ0*Cnot5*rho_input{ii}*Cnot5'*CZ0',[2,4,6,8,10]));
   p0(ii)=trace(rf);
   rho_0{ii}=(rf/p0(ii));
  %disp( rho_0{ii})
   fid0(ii)=trace(rho_0{ii}*rsi5);
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
title(['GHZ^5 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);
legend({'Fid_{0}^{k}'});
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})


figure(2)
L11=plot(p1r,'k--o','LineWidth',2); hold on
L22=plot(p0r,'r-.o','LineWidth',2); hold on
legend({'P_{1}^{k}','P_{0}^{k}'});
xlim([2 Inf])
grid on
xlabel('k^{th} round of distillation','FontSize',12);
ylabel('Net probability','FontSize',12);
title(['GHZ^5 (r=' ,num2str(r),', ','F_{max}=' ,num2str(f_max),')']);
xticks([1 2 3 4])
xticklabels({'0','1','2','3'})
