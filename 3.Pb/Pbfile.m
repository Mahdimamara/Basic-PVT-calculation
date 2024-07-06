clc;clear;
% EOS=input('choose the EOS function:\n SRK(1),SRK&G(2),PR76(3),PR78(4)\n');
% fluidtype=input('choose the fluid type:\n Condensate(1), Black oil(2),Volatile oil(3)\n');
%%
disp('buuble points diagram');
disp('Equation of state: SRK');
disp('Fluid type: Gas condensate');
EOS=1;
fluidtype=1;
[Tc, Pc, W, Mw, comp, K_ij, Nc, Mw_t, R] = FluidPT(fluidtype);

Fv=0;

[m,C,sigma1,sigma2,omegaa,omegab] = EQchar(W,EOS,Mw,Nc);
%%
countb=1;
Errorb=2;
Pb=zeros(1,500);
Tb=zeros(1,500);
Kb=zeros(500,Nc);
Tb(countb)=380;

count=1;

for i=1:35
    %%
        Pb(countb)=sum(comp.*Pc.*exp(5.37*(1+W).*(1-Tc/Tb(countb))))+800000.*(Tb(countb)-400)+20*10^6;
        Kb(countb,:)=(Pc./Pb(countb)).*exp(5.37*(1+W).*(1-Tc/Tb(countb)));
 
 
    %%
    [b,a,ac,alpha] = coglob(R,Tc,Pc,m,Tb(countb),omegaa,omegab);
    xb=zeros(1,Nc);
    yb=zeros(1,Nc);
    Tolb=1;
    %%
    while Tolb>10^-3
        xb=comp./(1+Fv*(Kb(countb,:)-1));
        yb=Kb(countb,:).*xb;
        [Sx,atx,btx,Ax,Bx] = coefficientcal(xb,Pb(countb),Tb(countb),Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
        [Sy,aty,bty,Ay,By] = coefficientcal(yb,Pb(countb),Tb(countb),Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
        [Zx] = solvecub(Ax,Bx,C,b,btx,Sx,atx,sigma1,sigma2,xb,R,Tb(countb));
        [Zy] = solvecub(Ay,By,C,b,bty,Sy,aty,sigma1,sigma2,yb,R,Tb(countb));
        [Phix] = FugacityCal(Zx,Ax,Bx,sigma1,sigma2,b,btx,Sx,atx);
        [Phiy] = FugacityCal(Zy,Ay,By,sigma1,sigma2,b,bty,Sy,aty);
        Kb(countb,:)=Phix./Phiy;
        Db=sum(comp.*(Kb(countb,:)-1)./(1+Fv*(Kb(countb,:)-1)));
        D=sum(comp./(Kb(countb,:)));
        Tolb=abs(Db);
        
        if isreal(Db)==0
            break
        end
        if Tolb<10^-5
            break
        else
            Pb(countb)=(D)*Pb(countb);
            
        end
    end
    %%
    
    a=2;
    if fluidtype==3
        a=a*0.5;
    end
    if Fv~=0
        a=a*10*Fv;
    end
    Errorb=sum(abs(Kb(countb,:)-1))<(a);
    countb=countb+1;
    Tb(countb)=Tb(countb-1)+5;
    
    
end
%%
Tb=Tb(Tb>0);
Pb=Pb(Pb>0);
for i=2:size(Tb)
    if abs(Pb(i)-Pb(i-1))>1000 && abs(Tb(i)-Tb(i-1))<50
        Pb(i)=Pb(i-1);
        Tb(i)=Tb(i-1);
    end
end
%%
T_b=Tb(1:end-2);
P_b=Pb(1:end-1)./10^5;
hold on
grid on

if Fv==0
    plot(T_b,P_b,'r*');
elseif Fv==0.1
    plot(T_b,P_b,'y.');
elseif Fv==0.2
    plot(T_b,P_b,'m.');
elseif Fv==0.5
    plot(T_b(1:end-4),P_b(1:end-4),'c.');
else
    plot(T_b,P_b);
end
xlabel('Temrature K');
ylabel('Pressure psia');


