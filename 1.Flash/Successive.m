
clc;clear all;
% EOS=input('choose the EOS function:\n SRK(1),SRK&G(2),PR(3),\n');
% fluidtype=input('choose the fluid type:\n Condensate(1), Black oil(2),Volatile oil(3)\n');

%%
EOS=1;
fluidtype=1;
disp('An example of flash resaults');
disp('Equation of state: SRK');
disp('Fluid type: Gas condensate');
disp('Pressure:3615 psia');
disp('Temprature: 305 F');
[P, T, Tc, Pc, W, Mw, comp, K_ij, Nc, Mw_t, R,Parachor] = Fluid(fluidtype);
[m,C,sigma1,sigma2,omegaa,omegab] = EQchar(W,EOS,Mw,Nc);
%%
[K,u] = stability(P,T,comp,Tc,Pc,W,R,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,fluidtype,EOS);

if sum(u)<1
    disp('the phase is stable')
    break
end
%%
Tol=1;
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab);
if fluidtype==1
    Fv=1;
else
    Fv=0.5;
end
count=0;
%%
while Tol>10^-13
    F=1;
    while F>10^-14
        F=sum((K-1).*comp./(1+Fv.*(K-1)));
        F_prime=sum((-(K-1).^2).*comp./(1+Fv.*(K-1)).^2);
        Fv=Fv-F./F_prime;
    end
    x=comp./(1+Fv.*(K-1));
    y=K.*x;
    [Sx,atx,btx,Ax,Bx] = coefficientcal(x,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zx] = solvecub(Ax,Bx,C,b,btx,Sx,atx,sigma1,sigma2,x,R,T);
    [Phix] = FugacityCal(Zx,Ax,Bx,sigma1,sigma2,b,btx,Sx,atx);
    [Sy,aty,bty,Ay,By] = coefficientcal(y,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zy] = solvecub(Ay,By,C,b,bty,Sy,aty,sigma1,sigma2,y,R,T);
    [Phiy] = FugacityCal(Zy,Ay,By,sigma1,sigma2,b,bty,Sy,aty);
    K=Phix./Phiy;
    fl=Phix.*x*P;
    fv=Phiy.*y*P;
    Tol=sum(abs(1-fl./fv));
    count=count+1;
end
%%
Mwl=sum(x.*Mw);
Mwv=sum(y.*Mw);
Mwz=sum(comp.*Mw);
liquid_density=(P.*Mwl/(Zx.*R*T))*10^-6;
Vapour_density=(P.*Mwv/(Zy.*R*T))*10^-6;
Mwa=sum(comp.*Mw);
surf_tension=(sum(Parachor.*(x.*(liquid_density/Mwl)-y.*(Vapour_density/Mwv))))^4;

disp(' ');
disp(' ');
n2='N2  ';co2='CO2 ';c1='C1  ';c2='C2  ';c3='C3  ';ic4='IC4 ';nc4='NC4 ';ic5='IC5 ';nc5='NC5 ';c6='C6  ';c7='C7  ';c8='C8  ';c9='C9  ';c10='C10 ';c11='C11 ';c12='C12+';
compo=[n2;co2;c1;c2;c3;ic4;nc4;ic5;nc5;c6;c7;c8;c9;c10;c11;c12];
empty=['   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';'   ';];
ccc=comp(1:15);
ccc(16)=comp(21);
xxx=x(1:15);
xxx(16)=x(21);
yyy=y(1:15);
yyy(16)=y(21);
KKK=K(1:15);
KKK(16)=K(21);
disp('name  mole fr.      x           y          K')
disp([compo,empty,num2str(ccc'),empty,num2str(xxx'),empty,num2str(yyy'),empty,num2str(KKK')])
