clc;clear;
% EOS=input('choose the EOS function:\n SRK(1),SRK&G(2),PR76(3),PR78(4)\n');
% fluidtype=input('choose the fluid type:\n Condensate(1), Black oil(2),Volatile oil(3)\n');
Fvd=1;
%%
disp('dew points diagram');
disp('Equation of state: SRK');
disp('Fluid type: Gas condensate');
EOS=1;
fluidtype=1;
[Tc, Pc, W, Mw, comp, K_ij, Nc, Mw_t, R] = FluidPT(fluidtype);
[m,C,sigma1,sigma2,omegaa,omegab] = EQchar(W,EOS,Mw,Nc);
%%
countd=1;
countdT=1;
Sloped=1;
Errord=2;
T_d=zeros(1,1000);
P_d=zeros(1,1000);
K_d=zeros(1000,Nc);
if fluidtype==2
    Td=475;
else
    Td=450;
end
%%
for i=1:10000
    if countd<3
        Td=Td+5;
        Pd=1/sum(comp./(Pc.*exp(5.37*(1+W).*(1-Tc/Td))));
        Kd=(Pc./Pd).*exp(5.37*(1+W).*(1-Tc/Td));
        
        [P_d(countd),K_d(countd,:)] =Pdcalculator1(R,Tc,Pc,m,Td,Pd,Kd,omegaa,omegab,sigma1,sigma2,K_ij,Nc,C,comp,Fvd);
        T_d(countd)=Td;
        countd=countd+1;
        
    else
        if Fvd==1
            a=10^4;
        else
            a=10^5;
        end
        if abs(Sloped)<(a)
            if Sloped>0
                Td=Td+3;
            else
                Td=Td;
            end
            
            
            Sloped = ((P_d(countd-1) - P_d(countd-2))/(T_d(countd-1) - T_d(countd-2)));
            
            Hd=log(K_d(countd-1,:));
            [JacKd,JacTd] = JacobianfuncKTd(P_d(countd-1),T_d(countd-1),comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K_d(countd-1,:),Fvd);
            Kd=exp(Hd-(log(Td)-log(T_d(countd-1))).*(JacTd/JacKd));
            E=log(P_d(countd-1));
            JacobPd = JacobianPd(P_d(countd-1),Td,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K_d(countd-1,:),Fvd);
            Dpd = sum((Kd-1).*comp./(1+Fvd.*(Kd-1)));
            Pd=exp(E-Dpd/JacobPd);
            
            [P_d(countd), K_d(countd,:)] =Pdcalculator(R,Tc,Pc,m,Td,Pd,Kd,omegaa,omegab,sigma1,sigma2,K_ij,Nc,C,comp,Fvd);
            T_d(countd)=Td;
            countd=countd+1;
        else
            
            Sloped = ((P_d(countd-1) - P_d(countd-2))/(T_d(countd-1) - T_d(countd-2)));
            
            Kd=K_d(countd-1,:)+((K_d(countd-1,:))-(K_d(countd-2,:)))./((P_d(countd-1))-(P_d(countd-2)))*((Pd)-(P_d(countd-1)));
            Td = T_d(countd-1)+(1/Sloped)*(Pd - P_d(countd-1));
            a=600;
            if countd>a
                if abs(Td-T_d(countd-1))<0.0001
                    break
                end
            end
            
            Pd=P_d(countd-1)+(P_d(countd-1)-P_d(countd-2));
            
            Kd=K_d(countd-1,:)+((K_d(countd-1,:))-(K_d(countd-2,:)))./((P_d(countd-1))-(P_d(countd-2)))*((Pd)-(P_d(countd-1)));
            Td = T_d(countd-1)+(1/Sloped)*(Pd - P_d(countd-1));
            
            [T_d(countd), K_d(countd,:)] = Tdcalculator(R,Tc,Pc,m,Td,Pd,Kd,omegaa,omegab,sigma1,sigma2,K_ij,Nc,C,comp,Fvd);
            P_d(countd)=Pd;
            countd=countd+1;
        end
    end
    if EOS==3 || EOS==4
        a=0.1;
    else
        a=2;
    end
    if fluidtype==3
        a=0.5*a;
    end
    if Fvd~=1
        a=a*2;
    end
    Errord=((sum(abs(K_d(countd-1,:)-1)))<(a));
    
end
T_d= T_d(T_d>0);
P_d=P_d(P_d>0)./101325.*14.7 ;   %psi
hold on
%%
if Fvd==1
    plot(T_d*1.8-460,P_d,'b*')
elseif Fvd==0.9
    plot(T_d*1.8-460,P_d,'K.')
elseif Fvd==0.8
    plot(T_d*1.8-460,P_d,'g.')
end
xlabel('Temrature K');
ylabel('Pressure psia');
grid on