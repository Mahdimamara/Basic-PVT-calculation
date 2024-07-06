function [Td, K] = Tdcalculator(R,Tc,Pc,m,Td,Pd,K,omegaa,omegab,sigma1,sigma2,K_ij,Nc,C,comp,Fvd)
xd=zeros(1,Nc);
yd=zeros(1,Nc);
Told=1;
Dd=1;
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,Td,omegaa,omegab);
while Dd>10^-5
    yd=comp.*K./(1+Fvd*(K-1));
    xd=yd./K;
    [Sx,atx,btx,Ax,Bx] = coefficientcal(xd,Pd,Td,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Sy,aty,bty,Ay,By] = coefficientcal(yd,Pd,Td,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zx] = solvecub(Ax,Bx,C,b,btx,Sx,atx,sigma1,sigma2,xd,R,Td);
    [Zy] = solvecub(Ay,By,C,b,bty,Sy,aty,sigma1,sigma2,yd,R,Td);
    [Phix] = FugacityCal(Zx,Ax,Bx,sigma1,sigma2,b,btx,Sx,atx);
    [Phiy] = FugacityCal(Zy,Ay,By,sigma1,sigma2,b,bty,Sy,aty);
    K=Phix./Phiy;
    Dd=sum((K-1).*comp./(1+Fvd.*(K-1)));
    if Dd<10^-5
        break
    else
        JacT = JacobianTd(Pd,Td,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,Fvd);
        H=log(Td);
        Td=exp(H-Dd/JacT);
    end
end
end

