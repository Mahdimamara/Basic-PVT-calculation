function [Pd, K] =Pdcalculator(R,Tc,Pc,m,Td,Pd,K,omegaa,omegab,sigma1,sigma2,K_ij,Nc,C,comp,Fvd)
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,Td,omegaa,omegab);
xd=zeros(1,Nc);
yd=zeros(1,Nc);
Dd=1;
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
    Dd=sum(comp.*(K-1)./(1+Fvd*(K-1)));
    if Dd<10^-5
        break
    else
        JacP = JacobianPd(Pd,Td,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,Fvd);
        E=log(Pd);
        Pd=exp(E-Dd/JacP);
    end
end
end

