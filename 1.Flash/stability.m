function [K,u] = stability(P,T,comp,Tc,Pc,W,R,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,fluidtype,EOS)
K=(Pc./P).*exp(5.37*(1+W).*(1-(Tc./T)));
Fv=0.5;
x=comp;
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab);
[Sx,atx,btx,Ax,Bx] = coefficientcal(x,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zx] = solvecub(Ax,Bx,C,b,btx,Sx,atx,sigma1,sigma2,x,R,T);
[Phix] = FugacityCal(Zx,Ax,Bx,sigma1,sigma2,b,btx,Sx,atx);
Tol=1;
D=1;
while Tol>10^-4
    
    switch fluidtype
        case 1
            u = x./K;
        case {2,3}
            u=K.*x;
    end
    U = u/sum(u);
    [SU,atU,btU,AU,BU] = coefficientcal(U,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [ZU] = solvecub(AU,BU,C,b,btU,SU,atU,sigma1,sigma2,U,R,T);
    [PhiU] = FugacityCal(ZU,AU,BU,sigma1,sigma2,b,btU,SU,atU);
    
    O=K;
    if fluidtype==1
        D = log(K)-log(PhiU)+log(Phix);
    else
        D = log(K)+log(PhiU)-log(Phix) ;
    end
    H=log(K);
    Jac = Jacobianfuncmatric(P,T,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,fluidtype,Fv);
    delta_Phsy=Jac\D';
    H=H-0.7*(delta_Phsy)';
    K=exp(H);
    Tol=sum(abs(K-O));

    
end
end

