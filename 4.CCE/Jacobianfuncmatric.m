function [Jac] = Jacobianfuncmatric(P,T,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,fluidtype,Fv)
Jac1=zeros(Nc);
for i=1:Nc
    Kplus=K;
    Kmin=K;
    Kplus(i)=exp(log(K(i))+0.01*log(K(i)));
    Kmin(i)=exp(log(K(i))-0.01*log(K(i)));
    xplus=comp./(1+Fv*(Kplus-1));
    xmin=comp./(1+Fv*(Kmin-1));
    
    yplus=Kplus.*xplus;
    ymin=Kmin.*xmin;
    [b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab);
    [Sxplus,atxplus,btxplus,Axplus,Bxplus] = coefficientcal(xplus,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zxplus] = solvecub(Axplus,Bxplus,C,b,btxplus,Sxplus,atxplus,sigma1,sigma2,xplus,R,T);
    [Phixplus] = FugacityCal(Zxplus,Axplus,Bxplus,sigma1,sigma2,b,btxplus,Sxplus,atxplus);
    
    [Sxmin,atxmin,btxmin,Axmin,Bxmin] = coefficientcal(xmin,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zxmin] = solvecub(Axmin,Bxmin,C,b,btxmin,Sxmin,atxmin,sigma1,sigma2,xmin,R,T);
    [Phixmin] = FugacityCal(Zxmin,Axmin,Bxmin,sigma1,sigma2,b,btxmin,Sxmin,atxmin);
    
    [Syplus,atyplus,btyplus,Ayplus,Byplus] = coefficientcal(yplus,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zyplus] = solvecub(Ayplus,Byplus,C,b,btyplus,Syplus,atyplus,sigma1,sigma2,yplus,R,T);
    [Phiyplus] = FugacityCal(Zyplus,Ayplus,Byplus,sigma1,sigma2,b,btyplus,Syplus,atyplus);
    
    [Symin,atymin,btymin,Aymin,Bymin] = coefficientcal(ymin,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
    [Zymin] = solvecub(Aymin,Bymin,C,b,btymin,Symin,atymin,sigma1,sigma2,ymin,R,T);
    [Phiymin] = FugacityCal(Zymin,Aymin,Bymin,sigma1,sigma2,b,btymin,Symin,atymin);
    Dplus = log(Kplus)+log(Phiyplus)-log(Phixplus);
    Dmin = log(Kmin)+log(Phiymin)-log(Phixmin);
    Jac1(i,:)=(Dplus-Dmin)./(0.02*log(K(i)));
end
Jac=Jac1';
end

