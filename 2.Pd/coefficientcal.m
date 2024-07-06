function [S,at,bt,A,B] = coefficientcal(comp,P,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha)
S=zeros(1,Nc);
for i=1:Nc;
    s=sum(comp.*(1-K_ij(i,:)).*sqrt(a));
    S(i)=sqrt(a(i))*s;
end
at=sum(comp.*S);
bt=sum(comp.*b);
A=(at*P)/(R*T)^2;
B=(bt*P)/(R*T);
end


