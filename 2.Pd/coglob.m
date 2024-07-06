function [b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab)
b=(omegab.*R)*(Tc./Pc);
ac=(omegaa.*((R*Tc).^2))./Pc;
Tr=T./Tc;
alpha=(1+m.*(1-sqrt(Tr))).^2;
a=ac.*alpha;
end

