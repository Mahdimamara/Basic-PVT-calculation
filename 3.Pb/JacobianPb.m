function JacP = JacobianPb(P,T,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,Fv)
x=comp./(1+Fv*(K-1));
Pplus=exp(1.01*log(P));
Pmin=exp(0.99*log(P));
y=x.*K;
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab);
[Slsplu,atlplus,btlplus,Alplus,Blplus] = coefficientcal(x,Pplus,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zlplus] = solvecub(Alplus,Blplus,C,b,btlplus,Slsplu,atlplus,sigma1,sigma2,x,R,T);
[Philplus] = FugacityCal(Zlplus,Alplus,Blplus,sigma1,sigma2,b,btlplus,Slsplu,atlplus);
[Slmin,atlmin,btlmin,Almin,Blmin] = coefficientcal(x,Pmin,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zlmin] = solvecub(Almin,Blmin,C,b,btlmin,Slmin,atlmin,sigma1,sigma2,x,R,T);
[Philmin] = FugacityCal(Zlmin,Almin,Blmin,sigma1,sigma2,b,btlmin,Slmin,atlmin);
[Svplus,atvplus,btvplus,Avplus,Bvplus] = coefficientcal(y,Pplus,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zvplus] = solvecub(Avplus,Bvplus,C,b,btvplus,Svplus,atvplus,sigma1,sigma2,y,R,T);
[Phivplus] = FugacityCal(Zvplus,Avplus,Bvplus,sigma1,sigma2,b,btvplus,Svplus,atvplus);
[Svmin,atvmin,btvmin,Avmin,Bvmin] = coefficientcal(y,Pmin,T,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zvmin] = solvecub(Avmin,Bvmin,C,b,btvmin,Svmin,atvmin,sigma1,sigma2,y,R,T);
[Phivmin] = FugacityCal(Zvmin,Avmin,Bvmin,sigma1,sigma2,b,btvmin,Svmin,atvmin);
Kpplus=Philplus./Phivplus;
Kpmin=Philmin./Phivmin;
Dpplus=sum((Kpplus-1).*comp./(1+Fv*(Kpplus-1)));
Dpmin=sum((Kpmin-1).*comp./(1+Fv*(Kpmin-1)));
JacP=(Dpplus-Dpmin)/(0.02*log(P));
end

