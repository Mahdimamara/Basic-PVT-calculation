function JacT = JacobianTd(P,T,comp,R,Tc,Pc,omegaa,omegab,sigma1,sigma2,K_ij,m,C,Nc,K,Fv)
y=comp.*K./(1+Fv*(K-1));
Tplus=exp(1.01*log(T));
Tmin=exp(0.99*log(T));
x=y./K;
[b,a,ac,alpha] = coglob(R,Tc,Pc,m,T,omegaa,omegab);
[Slsplu,atlplus,btlplus,Alplus,Blplus] = coefficientcal(x,P,Tplus,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zlplus] = solvecub(Alplus,Blplus,C,b,btlplus,Slsplu,atlplus,sigma1,sigma2,x,R,Tplus);
[Philplus] = FugacityCal(Zlplus,Alplus,Blplus,sigma1,sigma2,b,btlplus,Slsplu,atlplus);
[Slmin,atlmin,btlmin,Almin,Blmin] = coefficientcal(x,P,Tmin,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zlmin] = solvecub(Almin,Blmin,C,b,btlmin,Slmin,atlmin,sigma1,sigma2,x,R,Tmin);
[Philmin] = FugacityCal(Zlmin,Almin,Blmin,sigma1,sigma2,b,btlmin,Slmin,atlmin);
[Svplus,atvplus,btvplus,Avplus,Bvplus] = coefficientcal(y,P,Tplus,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zvplus] = solvecub(Avplus,Bvplus,C,b,btvplus,Svplus,atvplus,sigma1,sigma2,y,R,Tplus);
[Phivplus] = FugacityCal(Zvplus,Avplus,Bvplus,sigma1,sigma2,b,btvplus,Svplus,atvplus);
[Svmin,atvmin,btvmin,Avmin,Bvmin] = coefficientcal(y,P,Tmin,Pc,Tc,m,R,omegaa,omegab,Nc,K_ij,b,ac,a,alpha);
[Zvmin] = solvecub(Avmin,Bvmin,C,b,btvmin,Svmin,atvmin,sigma1,sigma2,y,R,Tmin);
[Phivmin] = FugacityCal(Zvmin,Avmin,Bvmin,sigma1,sigma2,b,btvmin,Svmin,atvmin);
KTplus=Philplus./Phivplus;
KTmin=Philmin./Phivmin;
DTplus=sum((KTplus-1).*comp./(1+Fv*(KTplus-1)));
DTmin=sum((KTmin-1).*comp./(1+Fv*(KTmin-1)));
JacT=(DTplus-DTmin)/(0.02*log(T));

end

