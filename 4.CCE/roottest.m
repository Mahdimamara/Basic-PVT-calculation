function Z = roottest(Z1,Z2,Z3,x,T,P)
global R sigma2 sigma1 b
[aa,bb,S,A,B,a1,a2,a3] = coefficient(x,T,P)
Zroot=[Z1,Z2,Z3];
if size(Zroot,2)==3
H=sort(Zroot);
if H(1)>0
    Zl=H(1);
    Zh=H(3);
    delta_G=R*T*sum(x.*((b/bb)*(Zh-Zl)-log((Zh-B)/(Zl-B))-(A/(B*(sigma2-sigma1)))*(2*S/aa-b/bb)*log((Zh+sigma2*B)*(Zl+sigma1*B)/((Zh+sigma1*B)*(Zl+sigma2*B)))));
if delta_G>0
    Z=Zl
else
   Z=Zh
end 
elseif H(1)<0 && H(2)>0
    Zl=H(2);
    Zh=H(3);
    delta_G=R*T*sum(x.*((b/bb)*(Zh-Zl)-log((Zh-B)/(Zl-B))-(A/(B*(sigma2-sigma1)))*(2*S/aa-b/bb)*log((Zh+sigma2*B)*(Zl+sigma1*B)/((Zh+sigma1*B)*(Zl+sigma2*B)))));
if delta_G>0
    Z=Zl
else
   Z=Zh
end 
else if H(2)<0&& H(3)>0
        Z=H(3);   
    end
end
else
    Z=Z1;
       
end
end

