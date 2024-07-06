function [Z] = solvecub(A,B,C,bi,b,S,a,sigma1,sigma2,x,R,T)
a1=-(1-B*C);
a2=(A-B*(1+C)-(B^2)*(1+2*C));
a3=-(A*B-C*(B^3+B^2));
Q=(3*(a2)-(a1)^2)/9;
J=(9*(a1)*(a2)-27*(a3)-2*(a1)^3)/54;
D=Q^3+J^2;
Z1=[];
Z2=[];
Z3=[];
if D>0
    Z1=nthroot((J+sqrt(D)),3)+nthroot((J-sqrt(D)),3)-(a1)/3;
elseif D<0
    tetha=acos(J/sqrt(-Q^3));
    Z1=2*sqrt(-Q)*cos(tetha/3)-(a1)/3;
    Z2=2*sqrt(-Q)*cos(tetha/3+2*pi/3)-(a1)/3;
    Z3=2*sqrt(-Q)*cos(tetha/3+4*pi/3)-(a1)/3;
else
    Z1=2*nthroot(J,3)-(a1)/3;
    Z2=-nthroot(J,3)-(a1)/3;
    Z3=Z2;
end
Zroot=[Z1 Z2 Z3];
if size(Zroot,2)==3
    Zsort=sort(Zroot);
    if Zsort(1)>0
        ZL=Zsort(1);
        ZH=Zsort(3);
        delta_G=R*T*sum(x.*(((bi/b)*(ZH-ZL)-log((ZH-B)/(ZL-B))-(A/(B*(sigma2-sigma1)))*(2*S/a-bi/b)*log((ZH+sigma2*B)*(ZL+sigma1*B)/((ZH+sigma1*B)*(ZL+sigma2*B))))));
        if delta_G<0
            Z=ZH;
        else
            Z=ZL;
        end
    elseif (Zsort(2)>0 && Zsort(1)<0)
        ZL=Zsort(2);
        ZH=Zsort(3);
        delta_G=R*T*sum(x.*(((bi/b)*(ZH-ZL)-log((ZH-B)/(ZL-B))-(A/(B*(sigma2-sigma1)))*(2*S/a-bi/b)*log((ZH+sigma2*B)*(ZL+sigma1*B)/((ZH+sigma1*B)*(ZL+sigma2*B))))));
        if delta_G<0
            Z=ZH;
        else
            Z=ZL;
        end
    elseif (Zsort(2)<0 && Zsort(3)>0)
        Z=Zsort(3);
    end
else
    Z=Z1;
end

end

