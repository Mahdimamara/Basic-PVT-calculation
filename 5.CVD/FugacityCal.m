function [Phi] = FugacityCal(Z,A,B,sigma1,sigma2,b,bt,S,at)
lnPhi=(b./bt)*(Z-1)-log(Z-B)- ...
        (A/(B*(sigma2-sigma1)))*(2*(S./at)-b./bt)* ...
        log((Z+sigma2*B)/(Z+sigma1*B));
    Phi=exp(lnPhi);
end

