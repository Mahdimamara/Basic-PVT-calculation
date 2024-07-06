function [m,C,sigma1,sigma2,omegaa,omegab] = EQchar(W,EOS,Mw,Nc)
switch EOS
    case 3
    m = 0.37464+1.54226*W-0.26992*(W.^2); 
    C = 1;    
    sigma1 = 1-sqrt(2);
    sigma2 = 1+sqrt(2);
    omegaa = 0.457235;
    omegab = 0.077796; 
    case 4
        a=(Mw > 134);
        for i=1:Nc
     if a(i)==1
        m(i) = 0.37964+1.48503*W(i)-0.164423*(W(i).^2);
    else
        m(i) = 0.37964+1.54226*W(i)-0.26992*(W(i).^2);
     end
        end
    C = 1;    
    sigma1 = 1-sqrt(2);
    sigma2 = 1+sqrt(2);
    omegaa = 0.457235;
    omegab = 0.077796;    
    case 1
    m = 0.480+1.574*W-0.176*(W.^2); 
    C = 0;    
    sigma1 = 0;
    sigma2 = 1;
    omegaa = 0.42747;
    omegab = 0.08664;
    case 2
    m = 0.48508+1.55171*W-0.15613*(W.^2); 
    C = 0;    
    sigma1 = 0;
    sigma2 = 1;
    omegaa = 0.42747;
    omegab = 0.08664;  
end

end

