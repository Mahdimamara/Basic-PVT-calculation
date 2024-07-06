% CCE Experiment
clc;clear;format bank
disp('CCE Test')
disp('Results from Pd calculation:')

fluidtype=1;
EOS=1;
TTd = 580; %K
[PPd,Vsat]=Pdcalc(TTd);
disp(['dew point pressure: ',num2str(PPd),'  and calculated Volume at this point: ',num2str(Vsat)]);
P=[5824 5624 5424 5224 5024 4824 4624 4424 4224 4024 3624 3424 3224 3024 2824 2624 2424 2224 2024 1824].*101325./14.7;   
XX=[5824 5624 5424 5224 5024 4824 4624 4424 4224 4024 3624 3424 3224 3024 2824 2624 2424 2224 2024 1824];
for i=1:10
    PP=P(i);   
    TT=580;    %K
    [x,y,Fv,K,Zx,Zy,Vv,Vl] = flash(PP,TT);
ZZy(i)=Zy;
VVv(i)=Zy*8.314*580/XX(i);

end
for i=11:20
    PP=P(i);
    TT=580;
    [x,y,Fv,K,Zx,Zy,Vv,Vl] = flash(PP,TT);
    
VV(i)=(8.314*580/XX(i))*(Fv*Zy+(1-Fv)*Zx);
end
X1=[5824 5624 5424 5224 5024 4824 4624 4424 4224 4024];
X2=[3624  3424 3224 3024 2824 2624 2424 2224 2024 1824];
X3=VV./Vsat;
X4=X3(11:20);
disp(' ');
disp('       Pressure    Relative volume    Z factor');
disp([X1',(VVv./Vsat)',ZZy']);
disp([PPd,1]);
disp([X2',X4']);


