function [Tc, Pc, W, Mw, comp, K_ij, Nc, Mw_t, R] = FluidPT(fluidtype)
R=8.314;
switch fluidtype
    case 1
P = 4000/14.7*101325;
%240 F
T = (240+459.67)/1.8;   %k

Tc =  [126.20, 304.20, 190.60, 305.40, 369.8, 408.10, 425.20, 460.40, 469.60, 507.50, 543.20, 570.50, 598.50, 622.10, 643.60, 663.90, 682.40 , 700.70, 718.60, 734.50, 804.229];  %K

Pc = [33.500, 72.800, 45.400, 48.200, 41.900, 36.000, 37.500, 33.400, 33.300, 32.460, 30.970, 29.120, 26.940, 25.010, 23.170, 21.630, 20.430, 19.330, 18.250, 17.150, 14.198]*101325;  %pa

W = [0.040000, 0.225000, 0.008000, 0.098000, 0.152000, 0.176000, 0.193000, 0.227000, 0.251000, 0.275040, 0.308301, 0.351327, 0.390781, 0.443774, 0.477482, 0.522263, 0.559558, 0.604823, 0.651235, 0.683728, 0.761414];

Mw = [28.0130, 44.0100, 16.0430, 30.0700, 44.0970, 58.1240, 58.1240, 72.1510, 72.1510, 86.0000, 96.0000, 107.0000, 121.0000, 134.0000, 147.0000, 161.0000, 175.0000, 190.0000, 206.0000, 222.0000, 265.0000];

comp = [0.0312, 0.0323, 0.6976, 0.0903, 0.0402, 0.0081, 0.0144, 0.006, 0.0055, 0.0096, 0.011, 0.0127, 0.0086, 0.0061, 0.0029, 0.00, 0.00, 0.00, 0.000, 0.000, 0.0234];

K_ij = [0	   -0.02	0.031	    0.042	    0.091	    0.095	    0.095	    0.095	    0.095	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12
       -0.02	0	    0.103	    0.13	    0.135	    0.13	    0.13	    0.125	    0.125	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15
        0.031	0.103	0	        0.002689	0.008537	0.015715	0.014749	0.020879	0.020641	0.025345	0.029573	0.033997	0.039305	0.044372	0.049494	0.054287	0.05852	    0.062667	0.067047	0.071547	0.085936
        0.042	0.13	0.002689	0	        0.001662	0.005486	0.004914	0.008734	0.008578	0.011748	0.014731	0.01796	    0.021954	0.025867	0.029906	0.033751	0.037195	0.040607	0.044248	0.048026	0.060324
        0.091	0.135	0.008537	0.001662	0	        0.001117	0.000866	0.002801	0.002712	0.00462	    0.006573	0.008807	0.011699	0.014638	0.017758	0.020794	0.023561	0.026339	0.02934	    0.03249	    0.042952
        0.095	0.13	0.015715	0.005486	0.001117	0	        0.000016	0.000382	0.00035	    0.0012	    0.002286	0.003679	0.005634	0.007741	0.010073	0.012413	0.014593	0.016821	0.019265	0.021866	0.03071
        0.095	0.13	0.014749	0.004914	0.000866	0.000016	0	        0.000554	0.000515	0.001492	0.002682	0.004176	0.006244	0.008452	0.010879	0.013303	0.015555	0.01785	    0.020362	0.023029	0.032071
        0.095	0.125	0.020879	0.008734	0.002801	0.000382	0.000554	0	        0.000001	0.000228	0.0008	    0.001694	0.003092	0.004703	0.006561	0.008482	0.010309	0.012205	0.014311	0.016578	0.024435
        0.095	0.125	0.020641	0.008578	0.002712	0.00035	    0.000515	0.000001	0	        0.000255	0.000849	0.001765	0.003187	0.004819	0.006698	0.008637	0.010479	0.01239	    0.014511	0.016793	0.024693
        0.12	0.15	0.025345	0.011748	0.00462	    0.0012	    0.001492	0.000228	0.000255	0	        0.000174	0.00068	    0.001643	0.002866	0.004354	0.005947	0.007496	0.00913	    0.01097	    0.012974	0.02005
        0.12	0.15	0.029573	0.014731	0.006573	0.002286	0.002682	0.0008	    0.000849	0.000174	0	        0.000166	0.000749	0.00163	    0.002793	0.004096	0.0054	    0.006803	0.008409	0.01018	    0.016561
        0.12	0.15	0.033997	0.01796	    0.008807	0.003679	0.004176	0.001694	0.001765	0.00068	    0.000166	0	        0.00021	    0.000757	0.0016	    0.002618	0.00368	    0.004856	0.006229	0.007769	0.013462
        0.12	0.15	0.039305	0.021954	0.011699	0.005634	0.006244	0.003092	0.003187	0.001643	0.000749	0.00021	    0	        0.00017	    0.000652	0.001348	0.002137	0.003054	0.004163	0.005441	0.010351
        0.12	0.15	0.044372	0.025867	0.014638	0.007741	0.008452	0.004703	0.004819	0.002866	0.00163	    0.000757	0.00017	    0	        0.000156	0.000562	0.001104	0.001787	0.002657	0.003697	0.007894
        0.12	0.15	0.049494	0.029906	0.017758	0.010073	0.010879	0.006561	0.006698	0.004354	0.002793	0.0016	    0.000652	0.000156	0	        0.000125	0.00043	    0.000887	0.001526	0.002337	0.005844
        0.12	0.15	0.054287	0.033751	0.020794	0.012413	0.013303	0.008482	0.008637	0.005947	0.004096	0.002618	0.001348	0.000562	0.000125	0	        0.000091	0.000346	0.000778	0.001382	0.004266
        0.12	0.15	0.05852	    0.037195	0.023561	0.014593	0.015555	0.010309	0.010479	0.007496	0.0054	    0.00368	    0.002137	0.001104	0.00043	    0.000091	0	        0.000082	0.000337	0.000764	0.003116
        0.12	0.15	0.062667	0.040607	0.026339	0.016821	0.01785	    0.012205	0.01239	    0.00913	    0.006803	0.004856	0.003054	0.001787	0.000887	0.000346	0.000082	0	        0.000086	0.000346	0.002189
        0.12	0.15	0.067047	0.044248	0.02934	    0.019265	0.020362	0.014311	0.014511	0.01097	    0.008409	0.006229	0.004163	0.002657	0.001526	0.000778	0.000337	0.000086	0	        0.000086	0.001407
        0.12	0.15	0.071547	0.048026	0.03249	    0.021866	0.023029	0.016578	0.016793	0.012974	0.01018	    0.007769	0.005441	0.003697	0.002337	0.001382	0.000764	0.000346	0.000086	0	        0.000797
        0.12	0.15	0.085936	0.060324	0.042952	0.03071	    0.032071	0.024435	0.024693	0.02005	    0.016561	0.013462	0.010351	0.007894	0.005844	0.004266	0.003116	0.002189	0.001407	0.000797	0];    
    
    case 2
P = 2500/14.7*101325;

T = (560+459.67)/1.8;

Pc = [33.500, 72.800, 45.400, 48.200 , 41.900, 36.000, 37.500, 33.400, 33.300, 32.460, 30.970, 29.120, 26.940, 25.010, 23.170, 21.630, 20.430, 19.330, 18.250, 17.150, 16.350, 15.650, 15.060, 10.394]*101325;

Tc = [126.200, 304.200, 190.600, 305.400, 369.8000, 408.100, 425.200, 460.400, 469.600, 507.500, 543.200, 570.500 ,598.500, 622.100, 643.600, 663.900, 682.400, 700.700, 718.600, 734.500, 749.200, 760.500, 771, 905.03400];

W = [0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.27504, 0.308301, 0.351327, 0.390781, 0.443774, 0.477482, 0.522263, 0.559558, 0.604823, 0.651235, 0.683728, 0.72857, 0.757409, 0.790075, 1.054091];

Mw = [28.0130, 44.0100, 16.0430, 30.0700, 44.0970, 58.1240, 58.1240, 72.1510, 72.1510, 86.0000, 96.0000, 107.0000, 121.0000, 134.0000, 147.0000, 161.0000, 175.0000, 190.0000, 206.0000, 222.0000, 237.0000, 251.0000, 263.0000, 392.0000];

comp = [0.0067, 0.0211, 0.3493, 0.07, 0.0782, 0.038, 0.0168, 0.0202, 0.0178, 0.0304, 0.0439, 0.0471, 0.0321, 0.0179, 0.0172, 0.0174, 0.0174, 0.0135, 0.0134, 0.0106, 0.0102, 0.01, 0.009, 0.0918];

K_ij = [0	   -0.02	0.031	    0.042	    0.091	    0.095	    0.095	    0.095	    0.095	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12	    0.12
       -0.02	0	    0.103	    0.13	    0.135	    0.13	    0.13	    0.125	    0.125	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15	    0.15
        0.031	0.103	0	        0.002689	0.008537	0.015715	0.014749	0.020879	0.020641	0.025345	0.029573	0.033997	0.039305	0.044372	0.049494	0.054287	0.05852	    0.062667	0.067047	0.071547	0.075197	0.078505	0.081501	0.108804
        0.042	0.13	0.002689	0	        0.001662	0.005486	0.004914	0.008734	0.008578	0.011748	0.014731	0.01796	    0.021954	0.025867	0.029906	0.033751	0.037195	0.040607	0.044248	0.048026	0.051116	0.053934	0.056501	0.080426
        0.091	0.135	0.008537	0.001662	0	        0.001117	0.000866	0.002801	0.002712	0.00462	    0.006573	0.008807	0.011699	0.014638	0.017758	0.020794	0.023561	0.026339	0.02934	    0.03249	    0.03509	    0.03748	    0.039669	0.060574
        0.095	0.13	0.015715	0.005486	0.001117	0	        0.000016	0.000382	0.00035	    0.0012	    0.002286	0.003679	0.005634	0.007741	0.010073	0.012413	0.014593	0.016821	0.019265	0.021866	0.024036	0.026048	0.027905	0.046118
        0.095	0.13	0.014749	0.004914	0.000866	0.000016	0	        0.000554	0.000515	0.001492	0.002682	0.004176	0.006244	0.008452	0.010879	0.013303	0.015555	0.01785	    0.020362	0.023029	0.025252	0.02731	    0.029207	0.047752
        0.095	0.125	0.020879	0.008734	0.002801	0.000382	0.000554	0	        0.000001	0.000228	0.0008	    0.001694	0.003092	0.004703	0.006561	0.008482	0.010309	0.012205	0.014311	0.016578	0.018488	0.020269	0.021922	0.038479
        0.095	0.125	0.020641	0.008578	0.002712	0.00035	    0.000515	0.000001	0	        0.000255	0.000849	0.001765	0.003187	0.004819	0.006698	0.008637	0.010479	0.01239	    0.014511	0.016793	0.018713	0.020505	0.022167	0.038797
        0.12	0.15	0.025345	0.011748	0.00462	    0.0012	    0.001492	0.000228	0.000255	0	        0.000174	0.00068	    0.001643	0.002866	0.004354	0.005947	0.007496	0.00913	    0.01097	    0.012974	0.014677	0.016276	0.017769	0.033006
        0.12	0.15	0.029573	0.014731	0.006573	0.002286	0.002682	0.0008	    0.000849	0.000174	0	        0.000166	0.000749	0.00163	    0.002793	0.004096	0.0054	    0.006803	0.008409	0.01018	    0.011699	0.013137	0.014486	0.028547
        0.12	0.15	0.033997	0.01796	    0.008807	0.003679	0.004176	0.001694	0.001765	0.00068	    0.000166	0	        0.00021	    0.000757	0.0016	    0.002618	0.00368	    0.004856	0.006229	0.007769	0.009107	0.010384	0.011592	0.024482
        0.12	0.15	0.039305	0.021954	0.011699	0.005634	0.006244	0.003092	0.003187	0.001643	0.000749	0.00021	    0	        0.00017	    0.000652	0.001348	0.002137	0.003054	0.004163	0.005441	0.006573	0.007668	0.008714	0.020265
        0.12	0.15	0.044372	0.025867	0.014638	0.007741	0.008452	0.004703	0.004819	0.002866	0.00163	    0.000757	0.00017    	0	        0.000156	0.000562	0.001104	0.001787	0.002657	0.003697	0.004642	0.005571	0.00647	    0.016799
        0.12	0.15	0.049494	0.029906	0.017758	0.010073	0.010879	0.006561	0.006698	0.004354	0.002793	0.0016	    0.000652	0.000156	0	        0.000125	0.00043	    0.000887	0.001526	0.002337	0.0031	    0.003869	0.004625	0.013765
        0.12	0.15	0.054287	0.033751	0.020794	0.012413	0.013303	0.008482	0.008637	0.005947	0.004096	0.002618	0.001348	0.000562	0.000125	0	        0.000091	0.000346	0.000778	0.001382	0.001982	0.002606	0.003234	0.011298
        0.12	0.15	0.05852	    0.037195	0.023561	0.014593	0.015555	0.010309	0.010479	0.007496	0.0054	    0.00368	    0.002137	0.001104	0.00043	    0.000091	0	        0.000082	0.000337	0.000764	0.001225	0.001725	0.002243	0.009384
        0.12	0.15	0.062667	0.040607	0.026339	0.016821	0.01785	    0.012205	0.01239	    0.00913	    0.006803	0.004856	0.003054	0.001787	0.000887	0.000346	0.000082	0	        0.000086	0.000346	0.000673	0.001056	0.001468	0.007726
        0.12	0.15	0.067047	0.044248	0.02934	    0.019265	0.020362	0.014311	0.014511	0.01097	    0.008409	0.006229	0.004163	0.002657	0.001526	0.000778	0.000337	0.000086	0	        0.000086	0.000277	0.000538	0.000843	0.00619
        0.12	0.15	0.071547	0.048026	0.03249	    0.021866	0.023029	0.016578	0.016793	0.012974	0.01018	    0.007769	0.005441	0.003697	0.002337	0.001382	0.000764	0.000346	0.000086	0	        0.000054	0.000193	0.00039	    0.004822
        0.12	0.15	0.075197	0.051116	0.03509	    0.024036	0.025252	0.018488	0.018713	0.014677	0.011699	0.009107	0.006573	0.004642	0.0031	    0.001982	0.001225	0.000673	0.000277	0.000054	0	        0.000043	0.000153	0.003859
        0.12	0.15	0.078505	0.053934	0.03748	    0.026048	0.02731	    0.020269	0.020505	0.016276	0.013137	0.010384	0.007668	0.005571	0.003869	0.002606	0.001725	0.001056	0.000538	0.000193	0.000043	0	        0.000034	0.003091
        0.12	0.15	0.081501	0.056501	0.039669	0.027905	0.029207	0.021922	0.022167	0.017769	0.014486	0.011592	0.008714	0.00647	    0.004625	0.003234	0.002243	0.001468	0.000843	0.00039	    0.000153	0.000034	0	        0.002478
        0.12	0.15	0.108804	0.080426	0.060574	0.046118	0.047752	0.038479	0.038797	0.033006	0.028547	0.024482	0.020265	0.016799	0.013765	0.011298	0.009384	0.007726	0.00619	    0.004822	0.003859	0.003091	0.002478	0];    
   
    case 3

P = 1800/14.7*101325;

T = (600+459.67)/1.8;

Pc = [33.500, 72.800, 45.400, 48.200, 41.900, 36.000, 37.500, 33.400, 33.300, 32.460, 18.016]*101325;

Tc = [126.200, 304.20, 190.60, 305.40, 369.800, 408.100, 425.200, 460.400, 469.600, 507.500, 708.800];

W = [0.04, 0.225, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.27504, 0.565824];

Mw=[28.0130 44.0100 16.0430 30.0700 44.0970 58.1240 58.1240 72.1510 72.1510 86.0000 190.0000];

comp=[0.0051, 0.0119, 0.4521, 0.0709, 0.0461, 0.0169, 0.0281, 0.0155, 0.0201, 0.0442, 0.2891];

K_ij=[ 0	   -0.02	0.031	    0.042	    0.091	    0.095	    0.095	    0.095	    0.095	    0.12	  0.12
      -0.02	    0	    0.103	    0.13	    0.135	    0.13	    0.13	    0.125	    0.125	    0.15	  0.15
       0.031	0.103	0	        0.002689	0.008537	0.015715	0.014749	0.020879	0.020641	0.025345  0.066514
       0.042	0.13	0.002689	0	        0.001662	0.005486	0.004914	0.008734	0.008578	0.011748  0.043803
       0.091	0.135	0.008537	0.001662	0	        0.001117	0.000866	0.002801	0.002712	0.00462	  0.028971
       0.095	0.13	0.015715	0.005486	0.001117	0	        0.000016	0.000382	0.00035	    0.0012	  0.018963
       0.095	0.13	0.014749	0.004914	0.000866	0.000016	0	        0.000554	0.000515	0.001492  0.020051
       0.095	0.125	0.020879	0.008734	0.002801	0.000382	0.000554	0	        0.000001	0.000228  0.014049
       0.095	0.125	0.020641	0.008578	0.002712	0.00035	    0.000515	0.000001	0	        0.000255  0.014248
       0.12	    0.15	0.025345	0.011748	0.00462	    0.0012	    0.001492	0.000228	0.000255	0	      0.01074
       0.12	    0.15	0.066514	0.043803	0.028971	0.018963	0.020051	0.014049	0.014248	0.01074	  0];
end
Nc = size(comp,2); 
mw0 = Mw.*comp;
Mw_t = sum(mw0);


end

