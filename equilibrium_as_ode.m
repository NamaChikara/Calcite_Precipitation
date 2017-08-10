function SRU=equilibrium_as_ode(n,b)
% n=t, b=x
SRU=zeros(11,1);  % % create zero column to store outputs
s1=2.3148e-6; % urea disollution constant=0.2(day^(-1)) in seconds^(-1)
k1=10;        % fast reactions (see pdf image in Calcite_Precipitation folder
kn1=10^7.35;
k2=10;
kn2=10^11.33;
k3=10^11;
kn3=10^1.75;
SRU(1)=-s1*b(1);  % b(1)=CO(NH2)2 ODE
SRU(2)=2*s1*b(1)-k1*2*b(2)+kn1*2*b(4)*b(5); % b(2)=NH3 ODE
SRU(3)=s1*b(1)-k2*b(3)*2*b(5)+kn2*b(6)*b(5);   % b(3)=H2CO3 ODE
SRU(4)=k1*2*b(2)-kn1*2*b(4)*b(5);   % b(4)=NH4^+
SRU(5)=k1*2*b(2)-kn1*2*b(4)*b(5)-k2*2*b(3)*b(5)+(kn2-k3)*b(6)*b(5);   % b(5)=OH^-
SRU(6)=k2*b(3)*2*b(5)-k3*b(6)*b(5);   % b(6)=HCO_3^-
SRU(7)=k3*b(6)*b(5)-kn3*b(7)+calcite_precipitation_rate(b(7),b(8));  % b(7)=CO3(2-) ODE
SRU(8)=calcite_precipitation_rate(b(7),b(8));  % b(8)=Ca(2+) ODE
SRU(11)=-calcite_precipitation_rate(b(7),b(8));  % calcite formation is equal to loss  
                           % of precipitated ions 7 and 8