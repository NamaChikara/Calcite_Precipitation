function SRU=calcite_slow_ode(n,b)
% n=t, b=x
SRU=zeros(11,1);  % % create zero column to store outputs
ku=2.3148e-6; % urea disollution constant=0.2(day^(-1)) in seconds^(-1)
SRU(1)=-ku*b(1);  % x(1)=CO(NH2)2 ODE
SRU(2)=2*ku*b(1); % x(2)=NH3 ODE
SRU(3)=ku*b(1);   % x(3)=H2CO3 ODE
SRU(4)=0;   % These don't participate in the slow reactions
SRU(5)=0;
SRU(6)=0;
SRU(9)=0;
SRU(10)=0;
SRU(7)=deq2(b(7),b(8));  % x(7)=CO3(2-) ODE
SRU(8)=deq2(b(7),b(8));  % x(8)=Ca(2+) ODE
SRU(11)=-deq2(b(7),b(8));  % calcite formation is equal to loss  
                           % of precipitated ions 7 and 8