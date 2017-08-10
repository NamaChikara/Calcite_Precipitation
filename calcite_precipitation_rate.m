function rate=calcite_precipitation_rate(a,b)  % x(7),x(8)
Sc=50;          % critical concentration
Kso=3.8e-9;
%kp=0.2;          % constant in units of mmol/L*day
kp=2.3148e-9;    % constant in units of mol/L*s
%kp=2e-4;         % constant in units of mmol/L*day
for S=(a*b)/Kso;   % check value of calcite saturation state
    if S>=Sc       % if exceeds critical concentration...
        rate=-kp*(S-1)^2; % F(7)=F(8) value
    else rate=0;
    end
end
end