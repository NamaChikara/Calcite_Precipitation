function blahblah=calcite_newton_shell(a,b,c,d,e)
% a=x(8),b=x(9),c=NT, d=CT, e=previous H+ concentration
K1=10^(-6.35);  % Fast reaction rate constants
K2=10^(-10.33); % ^
K3=10^(9.25);   % ^^
Kw=10^(-14);    % Water dissasociation constant
[g,eg,eeg]=calcite_newton(@(C) (c*K3*C)/(1+K3*C)+C-Kw/C-d/((C/K1)...    
    +1+(K2/C))-2*d/((C^2)/(K1*K2)+(C/K2)+1)+2*a-b,...
    @(C) (c*K3)/(1+K3*C)^2+1+Kw/C^2 ...
    +d*(1/K1-K2/C^2)/(C/K1+1+K2/C)^2 ...
    +2*d*((2*C)/(K1*K2)+1/K2)/((C)^2/(K1*K2)+C/K2+1)^2,e,1e-9,40);
                % Equation for f([H+]), df([H+]), initial
                %  approximation, tolerance, and max iterations.
blahblah=g;     % Return the final iteration of the method
                %  (the most accurate value for [H+]).
