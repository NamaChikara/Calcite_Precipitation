function [x,t,l]=bigsolver2(n,t0,fin,x0)
% x0=[6E-3,0,0,0,1E-7,0,0,2E-3,4E-3,1E-7,0];
tic
% note that x0 should have 11 values, the 11th being for precipitate
% n is the number of grid points
% t0 and fin are the start and end times, respectively
% x=zeros(11,n);  % create x matrix to store n time steps
x(:,1)=x0;      % store initial x value
K1=10^(-6.35);  % Fast reaction rate constants
K2=10^(-10.33); % ^
K3=10^(9.25);  % ^^
Kw=10^(-14);    % Water dissasociation constant
t=zeros(1,n);   % Length of t     
t(1)=t0;        % Initial time value
mesh=fin/n;     % Step Size=time length/# grid points
for j=2:n
    t(j)=t(j-1)+mesh;               % current time step
    NT=x(2,j-1)+x(4,j-1);           % Total dissolved nitrogen (step 2).
    CT=x(3,j-1)+x(6,j-1)+x(7,j-1);  % Total dissolved carbon (step 2).
    x(10,j)=calcite_newton_shell(x(8,j-1),x(9,j-1),NT,CT,x(10,j-1));
    x(2,j)=NT/(K3*x(10,j)+1);                    % NH3 concentration.
    x(4,j)=(NT*K3*x(10,j))/(K3*x(10,j)+1);      % NH4+ concentration.
    c=x(10,j)^2+x(10,j)*K1+K1*K2;   % Variable to simplify carbon equations
    x(3,j)=CT/(c/x(10,j)^2);        % H2CO3 concentration.
    x(6,j)=(CT*K1)/(c/x(10,j));     % HCO3- concentration.
    x(7,j)=(CT*K1*K2)/c;            % CO3(2-) concentration.
    x(1,j)=x(1,j-1);                % > Not involved in fast reactions
    x(8,j)=x(8,j-1);                % > so we restore their previous
    x(9,j)=x(9,j-1);                % > time values as they do not
    x(11,j)=x(11,j-1);              % > change yet
    x(5,j)=x(4,j)+x(10,j)+2*x(8,j)-x(6,j)-2*x(7,j)-x(9,j); % OH- concentration.
    SlowUpd=calcite_slow_ode_shell(x(:,j),t(j),mesh); % Returns row vector with 11 
                                                % entries for x(i) at t(j+1)
    for i=1:11
        x(i,j)=SlowUpd(1,i);        % Update x column after slow ODEs
    end
end
toc
l=toc;
plot(t,x(1,:),'b',t,x(2,:),'r',t,x(3,:),'g',t,x(4,:),'b--',t,x(5,:),'r--',t,x(6,:),'g--',t,x(7,:),'b:',t,x(8,:),'r:',t,x(9,:),'g:',t,x(10,:),t,x(11,:))
legend('CO(NH_2)_2','NH_3','H_2CO_3','NH_4^+','OH^-','HCO_3^-','CO_3^{2-}','Ca^{2+}','Cl^-','H^+','CaCO_3')
end

    
    