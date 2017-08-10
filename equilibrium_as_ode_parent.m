function [x,t,l]=equilibrium_as_ode_parent(n,t0,fin,x0)
% x0=[6E-3,0,0,0,1E-7,0,0,2E-3,4E-3,1E-7,0];
tic
% note that x0 should have 11 values, the 11th being for precipitate
% n is the number of grid points
% t0 and fin are the start and end times, respectively
% x=zeros(11,n);  % create x matrix to store n time steps
x(:,1)=x0;      % store initial x value
t=zeros(1,n);   % Length of t     
t(1)=t0;        % Initial time value
mesh=fin/n;     % Step Size=time length/# grid points
for j=2:n
    t(j)=t(j-1)+mesh;               % current time step
    SlowUpd=equilibrium_as_ode_shell(x(:,j-1),t(j-1),mesh); % Returns row vector with 11 
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

    
    