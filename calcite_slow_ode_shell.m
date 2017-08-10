function SU=calcite_slow_ode_shell(a,tj,msh)
% a=x(:,j), tj=t(j), msh=mesh 
%a=[0.006,0,0,0,1e-7,0,0,2e-3,4e-3,1e-7,0]';
%tj=0.01;
%msh=0.01;
tspan=[tj:msh:tj+2*msh];
[~,xU]=ode45(@calcite_slow_ode,tspan,a); % Can't do [tj-msh,tj]. Why?
    % We played with Solve_Calc_prec.m using tspan=[0,1], and
    % it would return a lengthy set of values instead of just
    % x at two time steps.  However, [0:0.5:1] would return 3 values.
    % xU looks like [x1(0)&x2(0)//x1(1)&x2(1)//...//x1(n)&x2(n)]
    % n=2 based on how we define tspan, and we want to retrieve the
    % final row, t(j+1)
SU=xU(2,1:11);  % Retrieves x(i) values at t(j+1) for 1<=i<=11

