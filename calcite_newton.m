function [root,error,iter]=calcite_newton(f,g,x0,tol,max)
% f - anonymous function, g - anonymous derivative, x0 - initial guess, 
%   tol - tolerance, max - max iterations

x(1)=x0;                             % Initialize vector of x values.
x(2)=x(1)-f(x(1))/g(x(1));           % Manually compute x(2).
for i=3:max                          % Only go to max iterations.
    if abs(x(i-1)-x(i-2))>=tol       % If tolerance condition isn't met...
    x(i)=x(i-1)-f(x(i-1))/g(x(i-1)); %   compute next root estimate
    else                             % Otherwise...
        break                        %   leave nested loop.
    end
end

if length(x)==max           % Display error if method didn't succeed.
    error('Max iteration reached before error<tol.',class(n))
end

a=length(x);                % Determine # of root guesses.
e(1)=0;                     % Difference between successive errors.
for i=2:a                   % ^
    e(i)=abs(x(i-1))-abs(x(i));
end
%negative_E_diff=min(e);     % Check to see if error did not improve
%if negative_E_diff<0        %   at any iteration.  Error if so.
%    error('Error. Method is not improving root guess.',class(n))
%end
error=abs(x(a)-x(a-1));     % Compute difference between final 2 roots.
root=x(a);                  % Record final root value.
iter=a;                     % Record # of iterations.
    