%Rohit Thirumala
%21036098

clc
clear all

[solution,iteration]=newtonFD(@myfun,@myjacFDF,[1;1])
[solution,iteration]=newton(@myfun,@myjac,[1;1])


function [jacobian] = myjacFDF(functions,vectorx, delta)
    jacobian_column1 = (functions([vectorx(1)*(1+delta),vectorx(2)])-functions([vectorx(1),vectorx(2)]))/(delta*vectorx(1));
    jacobian_column2 = (functions([vectorx(1),vectorx(2)*(1+delta)])-functions([vectorx(1),vectorx(2)]))/(delta*vectorx(2));
    jacobian = [jacobian_column1,jacobian_column2];
end

function F = myfun(x)
    %the input is a vector now 
    %x=[x1,x2]
	f1 = x(1)^3 -x(2)^2 - 55;
	f2 = x(1)*x(2) -12;
	F = [f1;f2];
end


function J = myjac(x)
    %x=[x1,x2]
    %for this specific problem
    %J=[3x1^2 -2x2;x2 x1]
    J=[3*x(1)^2 -2*x(2);x(2) x(1)];
	
end
function [sol,iter] = newton(f,J,x0)
    % Solve a system of nonlinear equations f(x) = 0, 
    %Jac is the Jacobian 
    % x0 is a vector of the initial guesses
    maxiter = 10; % maximum number of iterations before exiting the program
    tol = 1e-4; % tolerance for convergence
    eps = 1; % initial epsilon, any value OK as long as greater than the tolerance
    xold = x0; % initial guess
    i = 0; % initialize the number of iterations
 
    while eps > tol && i <= maxiter
        i = i+1; % increase the number of iterations
        Ja = J(xold); % calculate the Jacobian at {xi} (i.e. xold)
        xnew = xold - inv(Ja)*f(xold); % calculate {xi+1} (i.e. xnew)
        eps = max(abs((xnew-xold)./xnew)); % determine the updated epsilon
        xold = xnew;
    end
    if eps > tol
        disp('Maxium number of iterations reached, cannot find solution')
        disp('Change the initial guesses')
        sol = [];
        iter = maxiter;
    else
        sol = xnew;
        iter = i;
        fprintf("The number of iterations for the analytical method are %2.0f",iter)
    end
end


function [sol,iter] = newtonFD(f,J,x0)
    % Solve a system of nonlinear equations f(x) = 0, 
    %Jac is the Jacobian 
    % x0 is a vector of the initial guesses
    maxiter = 100; % maximum number of iterations before exiting the program
    tol = 1e-4; % tolerance for convergence
    eps = 1; % initial epsilon, any value OK as long as greater than the tolerance
    xold = x0; % initial guess
    i = 0; % initialize the number of iterations
    delta = 0.1;
    while eps > tol && i <= maxiter
        i = i+1; % increase the number of iterations
        Ja = J(f,xold,delta); % calculate the Jacobian at {xi} (i.e. xold)
        xnew = xold - inv(Ja)*f(xold); % calculate {xi+1} (i.e. xnew)
        eps = max(abs((xnew-xold)./xnew)); % determine the updated epsilon
        xold = xnew;
    end
    if eps > tol
        disp('Maxium number of iterations reached, cannot find solution')
        disp('Change the initial guesses')
        sol = [];
        iter = maxiter;
    else
        sol = xnew;
        iter = i;
        fprintf("The number of iterations for the finite difference method are %2.0f",iter)
    end
end
%The analytical method is 2 iterations faster than the finite difference method
