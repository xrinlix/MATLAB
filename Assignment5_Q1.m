%Rohit Thirumala
%21036098
clc
clear all
f = @(x) (4-x)*exp(-0.5*x)-2;
fprime = @(x) (((x-6)*exp(-x/2))/2);
x0 = 2;
[sol,iter] = my_NR(f,2,fprime);
fprintf("for an intial guess of %2.0f, root is %2f, which is found after %2.0f, iterations \n ",x0, sol, iter)
x0 = 4;
[sol,iter] = my_NR(f,4,fprime);
fprintf("for an intial guess of %2.0f, root is %2f, which is found after %2.0f, iterations \n ",x0, sol, iter)
x0 = 6;
[sol,iter] = my_NR(f,6,fprime);
fprintf("for an intial guess of %2.0f, root is %2f, which is found after %2.0f, iterations \n ",x0, sol, iter)

%The newton raphson method works for initial guess's of x0 = 2, and x0 = 4
%It fails to work for x0 = 6 however.
function [sol,iter] = my_NR(f,x0,fprime)

    maxiter = 100; 
    tol = 1e-6; 
    eps = 1;
    xold = x0; 
    i = 0; 
 
    while eps > tol && i <= maxiter
        i = i+1; 
        xnew = xold - (f(xold))/fprime(xold); 
        eps = (abs((xnew-xold)/xnew))*100; 
        xold = xnew;
    end
    if eps > tol
        disp('Maxium number of iterations reached, cannot find solution');
        disp('Change the initial guesses');
        sol = [];
        iter = maxiter;
    else
        sol = xnew;
        iter = i;
        fprintf("The number of iterations for the newton raphson method are %2.0f \n",iter)
    end
end

