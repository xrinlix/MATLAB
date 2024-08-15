%Rohit Thirumala
%21036098
clc
clear all

%System 1
f = @(x) x.^2 - 10;
g = @(x) 15 - x;

x = linspace(-10,10);
plot(x,f(x))
hold on 
plot(x,g(x))

%System 2

w = @(x) (26-x.^2).^0.5
v = @(x) ((100 - 3.*x.^2)/25).^0.5

plot(x,w(x))
hold on 
plot(x,v(x))

root1 = newton(@(x) myfun(x),@(x) myjac(x),[-1;1])
root2 = newton(@(x) myfun(x),@(x) myjac(x),[1;1])
root3 = newton(@(x) myfun(x),@(x) myjac(x),[-1;-1])
root4 = newton(@(x) myfun(x),@(x) myjac(x),[1;-1])

function F = myfun(x)

	f1 = x(1)^2 + x(2)^2 - 26;
	f2 = 3*x(1)^2+ 25*x(2)^2 - 100;
	F = [f1;f2];
end


function J = myjac(x)

    J=[2*x(1) 2*x(2);6*x(1) 50*x(1)];
	
end

function [sol,iter] = newton(f,J,x0)
    
    maxiter = 10;
    tol = 0.001; 
    eps = 1; 
    xold = x0; 
    i = 0; 
 
    while eps > tol && i <= maxiter
        i = i+1; 
        Ja = J(xold);
        xnew = xold - inv(Ja)*f(xold); 
        eps = max(abs((xnew-xold)./xnew)); 
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
    end
end
