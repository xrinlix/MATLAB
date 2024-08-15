%Rohit Thirumala
%21036098

clc
clear all

[solution,iteration]=newton_conc(@myfun_conc,@myjac_conc,[1;1])

x1 = solution(1);
x2 = solution(2);

K1 = 0.0004;
K2 = 0.037;
ca0 = 50;
cb0 = 20;
cc0 = 5;
cd0 = 10;
ca = ca0 - 2*x1 - x2;
cb = cb0 - x1;
cc = cc0 +x1 + x2;
cd = cd0 -x2;

fprintf("The values are: x1 = %g, x2 = %g, ca = %g, cb = %g, cc = %g, cd = %g", x1, x2, ca, cb, cc, cd)
fsolve(@myfun_conc,[1;1])

%The values x1,x2 are extremely close in value
function [jacobian] = myjac_conc(functions,vectorx, delta)
    jacobian_column1 = (functions([vectorx(1)*(1+delta),vectorx(2)])-functions([vectorx(1),vectorx(2)]))/(delta*vectorx(1));
    jacobian_column2 = (functions([vectorx(1),vectorx(2)*(1+delta)])-functions([vectorx(1),vectorx(2)]))/(delta*vectorx(2));
    jacobian = [jacobian_column1,jacobian_column2];
end



function [sol,iter] = newton_conc(f,J,x0)
    % Solve a system of nonlinear equations f(x) = 0, 
    %Jac is the Jacobian 
    % x0 is a vector of the initial guesses
    maxiter = 100; % maximum number of iterations before exiting the program
    tol = 1e-4; % tolerance for convergence
    eps = 1; % initial epsilon, any value OK as long as greater than the tolerance
    xold = x0; % initial guess
    i = 0; % initialize the number of iterations
    delta = 1*(10^-6);
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



function f = myfun_conc(x)
    K1 = 0.0004;
    K2 = 0.037;
    ca0 = 50;
    cb0 = 20;
    cc0 = 5;
    cd0 = 10;
    
    ca = ca0 - 2*x(1) - x(2);
    cb = cb0 - x(1);
    cc = cc0 +x(1) + x(2);
    cd = cd0 -x(2);
    
    f1 = K1- cc/((ca^2)*cb);
    f2 = K2 - cc/(ca*cd);
    
    f = [f1;f2];

end
