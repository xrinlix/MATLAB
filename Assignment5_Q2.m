%Rohit Thirumala
%21036098
clc
clear all

%Set #3 cannot be solved using the Gauss-Seidel method since the diagonal
%is less than the sum of the row, (not a diagonally dominant matrix)

%Set number 3
A = [9 3 1
    -6 0 8
    2 5 -1]
b = [13
    2
    6]
x = [A,b]

val = Inf;
max_iter=1000; 
x=linspace(0,0,length(A))';
n=size(x,1);
tol=1e-3; 
iter=0;
while val>tol && iter<max_iter
    x_old=x;
    for i=1:n
        initial_guess=0;
        for j=1:i-1
            initial_guess=initial_guess+A(i,j)*x(j);
        end
        for j=i+1:n
            initial_guess=initial_guess+A(i,j)*x_old(j);
        end
        x(i)=(1/A(i,i))*(b(i)-initial_guess);
    end
    iter=iter+1;
    val=norm(x_old-x);
    
end
fprintf('Solution of the system is : %2.0f\n',x);

%Clearly the GS method is not converging since after 349 iterations the
%Result is Not a number(NaN)