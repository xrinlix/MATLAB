%Rohit Thirumala
%21036098
clc
clear all
x=[10 20 30 40 50 60 70 80];
y=[25 70 380 550 610 1220 830 1450];
[a, r2, syx] = linregr2(x,y)


y_main = a(1) + a(2)*linspace(10,80,8);
subplot(2,1,1)
plot(x,y_main)
title("Linear Fit Plot")
legend("x","function")
hold on
plot(x,y,"*")

subplot(2,1,2)

plot(x,y_main-y)
legend("Data points", "Fitted Values")
title("Residuals")


function [a, r2, syx] = linregr2(x,y)
    if length(x) ~= length(y)
        error("Lengths dont match")
    end
    n = length(x);
    a1 = (n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.^2)-(sum(x)^2));
    a0 = sum(y)./n - a1.*(sum(x)./n);
    r2 = (((n.*(sum(x.*y)) - sum(x).*sum(y)))/(sqrt(n.*sum(x.^2)-(sum(x).^2)).*sqrt(n.*sum(y.^2)-(sum(y).^2))))^2;
    Sr = sum((y-a0-a1.*x).^2);
    syx = sqrt(Sr./(n-2));
    a = [a0,a1];

end