%Rohit Thirumala
%21036098
clc
clear all
x = [0.1,0.2,0.4,0.6,0.9,1.3,1.5,1.7,1.8];
y = [0.75,1.25,1.45,1.25,0.85,0.55,0.35,0.28,0.18];


yreg = log(y./x);

p = polyfit(x,yreg,1);
a = exp(p(2));
b = p(1);

polyval(p,x);
l = log(a+b*x);
figure(1)

plot(x,yreg,"*",x,polyval(p,x))
title("Linear fit - Data")
legend("Data", "Linear Fit")
figure(2)

plot(x,y,"*", x,a.*x.*exp(b.*x))
title("Original data - Equation")
legend("Data","Equation")