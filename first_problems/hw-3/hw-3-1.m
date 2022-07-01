# Ньютон-Котес
# Вариант 12
# y(x) = sin(x) + x
# (0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi) среднее значение x=0.4*pi
clear all;
format long;

dx = (2*pi)/10;
x = [0:dx:pi];
y = sin(x) + x;
n = 6;

for i = 1:1:n
    y(i) = sin(x(i)) + x(i);
end

y0=dx * sum(y(1:1:(n - 1)))
y1=(dx/2) * (y(1) + y(n) + 2 * sum(y(2:1:(n - 1))))
y2=(dx/3) * (y(1) + y(n) + 4 * sum(y(2:2:n-1)) + 2 * sum(y(3:2:n-2)))
y3=((3 * dx)/8) * (y(1) + y(n) - sum(y(4:3:n-3)) + 3 * sum(y(2:1:(n-1))))
y4=((2 * dx)/45) * (7 * (y(1) + y(n)) + 12*sum(y(2:4:n-2)) + 14 * sum(y(4:4:n-4)) + 32 * sum(y(1:2:n-1)))