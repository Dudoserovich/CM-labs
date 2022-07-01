# Полином Ньютона
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

dy(1:1:n, 1:1:n) = 0;
dy(:,1) = y;
for j = 2:1:n
  for i = 1:1:(n-j+1)
    dy(i,j)=(dy(i+1,j-1)-dy(i,j-1));
  end;
end;

dy;

y = 1/dx*(dy(1,2) + 3*dy(1,3)/2 + 2*dy(1,4)/6 - 2*dy(1,5)/24)