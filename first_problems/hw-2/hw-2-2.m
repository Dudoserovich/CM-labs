# Разностные схемы
# Вариант 12
# y(x) = sin(x) + x
# (0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi) среднее значение x=0.4*pi

clear all;
format long;

dx = (2*pi)/10;

#f'(x) = (f(x+dx) - f(x))/dx + O(dx)
y1 = ((sin(0.4*pi + dx) + 0.4*pi + dx) - (sin(0.4*pi)+0.4*pi))/dx
#f'(x) = (f(x) - f(x-dx))/dx +O(dx)
y2 = ((sin(0.4*pi)+0.4*pi) - (sin(0.4*pi - dx) + (0.4*pi - dx)))/dx
#f'(x) = (f(x+dx) - f(x-dx))/2dx + O(dx^2)
y3 = ((sin(0.4*pi + dx) + 0.4*pi + dx) - (sin(0.4*pi - dx) + (0.4*pi - dx)))/(2*dx)
# f''(x) = (f(x+dx) - 2f(x) + f(x-dx))/(dx)^2
y4 = ((sin(0.4*pi + dx) + 0.4*pi + dx) - 2*((sin(0.4*pi) + 0.4*pi)) + (sin(0.4*pi - dx) + (0.4*pi - dx)))/(dx*dx)