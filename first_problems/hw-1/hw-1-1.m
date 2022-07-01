% Построить интерполяционный полином по значениям x и y(x) различными методами (прямой метод, Лагранжа, Ньютона).
% Построить интерполяционные сплайны второй степени. Методом наименьших квадратов построить аппроксимирующую функцию
% первой, второй и третьей степени.

% ПРЯМОЙ МЕТОД ДЕЙСТВИЯ

# Вариант 12
# y(x) = sin(x) + x
# (0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi) среднее значение x=0.4*pi

clear all;

% Входные данные
x = [0, 0.2 * pi, 0.4 * pi, 0.6 * pi, 0.8 * pi, pi]';
n = length(x);
y = [];
for i = 1:1:n
    y = [y (sin(x(i))+x(i))];
end
y = y';

for i=1:1:n
    for j=1:1:n
        a(i, j) = x(i) ^ (j - 1);
    end
end

k = a ^ (-1) * y;
x2 = [x(1):0.01:x(n)]';
n2 = length(x2);
y2 = zeros(n2);

for i=1:1:n
    y2 = y2 + k(i).* x2.^(i - 1);
end

[x3, y3] = least_square_method(1, x, y);
[x4, y4] = least_square_method(2, x, y);
[x5, y5] = least_square_method(3, x, y);

plot(x, y, 's', 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
hold on;
plot(x2, y2, 'LineWidth', 1);
plot(x3, y3, 'LineWidth', 1);
plot(x4, y4, 'LineWidth', 2);
plot(x5, y5, 'LineWidth', 1)
hold off;
grid on;
h = legend({'y(x)=sin(x)+x', 'Direct method of action', 'approx(1)', 'approx(2)', 'approx(3)'},
 'location', 'northeastoutside');
title('Direct method of action');
set(h, "fontsize", 12);