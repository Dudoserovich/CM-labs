% Построить интерполяционный полином по значениям x и y(x) различными методами (прямой метод, Лагранжа, Ньютона).
% Построить интерполяционные сплайны второй степени. Методом наименьших квадратов построить аппроксимирующую функцию
% первой, второй и третьей степени.

% МЕТОД ЛАГРАНЖА

# Вариант 12
# y(x) = sin(x) + x
# (0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi) среднее значение x=0.4*pi

clear all;

function [l] = lagrang(x, X, Y)
    n = size(X);
    l = 0;
    for k = 1:1:n
        mult = 1;
        for i = 1:1:(k-1)
            mult = mult.*(x-X(i))./(X(k)-X(i));
        end
        for i = (k+1):1:n
            mult = mult.*(x-X(i))./(X(k)-X(i));
        end
        l = l + mult.*Y(k);
    end
end

% Входные данные
x = [0, 0.2 * pi, 0.4 * pi, 0.6 * pi, 0.8 * pi, pi]';
n = length(x);
y = [];
for i = 1:1:n
    y(i) = sin(x(i))+x(i);
end

x2 = [x(1):0.01:x(n)]';
y2 = [];
for i = 1:1:length(x2)
    y2(i) = lagrang(x2(i), x, y);
end

[x3, y3] = least_square_method(1, x', y);
[x4, y4] = least_square_method(2, x', y);
[x5, y5] = least_square_method(3, x', y);

plot(x, y, 's', 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
hold on;
plot(x2, y2, 'b-', 'LineWidth', 2);
plot(x3, y3, 'g', 'LineWidth', 3);
plot(x4, y4, 'k:', 'LineWidth', 4);
plot(x5, y5, 'r- -', 'LineWidth', 2)
hold off;
grid on;
h = legend({'y(x)=sin(x)+x', 'Direct method of action', 'approx(1)', 'approx(2)', 'approx(3)'},
 'location', 'northeastoutside');
title('Direct method of action');
set(h, "fontsize", 12);