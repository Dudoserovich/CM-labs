% Построить интерполяционный полином по значениям x и y(x) различными методами (прямой метод, Лагранжа, Ньютона).
% Построить интерполяционные сплайны второй степени. Методом наименьших квадратов построить аппроксимирующую функцию
% первой, второй и третьей степени.

% ИНТЕРПОЛЯЦИОННЫЙ ПОЛИНОМ НЬЮТОНА

# Вариант 12
# y(x) = sin(x) + x
# (0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi) среднее значение x=0.4*pi

clear all;

function s = newton (x, y, t)
    n = length(x) ; % Вычисление количества точек в массивах x и y
    for i = 1:n % Запись в первый столбец C разделённых разностей y
        C(i, 1) = y(i);
    end
    for i = 2:n % Формирование матрицы разделённых разностей
        for j = 2:n
            if (i < j)
                C(i, j) = 0;
            else
                C(i, j) = (C(i, j - 1) - C(j - 1, j - 1)) / (x(i) - x(j - 1));
            end
        end
    end
    for i = 1:n % Формирование массива коэффициентов полинома Ньютона
        A(i) = C(i, i);
    end

    s = 0; % Расчёт значения полинома в точке t по формуле (12.5)
    for i = 1:n
         p = 1;
         for j = 1:i-1
             p = p * (t - x(j));
         end
         s = s + A(i) * p;
    end
end

function [sum] = PolinNewton(x, X, Y)
    n = size(X);
    dx = 2;
    f_razn(1:1:n,1:1:n)=0;
    f_razn(:,1)=X;
    f_razn(:,2)=Y;

    for i=1:1:(n-1)
        for j=1:1:(n-i)
            f_razn(j,i+2)=f_razn(j+1,i+1)-f_razn(j,i+1);
        end
    end

    sum = 0;
    for i=1:1:n
        pr = 1;
        for j=1:1:(i-1)
            pr=pr.*(x-X(j));
        end
        sum = sum + (f_razn(1, i+1) / (factorial(i-1) * dx^(i-1))) * pr;
    end
end

x = [0, 0.2 * pi, 0.4 * pi, 0.6 * pi, 0.8 * pi, pi]';
n = length(x);
y = [];
for i = 1:1:n
    y = [y sin(x(i)) + x(i)];
end

x2 = [x(1):0.01:x(n)]';
y2 = [];
for t = x(1):0.01:x(n)
    % не работает
    % y2 = [y2 PolinNewton(t, x, y)];
    y2 = [y2 newton(x, y, t)];
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