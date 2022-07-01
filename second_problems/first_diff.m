clear all;

% Хмелевский Е.Д., Б9119-02.03.03техпро
% diff №1(12): y' = 2x^2+y, y(0) = 2

% Полученная после дифференцирования функции является пораболой ветвями в верх,
% на числе 10 она возрастает уже очень сильно,
% поэтому на [0..10] мы видим график не особо напоминающий пораболу

Dif1 = @(X, Y)(2*X^2+Y);
Dif1Init = 2;
dF1dx = @(X, Y)(4*X);
dF1dy = @(X, Y)(1);
dx = 0.1;
C = 6;
points = [0:dx:10];
Dif1Solution = @(X)(2*X.^2 - 4*X + C*e.^X - 4);

function res = Yavniy(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  for I = 2:Size
    res(I) = res(I - 1) + Func(points(I - 1), res(I - 1)) * dx;
  end
end

function res = NeYavniy(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  for I = 2:Size
    res(I) = res(I - 1) + Func(points(I - 1), res(I - 1)) * dx;
    res(I) = res(I - 1) + Func(points(I - 1), res(I)) * dx;
  end
end

function res = Trepez(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  for I = 2:Size
    res(I) = res(I - 1) + dx/2*(Func(points(I - 1), res(I - 1)) + Func(points(I), res(I - 1) + Func(points(I - 1), res(I - 1)) * dx));
  end
end

function res = Utoncheniy(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  res(2) = res(1) + Func(points(1), res(1)) * dx;
  for I = 3:Size
    res(I) = res(I - 2) + 2*(Func(points(I - 1), res(I - 1)))*dx;
  end
end

function res = Ispravleniy(Func, points, Y, dFdx, dFdy)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  for I = 2:Size
    X = points(I - 1);
    Y = res(I - 1);
    res(I) = Y + Func(X, Y)*dx + dx*dx/2*(dFdx(X, Y) + dFdy(X, Y)*Func(X, Y));
  end
end

function res = RungeKutta(Func, points, Y, Alpha)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  for I = 2:Size
    FirstTerm = (1 - Alpha)*Func(points(I - 1), res(I - 1));
    SecondTerm = Alpha*Func(points(I - 1) + dx/(2*Alpha), points(I - 1) + dx/(2*Alpha)*Func(points(I - 1), res(I - 1)));
    res(I) = res(I - 1) + dx*(FirstTerm + SecondTerm);
  end
end

function res = AdamsOrder2(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  res(2) = res(1) + Func(points(1), res(1)) * dx;
  for I = 3:Size
    res(I) = res(I-1) + dx*(3/2*Func(points(I-1), res(I-1)) - 1/2*Func(points(I-2), res(I-2)));
  end
end

function res = AdamsOrder3(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  res(2) = res(1) + Func(points(1), res(1)) * dx;
  res(3) = res(2) + dx*(3/2*Func(points(2), res(2)) - 1/2*Func(points(1), res(1)));
  for I = 4:Size
    res(I) = res(I-1) + dx*(23/12*Func(points(I-1), res(I-1)) - 16/12*Func(points(I-2), res(I-2)) + 5/12*Func(points(I-3), res(I-3)));
  end
end

function res = AdamsOrder4(Func, points, Y)
  Size = size(points)(2);
  dx = points(2) - points(1);
  res(Size) = 0;
  res(1) = Y;
  res(2) = res(1) + Func(points(1), res(1)) * dx;
  res(3) = res(2) + dx*(3/2*Func(points(2), res(2)) - 1/2*Func(points(1), res(1)));
  res(4) = res(3) + dx*(23/12*Func(points(3), res(3)) - 16/12*Func(points(2), res(2)) + 5/12*Func(points(1), res(1)));
  for I = 5:Size
    res(I) = res(I-1) + dx*(55/24*Func(points(I-1), res(I-1)) - 59/24*Func(points(I-2), res(I-2)) + 37/24*Func(points(I-3), res(I-3)) - 9/24*Func(points(I-4), res(I-4)));
  end
end

printf("First eq, Euler\n")
printf("Explicit Euler ");
mean(abs(Dif1Solution(points) - Yavniy(Dif1, points, Dif1Init)))
printf("Implicit Euler ");
mean(abs(Dif1Solution(points) - NeYavniy(Dif1, points, Dif1Init)))
printf("Trapezoid Euler ");
mean(abs(Dif1Solution(points) - Trepez(Dif1, points, Dif1Init)))
printf("Refined Euler ");
mean(abs(Dif1Solution(points) - Utoncheniy(Dif1, points, Dif1Init)))
printf("Corrected Euler ");
mean(abs(Dif1Solution(points) - Ispravleniy(Dif1, points, Dif1Init, dF1dx, dF1dy)))

figure(1);
plot(points, Yavniy(Dif1, points, Dif1Init), "linewidth", 2,
     points, NeYavniy(Dif1, points, Dif1Init), "linewidth", 2,
     points, Trepez(Dif1, points, Dif1Init), "linewidth", 2,
     points, Utoncheniy(Dif1, points, Dif1Init), "linewidth", 2,
     points, Ispravleniy(Dif1, points, Dif1Init, dF1dx, dF1dy), "linewidth", 2,
     points, Dif1Solution(points), "linewidth", 2
     );
legend("Explicit Euler", "Implicit Euler", "Trapezoid Euler", "Refined Euler", "Corrected Euler", "Real Solution");

printf("\n")
printf("First eq, RungeKutta\n")

printf("Alpha = 0.2 ")
mean(abs(Dif1Solution(points) - RungeKutta(Dif1, points, Dif1Init, 0.2)))
printf("Alpha = 0.1 ")
mean(abs(Dif1Solution(points) - RungeKutta(Dif1, points, Dif1Init, 0.1)))
printf("Alpha = 0.075 ")
mean(abs(Dif1Solution(points) - RungeKutta(Dif1, points, Dif1Init, 0.075)))
printf("Alpha = 0.05 ")
mean(abs(Dif1Solution(points) - RungeKutta(Dif1, points, Dif1Init, 0.05)))

figure(2);
plot(points, RungeKutta(Dif1, points, Dif1Init, 0.2), "linewidth", 2,
  points, RungeKutta(Dif1, points, Dif1Init, 0.1), "linewidth", 2,
  points, RungeKutta(Dif1, points, Dif1Init, 0.075), "linewidth", 2,
  points, RungeKutta(Dif1, points, Dif1Init, 0.05), "linewidth", 2,
  points, Dif1Solution(points), "linewidth", 2);
legend("0.2", "0.1", "0.075", "0.05", "RealSolution");

printf("\n");
printf("First eq, Adams-Bashforth\n");
printf("Second order ")
mean(abs(Dif1Solution(points) - AdamsOrder2(Dif1, points, Dif1Init)))
printf("Third order ")
mean(abs(Dif1Solution(points) - AdamsOrder3(Dif1, points, Dif1Init)))
printf("Fourth order ")
mean(abs(Dif1Solution(points) - AdamsOrder4(Dif1, points, Dif1Init)))