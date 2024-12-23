clc
close all
clear all

file = fopen('data6.txt');
array = fscanf(file, '%f', [3 inf])';
fclose(file);
N = size(array);
N = N(1);
x = array(:, 1);
y = array(:, 2);

figure;
plot(x, y);
hold on;
grid on;
title('y(x)')
xlabel('x')
ylabel('y')

disp('1 - ф-я Вейбулла, 2 - ф-я Хи2, 3 - ф-я Бета');
disp('-------------------------------------------------------');

q1 = [1 1]; % Начальное приближение
res1 = nlinfit(x, y, @weibull, q1); % Нелинейное регрессионное моделирование для функции Вейбулла
f1 = weibull(res1, x); % Вычисляет значения функции Вейбулла для заданных значений x
disp('Аппроксимированные начальные значения для ф-ции Вейбулла:');
disp(res1);

q2 = 2;
res2 = nlinfit(x, y, @Xi2, q2);
f2 = Xi2(res2, x);
disp('Аппроксимированные начальные значения для ф-ции Хи2:');
disp(res2);

q3 = [0.14 0.35];
res3 = nlinfit(x, y, @beta, q3);
f3 = beta(res3, x);
disp('Аппроксимированные начальные значения для ф-ции Бета:');
disp(res3);

figure;
hold on;
grid on;
plot(x, y, '.');
plot(x, f1, 'r'); % weibull
plot(x, f2, 'g'); % Xi2
plot(x, f3, 'b'); % Beta
title('f1, f2, f3')
xlabel('x')
ylabel('y')

% u1, u2, u3 - степени свободы
% N - число параметров - 1
u1 = N - 2 - 1;   
u2 = N - 1 - 1;
u3 = N - 2 - 1;
x_kv_1 = 0;
x_kv_2 = 0;
x_kv_3 = 0;

figure;
plot(x, y);
hold on;
grid on;
plot(x, f3, 'm');
title('Наилучшая аппроксимация Бета-функцией')
xlabel('x')
ylabel('y')

% Среднеквадратическое отклонение
sygm = array(:,3);

% Сумма квадратов взвешенных отклонений значений
for i = 1:N
    x_kv_1 = x_kv_1 + ((f1(i) - y(i))/sygm(i))^2;
    x_kv_2 = x_kv_2 + ((f2(i) - y(i))/sygm(i))^2;
    x_kv_3 = x_kv_3 + ((f3(i) - y(i))/sygm(i))^2;
end;

% Нормировка
x_kv_1 = x_kv_1/u1;
x_kv_2 = x_kv_2/u2;
x_kv_3 = x_kv_3/u3;

disp('-------------------------------------------------------');
%disp(newline);

disp('Нормированный критерий Xi2 для фунуции Вейбулла:');
disp(x_kv_1);

disp('Нормированный критерий Xi2 для фунуции Хи квадрат:');
disp(x_kv_2);

disp('Нормированный критерий Xi2 для фунуции Бета:');
disp(x_kv_3);

% Взвешенные остатки
R1 = zeros(1, N);
R2 = zeros(1, N);
R3 = zeros(1, N);

for i = 1:N
    R1(i) = (y(i) - f1(i))/sygm(i);
    R2(i) = (y(i) - f2(i))/sygm(i);
    R3(i) = (y(i) - f3(i))/sygm(i);
end;

figure;
plot(R1);
title('Взвешенные остатки 1');
figure;
plot(R2);
title('Взвешенные остатки 2');
figure;
plot(R3);
title('Взвешенные остатки 3');
histogram(R3);

[res3, r, J, covb] = nlinfit(x, y, @beta, q3);

% Автокорреляционная функция взвешенных остатков
A1 = zeros(N/2);
A2 = zeros(N/2);
A3 = zeros(N/2);

% Вычисление автокорреляционной функции:
for k = 1:(N/2)
    S1 = 0;
    for i = 1:(N-k+1)
        S1 = S1 + R1(i)*R1(i+k-1);
    end;
    S2 = 0;
    for i = 1:N
        S2 = S2 + (R1(i))^2;
    end;
    S2 = S2/N;
    A1(k) = (1/(N-k+1))*S1/S2;
    
    S1 = 0;
    for i = 1:(N-k+1)
        S1 = S1 + R2(i)*R2(i+k-1);
    end;
    S2 = 0;
    for i = 1:N
        S2 = S2 + (R2(i))^2;
    end;
    S2 = S2/N;
    A2(k) = (1/(N-k+1))*S1/S2;
    
    S1 = 0;
    for i = 1:(N-k+1)
        S1 = S1 + R3(i)*R3(i+k-1);
    end;
    S2 = 0;
    for i = 1:N
        S2 = S2 + (R3(i))^2;
    end;
    S2 = S2/N;
    A3(k) = (1/(N-k+1))*S1/S2;
end;    

k = 1:1:N/2;
figure;
plot(k, A1);
title('Автокорреляционная функция 1');
figure;
plot(k, A2);
title('Автокорреляционная функция 2');
figure;
plot(k, A3);
title('Автокорреляционная функция 3');

% Доверительный интервал - 1
[res3, r, J, covb] = nlinfit(x, y, @beta, q3);

disp('-------------------------------------------------------');

% Доверительный интервал по t-статистике

alpha = 0.32; % Уровень значимости

% Проверка корректности
if ~isnumeric(x(2)) || ~isnumeric(q3(2))
    error('x(2) or q3(2) is not a number!');
end

df = max(length(x) - length(q3) - 1, 1); % Степень свободы

t = tinv(1-alpha/2, df); % Критическое значение t-статистики

se = sqrt(diag(covb)); % Стандартная ошибка оценки

lower_bound = res3 - t * se'; % Нижняя граница
upper_bound = res3 + t * se'; % Верхняя граница

% Вывод информации
disp('Теоретический доверительный интервал:');
disp([lower_bound, upper_bound]);





