clear all
close all

disp("Задание 13");
% Создайте файл test_cycle.m, в котором сформируйте матрицу Q произвольного
% размера и заполните ее случайными целыми числами. Выведите в главном окне
% MATLAB (Command Window) матрицу Q. (добавьте в код команду disp(Q)).

n=5; % Количество строк
m=5; % Количество столбцов

% Формируем нулевую матрицу Q размером n x m
Q = zeros(n,m);

% В цикле заполняем матрицу Q новыми значениями
for k = 1:n
 for j = 1:m
    Q(k,j) = round(10*rand+5);
 end
end

disp(Q);

disp("Задание 14"); 
% В файле test_cycle.m реализуйте цикл, позволяющий вычислять сумму 
% элементов матрицы Q. Весь необходимый код напишите в файле test_cycle.m.

sum = 0;
for k = 1:n
    for j = 1:m
        sum = sum + Q(k,j);
    end
end

disp(['Sum = ', num2str(sum)]);