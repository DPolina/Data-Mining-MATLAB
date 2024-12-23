close all
clear all

disp("Задание 11"); 
% Создайте файл test_my_func.m, в котором наберите код, демонстрирующий 
% работу функции my_func. Запустите программу и просмотрите результат ее 
% работы. Создайте заголовки гистограмм и подписи по осям с помощью функций 
% title, xlabel, ylabel, добавив их после функции bar.

% Задание диапазона изменения X
X_left=-2;
X_right=2;

% Задание диапазона изменения Y
Y_left=-3;
Y_right=3;

% Задание количества сгенерированных точек
N=1000;

% Вызов функции
[X,Y]=my_func(X_left, X_right, Y_left, Y_right, N);

% Построение графика функции
plot(X,Y,'*')

%------------------------

% Инициализация гистограммы
BinNumber=10;
k=0:BinNumber;

% Вычисление границ карманов на оси X
X_bins=X_left + k*(X_right - X_left)/BinNumber;

% Вычисление границ карманов на оси Y
Y_bins=Y_left + k*(Y_right - Y_left)/BinNumber;

% Вычисление гистограммы для X
 N_X = histc(X,X_bins);

% Вычисление гистограммы для Y
 N_Y = histc(Y,Y_bins);

%------------------------

figure;
bar(X_bins, N_X);

figure;
bar(Y_bins, N_Y);