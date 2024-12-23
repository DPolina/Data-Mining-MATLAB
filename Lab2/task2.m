close all
clear all

    % Создание файла:

% Задаем количество точек
N = 10;

% Генерируем тестовые данные
x = linspace(0, 10, N);
y = sin(x) + randn(N, 1) * 0.2;
sigma = 0.2 * ones(N, 1);

% Задаем полный путь к папке
folder = 'C:\Users\dapol\OneDrive\Документы\MATLAB\LAB2';

% Создаем полный путь к файлу
fullpath = fullfile(folder, 'data.txt');

% Открываем файл для записи
fid = fopen('data.txt', 'w');

% Записываем данные в файл
for k = 1:N
    fprintf(fid, '%f,%f,%f\n', x(k), y(k), sigma(k));
end

% Закрываем файл
fclose(fid);

    % Загрузка данных:

% Задаем имя файла
filename = 'data.txt';

% Открываем файл для чтения
fid = fopen(filename, 'r');

% Читаем данные из файла
data = fscanf(fid, '%f,%f,%f');

% Закрываем файл
fclose(fid);

% Разбиваем данные на массивы
x = data(1:3:end);
y = data(2:3:end);
sigma = data(3:3:end);

    % Построение графика:

% Создаем фигуру
figure;

% Задаем заголовок
title('Экспериментальные данные');

% Отрисовываем точки с ошибками
errorbar(x, y, sigma, 'o');

% Подписываем оси
xlabel('x');
ylabel('y');

% Настраиваем сетку
grid on;

% Увеличиваем размер шрифта
set(gca, 'FontSize', 14);
