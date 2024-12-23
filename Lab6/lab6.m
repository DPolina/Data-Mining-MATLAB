clc
close all
clear all

% Задание 1

file = fopen('Cohonen_Data_Lab6/Cohonen_data6.txt', 'r');
A = fscanf(file, '%g', [2 inf]);
A1 = A';
fclose(file);

X = A(1,:);
X = X';
Y = A(2,:);
Y = Y';

% Матрица max и min значений для наблюдений
R(1,1) = min(A1(:,1));
R(1,2) = max(A1(:,1));
R(2,1) = min(A1(:,2));
R(2,2) = max(A1(:,2));

% Функция создания соревнующегося слоя Кохонена
net = newc(R',4);
disp(net);

net.trainParam.epochs = 100;
net = train(net, A);

% Веса входов
celldisp(net.IW)
% Преобразоване весов из ячейки в матрицу
neuron_coords = cell2mat(net.IW)

% Классификация векторов входов по сети
Y = sim(net, A);
% Конвертация объектов в индексы нейронов 
Yc = vec2ind(Y);

X = A(1,:);
X = X';
Y = A(2,:);
Y = Y';
hold on
gscatter(X, Y, Yc);

% Задание 2

file = fopen('SOM_Data_Lab6/Learning_data6.txt', 'r');
Learning = fscanf(file, '%f', [8 inf]);
fclose(file);

Min = min(Learning, [], 2); % Минимальный среди объектов по i-ому признаку
Max = max(Learning, [], 2);
T = [Min(1) Max(1); Min(2) Max(2); Min(3) Max(3); Min(4) Max(4); 
    Min(5) Max(5); Min(6) Max(6); Min(7) Max(7); Min(8) Max(8)]

% Создать нейронную сеть (например, 2x2) на основе карт Кохонена
size1 = 2;
size2 = 2;
net = newsom(T, [size1 size2]); % Для создания самоорганизующейся карты Кохонена
net.trainParam.epochs = 100; % Один этап обучения определяется как однократное представление всех входных векторов в сеть
net = train(net, Learning); % Обучение сети соответсвующим набором

% Классифицировать исходные объекты в кластеры с использованием разработанной нейронной сети
W = sim(net, Learning); % Simulate dynamic system, классификация векторов входов по разработанной нейросети
Clusters = vec2ind(W); % Конвертирует классифицированные объекты в индексы нейронов, какой объект к какому кластеру

file = fopen('SOM_Data_Lab6/PCA_data6.txt', 'r');
PCA = fscanf(file, '%f', [2 inf]); % Исходные данные в координатах первых двух главных компонент
fclose(file);

% Сгруппировать объекты в кластеры
Obj1 = find(Clusters==1);
s1 = length(Obj1);
Obj2 = find(Clusters==2);
s2 = length(Obj2);
Obj3 = find(Clusters==3);
s3 = length(Obj3);
Obj4 = find(Clusters==4);
s4 = length(Obj4);
C1 = zeros(s1, 2); C2 = zeros(s2, 2); C3 = zeros(s3, 2); C4 = zeros(s4, 2);
clust1 = zeros(s1, 8); clust2 = zeros(s2, 8); clust3 = zeros(s3, 8); clust4 = zeros(s4, 8);

for i = 1:s1 % Цикл по всем элементам i-го кластера
        C1(i,:) = PCA(:,Obj1(i)); % Переписать в новую матрицу только те объекты которые в p-ом кластере 
        clust1(i,:) = Learning(:,Obj1(i));
end 
for i = 1:s2 
        C2(i,:) = PCA(:,Obj2(i)); 
        clust2(i,:) = Learning(:,Obj2(i));
end 
for i = 1:s3 
        C3(i,:) = PCA(:,Obj3(i)); 
        clust3(i,:) = Learning(:,Obj3(i));
end 
for i = 1:s4 
        C4(i,:) = PCA(:,Obj4(i)); 
        clust4(i,:) = Learning(:,Obj4(i));
end 

% Центры кластеров данных из PCA
M_PCA(1,:) = mean(C1);
M_PCA(2,:) = mean(C2);
M_PCA(3,:) = mean(C3);
M_PCA(4,:) = mean(C4);

% Отображение центров кластеров
figure
gscatter(PCA(1,:), PCA(2,:), Clusters);
hold on
plot(M_PCA(:,1), M_PCA(:,2), '*');
hold off

% Cредние значения по каждому признаку в кластерах
M(1,:) = mean(clust1);
M(2,:) = mean(clust2);
M(3,:) = mean(clust3);
M(4,:) = mean(clust4);
M

% График средних значений признаков объектов попавших в кластер
figure
plot(M(1, :), 'r');
axis([1 8 1 10]); % Границы по осям [x1 x2 y1 y2] 
hold on
plot(M(2, :), 'g');
plot(M(3, :), 'c');
plot(M(4, :), 'b');
title('График средних значений признаков');
