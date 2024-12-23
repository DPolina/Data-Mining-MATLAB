clc
close all
clear all

% Вариант 6
% 1 - Eвклидово расстояние, 3 - метрика города, 5 - расстояние Минковского.
% a - Метод ближнего соседа, d - центроидный метод, e - метод медианной
% связи.

file = fopen('data6.txt');
X = fscanf(file, '%f', [2 inf])';
fclose(file);

x = X(:,1);
y = X(:,2);
N = size(X);
K = N(2);
N = N(1);

figure;
scatter(x, y, 'filled');
title('Экспериментальные данные');
xlabel('X2');
ylabel('X1');

% Расстояния между объектами

% Функция pdist в MATLAB вычисляет парные расстояния между точками в наборе данных.
% Она принимает один аргумент:
% X: матрица, в которой каждый столбец представляет собой точку данных.
% Функция pdist возвращает вектор D, в котором каждый элемент представляет 
% собой расстояние между двумя точками данных. Порядок элементов в D 
% соответствует индексам строк в X.
d_evkl = pdist(X, 'euclidean');
d_town = pdist(X, 'cityblock'); 
% d_town_matrix = squareform(d_town, 'tomatrix') % В виде матрицы
d_mink = pdist(X, 'minkowski', 4); % 4 для Минковского

% Связывание или группировка объектов в бинарные иерархические деревья

% Функция linkage в MATLAB используется для создания иерархической кластерной структуры из набора данных.
% Она принимает два аргумента:
% - X: матрица, в которой каждый столбец представляет собой точку данных.
% - method: метод кластеризации, который может быть одним из следующих:
%   - 'single': Метод ближнего соседа.
%   - 'centroid': Этот метод объединяет два кластера, основываясь на расстоянии 
% между центроидами этих кластеров. Центроид - это средняя точка всех точек в кластере.
%   - 'median': Этот метод объединяет два кластера, основываясь на медианном расстоянии между 
% точками в этих кластерах. Медиана - это значение, которое делит набор данных 
% на две половины с равным числом точек.
% Функция linkage возвращает матрицу связей Z, где каждая строка представляет собой кластерное слияние.
link_evkl_sing = linkage(d_evkl, 'single');
link_evkl_centr = linkage(d_evkl, 'centroid');
link_evkl_media = linkage(d_evkl, 'median');
link_town_sing = linkage(d_town, 'single');
link_town_centr = linkage(d_town,'centroid');
link_town_media = linkage(d_town, 'median');
link_mink_sing = linkage(d_mink, 'single');
link_mink_centr = linkage(d_mink, 'centroid');
link_mink_media = linkage(d_mink, 'median');

% Анализ качества кластеризации:
disp('Кофенетические корреляционные коэффициенты');
coph_corr_coeff = zeros(3);
coph_corr_coeff(1,1) = cophenet(link_evkl_sing, d_evkl);
coph_corr_coeff(2,1) = cophenet(link_evkl_centr, d_evkl);
coph_corr_coeff(3,1) = cophenet(link_evkl_media, d_evkl); 

coph_corr_coeff(1,2) = cophenet(link_town_sing, d_town);
coph_corr_coeff(2,2) = cophenet(link_town_centr, d_town); % Max effective
coph_corr_coeff(3,2) = cophenet(link_town_media, d_town);

coph_corr_coeff(1,3) = cophenet(link_mink_sing, d_mink); % Min effective
coph_corr_coeff(2,3) = cophenet(link_mink_centr, d_mink); 
coph_corr_coeff(3,3) = cophenet(link_mink_media, d_mink);

coph_corr_coeff

disp('Самый эффективный: метрика города + центроидный метод');
max_effective = max(max(coph_corr_coeff))
disp('Самый неэффективный: расстояние Минковского + метод ближнего соседа');
min_effective = min(min(coph_corr_coeff))

figure;
dendrogram(link_town_centr);
title('Mетрика города + центроидный метод');

% Фиксированное число кластеров:
N_clust = 3;

clust = cluster(link_mink_media, 'maxclust', N_clust); % construct clusters from linkages

figure;
gscatter(x, y, clust); % scatter plot by group
hold on;

% Центры кластеров
clust_centr = zeros(N_clust, K); % Матрица x и y координат центров
% Расстояния между центрами кластеров
dist_betw_clust = zeros(N_clust, N_clust);
% Радиусы кластеров
clust_rad = zeros(N_clust, 1);
% Дисперсии кластеров
clust_disp = zeros(N_clust, 1);

disp('Расстояния от элементов кластера до центра');
for k = 1:N_clust
    disp(['№', num2str(k)]);
    % Находим все объекты кластера
    obj = find(clust == k); % Какие по счёту из clust относятся к текущему кластеру
    N_obj = length(obj); % Их кол-во
    singl_clust = zeros(N_obj, K); % Координаты всех объектов конкретного кластера (на итерации, 
    % после вычисления очередного кластера предыдущее значение теряется за ненадобностью), матрица x и y
    
    % Находим центр и элементы кластера
    for i = 1:N_obj
        singl_clust(i,:) = X(obj(i),:);
        clust_centr(k,:) = clust_centr(k,:) + singl_clust(i,:);
    end;
    clust_centr(k,:) = clust_centr(k,:)/N_obj; 
    
    % Расстояния (геометрические) до центра от всех элементов кластера
    clust_dist_centr = zeros(N_obj, 1);
    for i = 1:N_obj
        for l = 1:K
            clust_dist_centr(i) = clust_dist_centr(i) + (singl_clust(i,l) - clust_centr(k,l))^2;
        end;
        % Дисперсия кластера
        clust_disp(k) = clust_disp(k) + clust_dist_centr(i);
        clust_dist_centr(i) = sqrt(clust_dist_centr(i));
    end;
    clust_disp(k) = clust_disp(k)/N_obj;
    clust_dist_centr
    % Радиус кластера 
    clust_rad(k) = max(clust_dist_centr);
end;

% Геометрическое расстояние между центрами
for i = 1:N_clust
    for j = i+1:N_clust
        for l = 1:K
            dist_betw_clust(i,j) = dist_betw_clust(i,j) + (clust_centr(i,l) - clust_centr(j,l))^2; %
        end;
        dist_betw_clust(i,j) = sqrt(dist_betw_clust(i,j));
    end;
end;
disp('Расстояния между центрами кластеров');
dist_betw_clust
disp('Дисперсии кластеров');
clust_disp
disp('Радиусы кластеров');
clust_rad

% Центры кластеров
scatter(clust_centr(:,1), clust_centr(:,2), 20, 'k', 'filled');
t = 0:pi/180:2*pi;
for i = 1:N_clust
    x = clust_rad(i)*cos(t) + clust_centr(i,1);
    y = clust_rad(i)*sin(t) + clust_centr(i,2);
    plot(x, y, '-'); % Рисует овальные границы
end;
title('Найденные кластеры');
xlabel('X2');
ylabel('X1');
