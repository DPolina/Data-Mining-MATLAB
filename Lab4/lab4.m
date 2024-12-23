clc
close all
clear all

% ������� 6
% 1 - E�������� ����������, 3 - ������� ������, 5 - ���������� �����������.
% a - ����� �������� ������, d - ����������� �����, e - ����� ���������
% �����.

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
title('����������������� ������');
xlabel('X2');
ylabel('X1');

% ���������� ����� ���������

% ������� pdist � MATLAB ��������� ������ ���������� ����� ������� � ������ ������.
% ��� ��������� ���� ��������:
% X: �������, � ������� ������ ������� ������������ ����� ����� ������.
% ������� pdist ���������� ������ D, � ������� ������ ������� ������������ 
% ����� ���������� ����� ����� ������� ������. ������� ��������� � D 
% ������������� �������� ����� � X.
d_evkl = pdist(X, 'euclidean');
d_town = pdist(X, 'cityblock'); 
% d_town_matrix = squareform(d_town, 'tomatrix') % � ���� �������
d_mink = pdist(X, 'minkowski', 4); % 4 ��� �����������

% ���������� ��� ����������� �������� � �������� ������������� �������

% ������� linkage � MATLAB ������������ ��� �������� ������������� ���������� ��������� �� ������ ������.
% ��� ��������� ��� ���������:
% - X: �������, � ������� ������ ������� ������������ ����� ����� ������.
% - method: ����� �������������, ������� ����� ���� ����� �� ���������:
%   - 'single': ����� �������� ������.
%   - 'centroid': ���� ����� ���������� ��� ��������, ����������� �� ���������� 
% ����� ����������� ���� ���������. �������� - ��� ������� ����� ���� ����� � ��������.
%   - 'median': ���� ����� ���������� ��� ��������, ����������� �� ��������� ���������� ����� 
% ������� � ���� ���������. ������� - ��� ��������, ������� ����� ����� ������ 
% �� ��� �������� � ������ ������ �����.
% ������� linkage ���������� ������� ������ Z, ��� ������ ������ ������������ ����� ���������� �������.
link_evkl_sing = linkage(d_evkl, 'single');
link_evkl_centr = linkage(d_evkl, 'centroid');
link_evkl_media = linkage(d_evkl, 'median');
link_town_sing = linkage(d_town, 'single');
link_town_centr = linkage(d_town,'centroid');
link_town_media = linkage(d_town, 'median');
link_mink_sing = linkage(d_mink, 'single');
link_mink_centr = linkage(d_mink, 'centroid');
link_mink_media = linkage(d_mink, 'median');

% ������ �������� �������������:
disp('�������������� �������������� ������������');
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

disp('����� �����������: ������� ������ + ����������� �����');
max_effective = max(max(coph_corr_coeff))
disp('����� �������������: ���������� ����������� + ����� �������� ������');
min_effective = min(min(coph_corr_coeff))

figure;
dendrogram(link_town_centr);
title('M������ ������ + ����������� �����');

% ������������� ����� ���������:
N_clust = 3;

clust = cluster(link_mink_media, 'maxclust', N_clust); % construct clusters from linkages

figure;
gscatter(x, y, clust); % scatter plot by group
hold on;

% ������ ���������
clust_centr = zeros(N_clust, K); % ������� x � y ��������� �������
% ���������� ����� �������� ���������
dist_betw_clust = zeros(N_clust, N_clust);
% ������� ���������
clust_rad = zeros(N_clust, 1);
% ��������� ���������
clust_disp = zeros(N_clust, 1);

disp('���������� �� ��������� �������� �� ������');
for k = 1:N_clust
    disp(['�', num2str(k)]);
    % ������� ��� ������� ��������
    obj = find(clust == k); % ����� �� ����� �� clust ��������� � �������� ��������
    N_obj = length(obj); % �� ���-��
    singl_clust = zeros(N_obj, K); % ���������� ���� �������� ����������� �������� (�� ��������, 
    % ����� ���������� ���������� �������� ���������� �������� �������� �� �������������), ������� x � y
    
    % ������� ����� � �������� ��������
    for i = 1:N_obj
        singl_clust(i,:) = X(obj(i),:);
        clust_centr(k,:) = clust_centr(k,:) + singl_clust(i,:);
    end;
    clust_centr(k,:) = clust_centr(k,:)/N_obj; 
    
    % ���������� (��������������) �� ������ �� ���� ��������� ��������
    clust_dist_centr = zeros(N_obj, 1);
    for i = 1:N_obj
        for l = 1:K
            clust_dist_centr(i) = clust_dist_centr(i) + (singl_clust(i,l) - clust_centr(k,l))^2;
        end;
        % ��������� ��������
        clust_disp(k) = clust_disp(k) + clust_dist_centr(i);
        clust_dist_centr(i) = sqrt(clust_dist_centr(i));
    end;
    clust_disp(k) = clust_disp(k)/N_obj;
    clust_dist_centr
    % ������ �������� 
    clust_rad(k) = max(clust_dist_centr);
end;

% �������������� ���������� ����� ��������
for i = 1:N_clust
    for j = i+1:N_clust
        for l = 1:K
            dist_betw_clust(i,j) = dist_betw_clust(i,j) + (clust_centr(i,l) - clust_centr(j,l))^2; %
        end;
        dist_betw_clust(i,j) = sqrt(dist_betw_clust(i,j));
    end;
end;
disp('���������� ����� �������� ���������');
dist_betw_clust
disp('��������� ���������');
clust_disp
disp('������� ���������');
clust_rad

% ������ ���������
scatter(clust_centr(:,1), clust_centr(:,2), 20, 'k', 'filled');
t = 0:pi/180:2*pi;
for i = 1:N_clust
    x = clust_rad(i)*cos(t) + clust_centr(i,1);
    y = clust_rad(i)*sin(t) + clust_centr(i,2);
    plot(x, y, '-'); % ������ �������� �������
end;
title('��������� ��������');
xlabel('X2');
ylabel('X1');
