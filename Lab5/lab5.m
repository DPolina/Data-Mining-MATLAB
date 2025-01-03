clc
close all
clear all

file = fopen('data6.txt');
X = fscanf(file, '%f', [2 inf])';
fclose(file);

x = X(:,1);
y = X(:,2);
N = size(X);
K = N(2); % K = 2
N = N(1); % N = 100

figure;
scatter(x, y, '.');
title('����������������� ������');
xlabel('X2');
ylabel('X1');

k = 4; % ������� ���������
C = zeros(k, K); % ������ ���������
for i = 1:k
    C(i,:) = X(i,:); % ������ k �������� �������� ������
end;

U = zeros(N, 2); % U [������ ��������, ���������� �� ������� �� ������]
Q_prev = 1435; % �������� ����������� �������� ������������� �� ���������� ����
eps = 0.0001435; % ��������
m = 1; % ����� ��������
clust_rad = zeros(1, k); % ������ �������� ������� �� ���������

while (true)
    R = zeros(1, k);
    % ������ ���������� �� ������� ���������, ���������� ������� U
    for i = 1:N 
        for n = 1:k
        R(n) = pdist([X(i,:); C(n,:)], 'minkowski', 4); % ���������� i-�� ������� �� ������� �� �������
        end 
        
        [d, n] = min(R); 
        U(i,1) = n; % ������� � ��������, ���������� �� �������� ����������
        U(i,2) = R(n); % ���������� �� ������ ��������
    end
    
    Q_m = 0; % �������� ������������� �� ������� ���� ��� p-�� ��������
    QQ = 0; % �������� ������������� �� ������� ����
    for p = 1:k
        Obj = find(U(:,1) == p); % ����� �� �������� ��������� � p-�� �������� (�� ���������� �����)
        s = length(Obj); % ������� ����� ���������
        singl_clust = zeros(1, s);
        X_clust = zeros(s, 2);
        for i = 1:s % ���� �� ���� ��������� p-�� ��������
            singl_clust(i) = U(Obj(i), 2); % ���������� �� ������ i-�� ��������
            X_clust(i,:) = X(Obj(i),:); % ���������� � ����� ������� ������ �� ������� ������� � p-�� �������� 
        end;
        xex = pdist(X_clust); % Euclidean �� ���������, ���������� ����� ������� ��������
        Q_m = sum(xex); % ����� ����������
        clust_rad(p) = max(singl_clust); % ������ �������� ��� ������������ ���������� �� ������ �� ������-���� ��� �������
        QQ = QQ + Q_m;
    end;
    
    % ���� �� ���� �������� ����� �������� ������������� ����� �� ����������, �� ���������� �������� ���������������
    if abs(QQ - Q_prev) <= eps
        break;
    else % ���� �������� �������� �� ����������, �� ������������� ������ ��� ������� �� ������� ��������
        for l = 1:k 
             Obj = find(U(:,1)==l); % ����� �� �������� ��������� � l-�� ��������
             s = length(Obj); % ������� ����� ���������
             for j = 1:K % �� ������� ��������, �� ���� ������� �� ����� � ������� �� �������
                summa = 0;
                for i = 1:s 
                    summa = summa + X(Obj(i), j);
                end;
                C(l,j) = summa/s; % �������
             end;
        end;
        Q_prev = QQ; % �������� �������� �������������
        m = m + 1;
    end;
end;

% ����������� ��������� ��������� � �������
figure
gscatter(x, y, U(:,1), '', '.', 15);
hold on;
t = 0:(pi/180):2*pi;
for i = 1:k
  x = clust_rad(i)*cos(t) + C(i,1);
  y = clust_rad(i)*sin(t) + C(i,2);
  plot(x, y, 'k');
end;
scatter(C(:,1), C(:,2), 20, 'k', 'filled'); 
title('��������� �������� � �� ������');
xlabel('X2');
ylabel('X1');
