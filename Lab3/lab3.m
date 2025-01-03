clc
close all
clear all

% mean(x) - ��������� ������� �������� �������.
% var(x) - ��������� ��������� �������.
% chi2inv(p,df) - ������� �������� ������� ������������� ��-������� 
%   p - �����������, df - ����� �������� �������.
% [V,D] = eig(A) - ������� ����������� �������� � ������� �������.
%   A: ���������� �������.
%   V: �������, ������� ������� �������� ������������ ��������� A.
%   D: ������������ �������, ������������ �������� ������� �������� ������������ ���������� A.
% y = flipud(x) - �������� ������ �� ���������.
% y = fliplr(x) - �������� ������ �� �����������.
% scatter(x,y) - ������� ��������� ���������.

file = fopen('data6.txt');
X = fscanf(file, '%d', [8 inf])' % ������� � �� 8 ���������
fclose(file);

X_mat = mean(X) % ������� �������� ������� �� ��������
srednekv_otkl = std(X); % ������������������ ����������
N = size(X);
K = N(2);
N = N(1);

% ���������� ������� ������ �� �������� (��� ���������� ��������������
% �����������)
% ����������� ���������� ������� �� ���� ��������
X_norm = zeros(N, K);
for i = 1:N
    for j = 1:K
        X_norm(i,j) = (X(i,j) - X_mat(j))/srednekv_otkl(j);
    end;
end;

X_norm_ = mean(X_norm)
srednekv_otkl = std(X_norm)

% ������� ���������� (���) ��������� � �������������� ������ ������������
R = (X_norm'*X_norm)/(N-1)

% �������� ������� ���������� (������ ������� ���������� ��� ���������) ��
% ������� �� ���������
d = 0;
for i = 1:K
    for j = (i+1):K
        d = d + (R(i,j))^2;
    end;
end;
d = d*N % > x_kv 
x_kv = chi2inv(0.95, K*(K-1)/2)

[A, L] = eig(R); 
A = fliplr(A);
% A - ���. �������, L (lambda) - ������������ ������� � ���. ���������� 

% ����������� �� �������� (�������� �� ����������� � ���������)
L = flipud(L); 
L = fliplr(L);
Z = X_norm*A % �������� �������� �� ��. ����������

sum_comp = sum(var(Z)) % ����� ���������� ��������� �������� �������� �� ������� ����������
sum_prizn = sum(var(X_norm)) % ����� ���������� ��������� �������� ���������

otn_dol_razbr_j = var(Z) / sum_comp; % (3.7) - ���. ���� ��������, ������������ �� j-� ������� ����������

otn_dol_razbr_i = zeros(K,1); % (3.8) - ���. ���� ��������, ������������ �� i ������ ���������
   for i = 1:K
    otn_dol_razbr_i(i) = 0;
      for j = 1:i
        otn_dol_razbr_i(i) = otn_dol_razbr_i(i) + otn_dol_razbr_j(j);
      end
   end
otn_dol_razbr_i = otn_dol_razbr_i';

%otn_dol_razbr = otn_dol_razbr_j
otn_dol_razbr = otn_dol_razbr_i

% ������������� ���� �������� ��� ������ ���� ���������
otn_dol_razbr_1_2 = otn_dol_razbr(1) + otn_dol_razbr(2); 

cov_matrix = cov(Z) % ������� ���������� ��� �������� �������� �� ������� ����������

% ����������� ������ ���� 8, ����� 2
% ������ ��� ���������� ������� �� �������� ����������� ������������ ��������, 
% �� ���� ���������� ���������
figure;
scatter(Z(:,1),Z(:,2), 35, '.', 'r');
title('��������� ��������� (��� ������ ���� ���������)');
xlabel('Z1');
ylabel('Z2');
