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

disp('1 - �-� ��������, 2 - �-� ��2, 3 - �-� ����');
disp('-------------------------------------------------------');

q1 = [1 1]; % ��������� �����������
res1 = nlinfit(x, y, @weibull, q1); % ���������� ������������� ������������� ��� ������� ��������
f1 = weibull(res1, x); % ��������� �������� ������� �������� ��� �������� �������� x
disp('������������������ ��������� �������� ��� �-��� ��������:');
disp(res1);

q2 = 2;
res2 = nlinfit(x, y, @Xi2, q2);
f2 = Xi2(res2, x);
disp('������������������ ��������� �������� ��� �-��� ��2:');
disp(res2);

q3 = [0.14 0.35];
res3 = nlinfit(x, y, @beta, q3);
f3 = beta(res3, x);
disp('������������������ ��������� �������� ��� �-��� ����:');
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

% u1, u2, u3 - ������� �������
% N - ����� ���������� - 1
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
title('��������� ������������� ����-��������')
xlabel('x')
ylabel('y')

% �������������������� ����������
sygm = array(:,3);

% ����� ��������� ���������� ���������� ��������
for i = 1:N
    x_kv_1 = x_kv_1 + ((f1(i) - y(i))/sygm(i))^2;
    x_kv_2 = x_kv_2 + ((f2(i) - y(i))/sygm(i))^2;
    x_kv_3 = x_kv_3 + ((f3(i) - y(i))/sygm(i))^2;
end;

% ����������
x_kv_1 = x_kv_1/u1;
x_kv_2 = x_kv_2/u2;
x_kv_3 = x_kv_3/u3;

disp('-------------------------------------------------------');
%disp(newline);

disp('������������� �������� Xi2 ��� ������� ��������:');
disp(x_kv_1);

disp('������������� �������� Xi2 ��� ������� �� �������:');
disp(x_kv_2);

disp('������������� �������� Xi2 ��� ������� ����:');
disp(x_kv_3);

% ���������� �������
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
title('���������� ������� 1');
figure;
plot(R2);
title('���������� ������� 2');
figure;
plot(R3);
title('���������� ������� 3');
histogram(R3);

[res3, r, J, covb] = nlinfit(x, y, @beta, q3);

% ������������������ ������� ���������� ��������
A1 = zeros(N/2);
A2 = zeros(N/2);
A3 = zeros(N/2);

% ���������� ������������������ �������:
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
title('������������������ ������� 1');
figure;
plot(k, A2);
title('������������������ ������� 2');
figure;
plot(k, A3);
title('������������������ ������� 3');

% ������������� �������� - 1
[res3, r, J, covb] = nlinfit(x, y, @beta, q3);

disp('-------------------------------------------------------');

% ������������� �������� �� t-����������

alpha = 0.32; % ������� ����������

% �������� ������������
if ~isnumeric(x(2)) || ~isnumeric(q3(2))
    error('x(2) or q3(2) is not a number!');
end

df = max(length(x) - length(q3) - 1, 1); % ������� �������

t = tinv(1-alpha/2, df); % ����������� �������� t-����������

se = sqrt(diag(covb)); % ����������� ������ ������

lower_bound = res3 - t * se'; % ������ �������
upper_bound = res3 + t * se'; % ������� �������

% ����� ����������
disp('������������� ������������� ��������:');
disp([lower_bound, upper_bound]);





