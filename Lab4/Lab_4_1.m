clear all
close all
clc

% ������� 6
% 1 - E�������� ����������, 3 - ������� ������, 5 - ���������� �����������.
% a - ����� �������� ������, d - ����������� �����, e - ����� ���������
% �����.

fid = fopen('data1.txt', 'r');
X0 = fscanf(fid, '%g %g', [2 inf]); % 24x2 
fclose(fid);

X=X0';
figure
plot(X(:,1),X(:,2), '.');

K = zeros(3); % �������������� ������������

% ��������� �������
DEuclidean = pdist(X, 'euclidean');

EuclideanSingle = linkage(DEuclidean,'single'); % �������� "�������� ������".
EuclideanComplete = linkage(DEuclidean,'complete'); % �������� "�������� ������".
EuclideanAverage = linkage(DEuclidean,'average'); % �������� "������� �����". ������

K(1,1) = cophenet(EuclideanSingle, DEuclidean);
K(2,1) = cophenet(EuclideanComplete, DEuclidean);
K(3,1) = cophenet(EuclideanAverage, DEuclidean);

% ������������������� ��������� �������
DSeuclidean = pdist(X, 'seuclidean');

DSeuclideanSingle = linkage(DSeuclidean,'single'); %�������� "�������� ������". ������
DSeuclideanComplete = linkage(DSeuclidean,'complete'); %�������� "�������� ������".
DSeuclideanAverage = linkage(DSeuclidean,'average'); %�������� "������� �����".

K(1,2) = cophenet(DSeuclideanSingle, DSeuclidean);
K(2,2) = cophenet(DSeuclideanComplete, DSeuclidean);
K(3,2) = cophenet(DSeuclideanAverage, DSeuclidean);

% ������� ������
DCityblock = pdist(X, 'cityblock');

DCityblockSingle = linkage(DCityblock,'single'); %�������� "�������� ������".
DCityblockComplete = linkage(DCityblock,'complete'); %�������� "�������� ������".
DCityblockAverage = linkage(DCityblock,'average'); %�������� "������� �����".

K(1,3) = cophenet(DCityblockSingle, DCityblock);
K(2,3) = cophenet(DCityblockComplete, DCityblock);
K(3,3) = cophenet(DCityblockAverage, DCityblock);

figure
dendrogram(EuclideanAverage);
title('EuclideanAverage')

T = cluster(EuclideanAverage, 'maxclust', 4);

figure
gscatter(X(:,1),X(:,2),T);
title('����������� ��������')
hold on

%���������� �������� �� ��������
k1=1; k2=1; k3=1; k4=1;
for i=1:24
     if (T(i) == 1)
        clust1(k1,:) = X(i,:);
        k1 = k1+1;
     end
     if (T(i) == 2)
        clust2(k2,:) = X(i,:);
        k2 = k2+1;
     end
     if (T(i) == 3)
        clust3(k3,:) = X(i,:);
        k3 = k3+1;
     end
     if (T(i) == 4)
        clust4(k4,:) = X(i,:);
        k4 = k4+1;
     end
end

M = [mean(clust1);mean(clust2);mean(clust3);mean(clust4)];
VAR = [var(clust1);var(clust2);var(clust3);var(clust4)];

DCenters = pdist(M,'euclidean'); % ���������� ����� �������� ���������

Dclust1 = squareform(pdist([clust1;mean(clust1)],'euclidean'));
Dclust1 = Dclust1(1,:);

Dclust2 = squareform(pdist([clust2;mean(clust2)],'euclidean'));
Dclust2 = Dclust2(1,:);

Dclust3 = squareform(pdist([clust3;mean(clust3)],'euclidean'));
Dclust3 = Dclust3(1,:);

Dclust4 = squareform(pdist([clust4;mean(clust4)],'euclidean'));
Dclust4 = Dclust4(1,:);
plot(M(:,1),M(:,2), 'p');
hold off