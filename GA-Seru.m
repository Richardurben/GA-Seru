clear
clc
population_size = 100;      % ��Ⱥ��С
chromosome_size = 40;       % Ⱦɫ�峤��
generation_size = 400;      % ����������
cross_rate = 0.5;           % �������
mutate_rate = 0.1;         % �������
serus=7;                    % serus�ĸ���
% pitch_kinds=5;              % ��������
% workes=10;                  % ������Ŀ
% pitches=30;                 % ��������
workers_index=[0.92 0.96 1.04 1.09 1.20;
               0.95 0.97 1.09 1.12 1.18;
               0.99 1.01 1.05 1.09 1.21;
               1.03 1.07 1.09 1.12 1.25;
               0.96 1.02 1.05 1.10 1.18;
               1.01 1.10 1.10 1.15 1.23;
               1.04 1.07 1.09 1.17 1.24;
               0.98 1.02 1.10 1.11 1.20;
               0.97 1.03 1.12 1.19 1.26;
               0.98 1.06 1.13 1.18 1.28];
% pitches_index=[0,0,0,0,49,0,54,0,0,0,0,0,0,0,0,0,54,0,0,0,53,0,0,0,0,0,53,0,0,0;
%                0,0,0,0,0,0,0,48,48,0,46,0,0,0,0,0,0,0,54,0,0,0,0,0,45,0,0,0,53,0;
%                55,0,54,0,0,0,0,0,0,48,0,0,48,0,0,0,0,0,0,0,0,46,0,0,0,44,0,0,0,0;
%                0,0,0,49,0,55,0,0,0,0,0,58,0,52,0,0,0,57,0,0,0,0,45,0,0,0,0,47,0,0;
%                0,53,0,0,0,0,0,0,0,0,0,0,0,0,48,51,0,0,0,49,0,0,0,46,0,0,0,0,0,0];
pitches_index=[[3,55];[5,53];[3,54];[4,49];[1,49];[4,55];[1,54];[2,48];[2,48];[3,48];
               [2,46];[4,58];[3,48];[4,52];[5,48];[5,51];[1,54];[4,57];[2,54];[5,49];
               [1,53];[3,46];[4,45];[5,46];[2,45];[3,44];[1,53];[4,47];[2,53];[3,52]];
population=zeros(population_size,chromosome_size); %��ʼ������
population_new=zeros(population_size,chromosome_size);
fitness_sum=zeros(1,population_size);
fitness_value=zeros(1,population_size);
fitness_average=zeros(1,generation_size);
best_individual=zeros(1,chromosome_size);
best_fitness = 0.;
best_generation = 0;
elitism=0;       %�Ƿ���о�Ӣѡ��
for i=1:population_size
    for j=1:chromosome_size
        population(i,j) =round(1+(serus-1)*rand);  % rand����(0,1)֮����������round()���������뺯��
    end
end

for G=1:generation_size 
    %��Ⱥ���н�����
    for i=1:population_size
    a=population(i,1:10);
    [c,xout] = hist(a,[1:serus]);
    z=ismember(0,c);
        while z
            k=find(c==max(c));
            d=find(a==k(1));
            f=find(c==0);
            a(d(1))=f(1);            
            [c,xout] = hist(a,[1:serus]);
            z=ismember(0,c);
            population(i,1:10)=a;
        end
        a=[];k=[];d=[];f=[];c=[];
    end
for i=1:population_size
    a=population(i,11:40);
    [c,xout] = hist(a,[1:serus]);
    z=ismember(0,c);
        while z
            k=find(c==max(c));
            d=find(a==k(1));
            f=find(c==0);
            a(d(1))=f(1);
            [c,xout] = hist(a,[1:serus]);
            z=ismember(0,c);
            population(i,11:40)=a;   
        end
        a=[];k=[];d=[];f=[];c=[];
end 
%����ÿ���������Ӧ��
for i=1:population_size
    serus_value=zeros(1,serus);
    a=population(i,:);
    for k=1:serus
    b=find(a(1:10)==k);
    c=find(a(11:40)==k);
    d=pitches_index(c,1).';
    e=pitches_index(c,2).';
    for j=1:length(c)
        seru_sum=0;
        seru_sum=seru_sum+sum(workers_index(b,d(j)))*5/length(b)/length(b)*e(j);
    end
    serus_value(1,k)=seru_sum;
    end
 fitness_value(i)=1/(mean(serus_value)+var(serus_value)/100);
end
%�Ը��尴��Ӧ�ȴ�С�������򣬲��ұ�����Ѹ���
min_index = 1;
temp = 1;
temp_chromosome=zeros(1,chromosome_size);
% ������Ⱥ 
% ð������
% ���population(i)����Ӧ����i������������population(1)��С��population(population_size)���
for i=1:population_size
    min_index = i;
    for j = i+1:population_size
        if fitness_value(j) < fitness_value(min_index)
            min_index = j;
        end
    end
    
    if min_index ~= i
        % ���� fitness_value(i) �� fitness_value(min_index) ��ֵ
        temp = fitness_value(i);
        fitness_value(i) = fitness_value(min_index);
        fitness_value(min_index) = temp;
        % ��ʱ fitness_value(i) ����Ӧ����[i,population_size]����С
        
        % ���� population(i) �� population(min_index) ��Ⱦɫ�崮
        for k = 1:chromosome_size
            temp_chromosome(k) = population(i,k);
            population(i,k) = population(min_index,k);
            population(min_index,k) = temp_chromosome(k);
        end
    end
end
 
% fitness_sum(i) = ǰi���������Ӧ��֮��
for i=1:population_size
    if i==1
        fitness_sum(i) = fitness_sum(i) + fitness_value(i);    
    else
        fitness_sum(i) = fitness_sum(i-1) + fitness_value(i);
    end
end
 
% fitness_average(G) = ��G�ε��� �����ƽ����Ӧ��
% fitness_average(G) = fitness_sum(population_size)/population_size;

% ���������Ӧ�ȺͶ�Ӧ�ĵ���������������Ѹ���(��Ѹ������Ӧ�����)
if fitness_value(population_size) > best_fitness
    best_fitness = fitness_value(population_size);
    best_generation = G;
    for j=1:chromosome_size
        best_individual(j) = population(population_size,j);
    end
end
% ���̶�ѡ�����
for i=1:population_size
    r = rand * fitness_sum(population_size);  % ����һ�����������[0,����Ӧ��]֮��
    first = 1;
    last = population_size;
    mid = round((last+first)/2);
    idx = -1;
    
    % ���з�ѡ�����
    while (first <= last) && (idx == -1) 
        if r > fitness_sum(mid)
            first = mid;
        elseif r < fitness_sum(mid)
            last = mid;     
        else
            idx = mid;
            break;
        end
        mid = round((last+first)/2);
        if (last - first) == 1
            idx = last;
            break;
        end
    end
   
   % ������һ������
   for j=1:chromosome_size
        population_new(i,j) = population(idx,j);
   end
end
 % �Ƿ�Ӣѡ��
if elitism
    p = population_size-1;
else
    p = population_size;
end

for i=1:p
   for j=1:chromosome_size
       % �����Ӣѡ�񣬽�population��ǰpopulation_size-1��������£����һ�����Ÿ��屣��
       population(i,j) = population_new(i,j);
   end
end

% ����Ϊ2 ������Ⱥ
for i=1:2:population_size
    % rand<������ʣ������������Ⱦɫ�崮���н������
    if(rand < cross_rate)
        cross_position = round(rand * chromosome_size);
        if (cross_position == 0 || cross_position == 1)
            continue;
        end
        % �� cross_position��֮��Ľ��н���
        for j=cross_position:chromosome_size
            temp = population(i,j);
            population(i,j) = population(i+1,j);
            population(i+1,j) = temp;
        end
    end
end
for i=1:population_size
    if rand < mutate_rate
        mutate_position = round(rand*chromosome_size);  % ����λ��
        if mutate_position == 0
            % ������λ��Ϊ0��������
            continue;
        end
        population(i,mutate_position) = round(1+(serus-1)*rand);     
    end
end % �������
end

%��ӡ�㷨��������
x = 1:1:generation_size;
y = fitness_average;
%plot(x,y)
m = best_individual;    % �����Ѹ���
n = best_fitness;       % ��������Ӧ��
p = best_generation;    % �����Ѹ������ʱ�ĵ�������
serus_value=zeros(1,serus);
a=best_individual;
for k=1:serus
    b=find(a(1:10)==k);
    c=find(a(11:40)==k);
    d=pitches_index(c,1).';
    e=pitches_index(c,2).';
    seru_sum=0;
    for j=1:length(c)
        seru_sum=seru_sum+sum(workers_index(b,d(j)))*5/length(b)/length(b)*e(j);
    end
    serus_value(1,k)=seru_sum;
    fprintf('Seru%d��ʱΪ��\n',k);
    disp(seru_sum);
    fprintf('Seru%d�й����У�\n',k);
    disp(b);
    fprintf('Seru%d�з���������У�\n',k);
    disp(c);
end
disp('Seru��ֵΪ��');
disp(mean(serus_value));
disp('Seru����Ϊ��');
disp(var(serus_value));

