clear
clc
population_size = 100;      % 种群大小
chromosome_size = 40;       % 染色体长度
generation_size = 400;      % 最大迭代次数
cross_rate = 0.5;           % 交叉概率
mutate_rate = 0.1;         % 变异概率
serus=7;                    % serus的个数
% pitch_kinds=5;              % 批次种类
% workes=10;                  % 工人数目
% pitches=30;                 % 批次总数
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
population=zeros(population_size,chromosome_size); %初始化数组
population_new=zeros(population_size,chromosome_size);
fitness_sum=zeros(1,population_size);
fitness_value=zeros(1,population_size);
fitness_average=zeros(1,generation_size);
best_individual=zeros(1,chromosome_size);
best_fitness = 0.;
best_generation = 0;
elitism=0;       %是否进行精英选择；
for i=1:population_size
    for j=1:chromosome_size
        population(i,j) =round(1+(serus-1)*rand);  % rand产生(0,1)之间的随机数，round()是四舍五入函数
    end
end

for G=1:generation_size 
    %种群可行解修正
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
%计算每个个体的适应度
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
%对个体按适应度大小进行排序，并且保存最佳个体
min_index = 1;
temp = 1;
temp_chromosome=zeros(1,chromosome_size);
% 遍历种群 
% 冒泡排序
% 最后population(i)的适应度随i递增而递增，population(1)最小，population(population_size)最大
for i=1:population_size
    min_index = i;
    for j = i+1:population_size
        if fitness_value(j) < fitness_value(min_index)
            min_index = j;
        end
    end
    
    if min_index ~= i
        % 交换 fitness_value(i) 和 fitness_value(min_index) 的值
        temp = fitness_value(i);
        fitness_value(i) = fitness_value(min_index);
        fitness_value(min_index) = temp;
        % 此时 fitness_value(i) 的适应度在[i,population_size]上最小
        
        % 交换 population(i) 和 population(min_index) 的染色体串
        for k = 1:chromosome_size
            temp_chromosome(k) = population(i,k);
            population(i,k) = population(min_index,k);
            population(min_index,k) = temp_chromosome(k);
        end
    end
end
 
% fitness_sum(i) = 前i个个体的适应度之和
for i=1:population_size
    if i==1
        fitness_sum(i) = fitness_sum(i) + fitness_value(i);    
    else
        fitness_sum(i) = fitness_sum(i-1) + fitness_value(i);
    end
end
 
% fitness_average(G) = 第G次迭代 个体的平均适应度
% fitness_average(G) = fitness_sum(population_size)/population_size;

% 更新最大适应度和对应的迭代次数，保存最佳个体(最佳个体的适应度最大)
if fitness_value(population_size) > best_fitness
    best_fitness = fitness_value(population_size);
    best_generation = G;
    for j=1:chromosome_size
        best_individual(j) = population(population_size,j);
    end
end
% 轮盘赌选择操作
for i=1:population_size
    r = rand * fitness_sum(population_size);  % 生成一个随机数，在[0,总适应度]之间
    first = 1;
    last = population_size;
    mid = round((last+first)/2);
    idx = -1;
    
    % 排中法选择个体
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
   
   % 产生新一代个体
   for j=1:chromosome_size
        population_new(i,j) = population(idx,j);
   end
end
 % 是否精英选择
if elitism
    p = population_size-1;
else
    p = population_size;
end

for i=1:p
   for j=1:chromosome_size
       % 如果精英选择，将population中前population_size-1个个体更新，最后一个最优个体保留
       population(i,j) = population_new(i,j);
   end
end

% 步长为2 遍历种群
for i=1:2:population_size
    % rand<交叉概率，对两个个体的染色体串进行交叉操作
    if(rand < cross_rate)
        cross_position = round(rand * chromosome_size);
        if (cross_position == 0 || cross_position == 1)
            continue;
        end
        % 对 cross_position及之后的进行交换
        for j=cross_position:chromosome_size
            temp = population(i,j);
            population(i,j) = population(i+1,j);
            population(i+1,j) = temp;
        end
    end
end
for i=1:population_size
    if rand < mutate_rate
        mutate_position = round(rand*chromosome_size);  % 变异位置
        if mutate_position == 0
            % 若变异位置为0，不变异
            continue;
        end
        population(i,mutate_position) = round(1+(serus-1)*rand);     
    end
end % 变异概率
end

%打印算法迭代过程
x = 1:1:generation_size;
y = fitness_average;
%plot(x,y)
m = best_individual;    % 获得最佳个体
n = best_fitness;       % 获得最佳适应度
p = best_generation;    % 获得最佳个体出现时的迭代次数
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
    fprintf('Seru%d工时为：\n',k);
    disp(seru_sum);
    fprintf('Seru%d中工人有：\n',k);
    disp(b);
    fprintf('Seru%d中分配的批次有：\n',k);
    disp(c);
end
disp('Seru均值为：');
disp(mean(serus_value));
disp('Seru方差为：');
disp(var(serus_value));

