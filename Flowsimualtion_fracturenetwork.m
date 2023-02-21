clear
clc
tic
Value_left = 20;
Value_right = 0;
original = textread('case6_result.txt');
n = size(original,1);         %裂隙个数
% n = 2;
% a = geometric_condition(3) - geometric_condition(1)/4;   %左边界
% b = geometric_condition(3) + geometric_condition(1)/4;   %右边界
% c = geometric_condition(4) - geometric_condition(2)/4;   %下边界
% d = geometric_condition(4) + geometric_condition(2)/4;   %上边界
a = 10; b = 70; c = 10; d = 50;
%% Extraction starting and ending point
x = zeros(n,1); x1 = zeros(n,1);  y = zeros(n,1); y1 = zeros(n,1); 
for i = 1:n
    if cos(original(i,3) * pi / 180)>0
        x(i) = original(i,1) - original(i,4)*cos(original(i,3) * pi / 180)/2;
        x1(i) = original(i,1) + original(i,4)*cos(original(i,3) * pi / 180)/2; 
        y(i) = original(i,2) - original(i,4)*sin(original(i,3) * pi / 180)/2;  
        y1(i) = original(i,2) + original(i,4)*sin(original(i,3) * pi / 180)/2;
    elseif cos(original(i,3) * pi / 180)<0
        x(i) = original(i,1) + original(i,4)*cos(original(i,3) * pi / 180)/2;
        x1(i) = original(i,1) - original(i,4)*cos(original(i,3) * pi / 180)/2; 
        y(i) = original(i,2) + original(i,4)*sin(original(i,3) * pi / 180)/2;  
        y1(i) = original(i,2) - original(i,4)*sin(original(i,3) * pi / 180)/2;
    else
        x(i) = original(i,1);
        x1(i) = original(i,1);
        y(i) = original(i,2) + l/2*original(i,4);
        y1(i) = original(i,2) - l/2*original(i,4);
        printf('存在90度裂隙')
    end
end
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)');
% for i = 1 : n
%   plot([x(i),x1(i)], [y(i),y1(i)]);
%      hold on
% end
% plot([0,0],[0,50],'-- r'); plot([50,50],[10,50],'-- r');
% plot([0,50],[0,0],'-- r'); plot([0,50],[50,50],'-- r');
% plot([a,a],[c,d],'-- r'); plot([b,b],[c,d],'-- r'); 
% plot([a,b],[c,c],'-- r'); plot([a,b],[d,d],'-- r');


% Delete the fracture outside the simualtion domain
 w = ( x(:) >= b | x1(:) <= a ) | (( y(:) >= d & y1(:) >= d ) | ( y(:) <= c & y1(:) <= c ));  
 x = x(w ~= 1); x1 = x1(w ~= 1); y = y(w ~= 1); y1 = y1(w ~= 1);
 original = original(w ~= 1,:);
 n = size(x,1);
 disp('Extraction starting and ending point is end');
 %%
%  figure;
%  axis equal;
%  axis([a b c d]);
%  xlabel('(m)');
%  ylabel('(m)');
%  for i = 1 : n
%      plot([x(i),x1(i)], [y(i),y1(i)]);
%      hold on
%  end
%  plot([a,a],[c,d],'-- r'); plot([b,b],[c,d],'-- r');
%  plot([a,b],[c,c],'-- r'); plot([a,b],[d,d],'-- r');
%%
for i = 1 : n
    if (x(i)<b) && (x1(i)>b) && ((y(i)<c) && (y1(i)>c) || (y(i)>d) && (y1(i)<d))
        y1(i) = y1(i) + tan(original(i,3))*(b - x1(i)); x1(i) = b;
        if (y1(i) <=  c) || (y1(i) >=  d)
            
            x(i) = inf; x1(i) = inf; y(i) = inf; y1(i) = inf; original(i,:) = inf;
        end
    end
    if (x(i)<a) && (x1(i)>a) && ((y(i)<d) && (y1(i)>d) || (y(i)>c) && (y1(i)<c))
        y(i) = y(i) + tan(original(i,3))*(a - x(i)); x(i) = a;
        if (y(i) <=  c) || (y(i) >=  d)
            x(i) = inf; x1(i) = inf; y(i) = inf; y1(i) = inf; original(i,:) = inf;
        end
    end
end
x = x(x~=inf); x1 = x1(x1~=inf);
y = y(y~=inf); y1 = y1(y1~=inf);
original(all(original==inf,2),:) = [];
n = size(x,1);
neighbor = ones(n,n);
 
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)');
% for i = 1 : n
%          plot([x(i),x1(i)], [y(i),y1(i)]);
%      hold on
% end
% plot([a,a],[c,d],'-- r'); plot([b,b],[c,d],'-- r'); 
% plot([a,b],[c,c],'-- r'); plot([a,b],[d,d],'-- r');
% Reset the endpoint of fracture
 for i = 1 : n
    if x(i) < a 
        x(i) = a;
        y(i) = original(i,2) + tan(original(i,3) * pi / 180)*(a - original(i,1));
    end
    if x1(i) < a 
        x1(i) = a;
        y1(i) = original(i,2) + tan(original(i,3) * pi / 180)*(a - original(i,1));
    end
    if x1(i) > b
        x1(i) = b;
        y1(i) = original(i,2) + tan(original(i,3) * pi / 180)*(b - original(i,1));
    end
    if x(i) > b
        x(i) = b;
        y(i) = original(i,2) + tan(original(i,3) * pi / 180)*(b - original(i,1));
    end
    if y(i) < y1(i)
        
        if y(i) < c 
            y(i) = c;
            x(i) = original(i,1) + (c - original(i,2))/tan(original(i,3) * pi / 180);
        end
        
        if y1(i) > d 
            y1(i) = d;
            x1(i) = original(i,1) + (d - original(i,2))/tan(original(i,3) * pi / 180);
        end
    elseif y(i) > y1(i)
        
       if y1(i) < c 
            y1(i) = c;
            x1(i) = original(i,1) + (c - original(i,2))/tan(original(i,3) * pi / 180);
       end
       
        if y(i) > d
            y(i) = d;
            x(i) = original(i,1) + (d - original(i,2))/tan(original(i,3) * pi / 180);
        end 
        
    end
 end
% %导出主干裂隙通道 
% fid_dat = fopen('zone_case6_result.dat','a');
% fprintf(fid_dat,'Title ="fracfure"\nVARIABLES = "X""Y"\n');
%  for i = 1:size(x,1)
%     fprintf(fid_dat,'ZONE T="ZONE %g"\n',i);
%     fprintf(fid_dat,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
%     fprintf(fid_dat,'%g\t%g\n%g\t%g\n',x(i),y(i),x1(i),y1(i));
%  end
%  fclose(fid_dat);

all_x = x;
all_y = y;
all_x1 = x1;
all_y1 = y1;

figure;
axis equal;
axis([a b c d]);
xlabel('(m)');
ylabel('(m)');
for i = 1 : n
    plot([x(i),x1(i)], [y(i),y1(i)]);
    hold on

end
daspect([1 1 1]);
% plot([40,40],[40,70],'-- r'); plot([70,70],[40,70],'-- r'); 
% plot([40,70],[40,40],'-- r'); plot([40,70],[70,70],'-- r');
% plot([a,a],[c,d],'-- r'); plot([b,b],[c,d],'-- r'); 
% plot([a,b],[c,c],'-- r'); plot([a,b],[d,d],'-- r');
%% 快速排斥实验
for i = 1:n
    neighbor(i,i) = 0;
    for j = i+1:n
        if x1(i)<x(j)
            neighbor(i,j) = 0;
        elseif max(y(i),y1(i))<min(y(j),y1(j))
            neighbor(i,j) = 0;      
        end
        neighbor(j,i) = neighbor(i,j);
    end
end
disp('快速排斥实验结束');
%% 跨立实验
for i = 1:n
    for j = i+1:n
        if neighbor(i,j) == 1
            CA = [x(i)-x(j),y(i)-y(j)];
            CD = [x1(j)-x(j),y1(j)-y(j)];
            CB = [x1(i)-x(j),y1(i)-y(j)];
            AC = [x(j)-x(i),y(j)-y(i)];
            AB = [x1(i)-x(i),y1(i)-y(i)];
            AD = [x1(j)-x(i),y1(j)-y(i)];
            if ((CA(1)*CD(2)-CA(2)*CD(1))*...
                    (CB(1)*CD(2)-CB(2)*CD(1))<=0)&&...
                    ((AC(1)*AB(2)-AC(2)*AB(1))*(AD(1)*AB(2)-...
                    AD(2)*AB(1))<=0)
                neighbor(i,j) = 1;
            else
                neighbor(i,j) = 0;
            end
            neighbor(j,i) = neighbor(i,j);
        end
    end
end
neighbor_oo = neighbor;
disp('跨立试验结束');
%% 补充对角线
% neighbor = neighbor + neighbor';
%% 找出孤立裂隙
alone = [];
g = 1;
for i = 1:n
    if sum(neighbor(i,:)) == 0
        alone(g,1) = i;
        g = g + 1;
    end
end

%% Delete the alone fractures
m = size(alone,1);
for i = 1:m
    j = alone(i);
    x(j) = inf; x1(j) = inf; y(j) = inf; y1(j) = inf;
    neighbor(j,:) = inf; neighbor(:,j) = inf;
    original(j,:) = inf; 
end
        
x = x(x~=inf); x1 = x1(x1~=inf);
y = y(y~=inf); y1 = y1(y1~=inf); 
neighbor(:,all(neighbor==inf,1)) = [];
neighbor(all(neighbor==inf,2),:) = [];
original(all(original==inf,2),:) = [];
Connected_fracture = [x y x1 y1 original(:,5)];
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)')
% n = length(x);
% for i = 1:n
%     hold on
%     plot([x(i),x1(i)],[y(i),y1(i)],'k');
% end
%% 确定东西方向起点终点
t = 1;
s = 1;
for i = 1:size(x,1)
    if x(i) == a
        start_and_end(t,1) = i;
        t = t+1;
    end
    
    if x1(i) == b    
        start_and_end(s,2) = i;
        s = s+1; 
    end
end
t = t-1;
s = s-1;
%% 绘制主干网
f = 1;
totalorder = [];
for i = 1:t
    Circle_num = i
    for j = 1:s
        result = findPaths(neighbor,start_and_end(i,1),start_and_end(j,2));
        [~,lengthresult] = size(result);
        for ss = 1:lengthresult
            q = result{1,ss};
            for sss = 1:size(q,2)
                totalorder(f) = q(1,sss);
                f = f+1;
            end
        end
    end
end

totalorder = unique(totalorder);


% disp('绘制主干裂隙网络结束');
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)');
% for j = 1:length(frac_eliminated)
%     i = frac_eliminated(j);
%     hold on;
%     plot([all_x(i),all_x1(i)],[all_y(i),all_y1(i)],'k');
% end
% 剔除非联通裂隙
for i = 1:size(x,1)
    k = 1;
    for j = 1:size(totalorder,2)
        if i == totalorder(j)
            k = 0;
            break
        end
    end
    if k == 1
        x(i) = inf; y(i) = inf; x1(i) = inf; y1(i) = inf;
        original(i,:) = inf;
    end
end
x = x(x~=inf); x1 = x1(x1~=inf);
y = y(y~=inf); y1 = y1(y1~=inf); 
original(all(original==inf,2),:) = [];
Connected_fracture = [x y x1 y1 original(:,5)];
 

% save('Connected fracture1.txt','Connected_fracture','-ascii','-double');
%% 求内交点坐标 
n = size(x,1);
intersection_x = zeros(n, n);  %交点坐标x的值
intersection_y = zeros(n, n);  %交点坐标y的值
%直线斜截式方程系数
k = (y1 - y) ./ (x1 - x);
b0 = (x.*y1 - x1.*y) ./ (x - x1);
%求任意两条线段的交点
for i = 1:n
    for j = i+1 : n
        x0 = ( b0(j) - b0(i) ) ./ ( k(i) - k(j) );
        y0 = ( b0(i)*k(j) - b0(j)*k(i) ) / (k(j) - k(i));
        if (x0 >= max([x(i) x(j)])) && (x0 <= min([x1(i) x1(j)]))&&...
                (y0 >= min([y(i) y1(i)])) && (y0 <= max([y(i) y1(i)]))&&...
                (y0 >= min([y(j) y1(j)])) && (y0 <= max([y(j) y1(j)]))
            intersection_x(j,i) = x0;   %输出交点X坐标 下三角
            intersection_y(j,i) = y0;   %输出交点Y坐标 下三角
        end
    end
end
intersection_x = intersection_x + intersection_x';
intersection_y = intersection_y + intersection_y';

%%求边界交点
k = (y1 - y) ./ (x1 - x);
b0 = (x.*y1 - x1.*y) ./ (x - x1);
%左边界
Left_border_x = zeros(1, n); Left_border_y = zeros(1, n);
x0 = a;
y0 = k*x0 + b0;
for i = 1: n
    if x(i) == x0
        Left_border_x(1,i) = x0;
        Left_border_y(1,i) = y0(i);
    end
end

%右边界
Right_border_x = zeros(1, n); Right_border_y = zeros(1, n);
x0 = b;
y0 = k*x0 + b0;
for i = 1: n
    if x1(i) == x0
        Right_border_x(1,i) = x0;
        Right_border_y(1,i) = y0(i);
    end
end

%下边界
Down_border_x = zeros(1, n); Down_border_y = zeros(1, n);
y0 = c;
x0 = (y0 - b0) ./ k;
for i = 1: n
    if (( y(i) == y0 )&&(x(i)~= a))||(( y1(i) == y0 )&&(x1(i)~= b))
        
      Down_border_x(1,i) = x0(i);        
      Down_border_y(1,i) = y0;
      
    end
end

%上边界
Up_border_x = zeros(1,n); Up_border_y = zeros(1,n);
y0 = d;
x0 = (y0 - b0) ./ k;
for i = 1: n
    if (( y(i) == y0 )&&(x(i)~= a))|| (( y1(i) == y0 )&&(x1(i)~= b))
      Up_border_x(1,i) = x0(i);        
      Up_border_y(1,i) = y0;
    end
end

%%总交点矩阵
Total_intersection_x = [intersection_x; Left_border_x; Right_border_x;...
                           Down_border_x; Up_border_x];
Total_intersection_y = [intersection_y; Left_border_y; Right_border_y;...
                           Down_border_y; Up_border_y];

% save('Total intersection_x2.txt','Total_intersection_x','-ascii','-double');
% save('Total intersection_y2.txt','Total_intersection_y','-ascii','-double');

%%剔除死端裂隙
Total_intersection_xx = Total_intersection_x;
Total_intersection_yy = Total_intersection_y;
for i = 1:size(Total_intersection_x,2)
    if length( find(Total_intersection_x(:,i)) ) == 1 || ...
            (length( find(Total_intersection_x(:,i)) ) == 0)%裂隙交点个数为1，满足条件
        Total_intersection_xx(i,:) = inf; Total_intersection_xx(:,i) = inf;
        Total_intersection_yy(i,:) = inf; Total_intersection_yy(:,i) = inf;
        Connected_fracture(i,:) = inf;
        totalorder(i) = inf;
    end
end
Total_intersection_x(:,all(Total_intersection_xx==inf,1)) = [];
Total_intersection_x(all(Total_intersection_xx==inf,2),:) = [];
Total_intersection_y(:,all(Total_intersection_yy==inf,1)) = [];
Total_intersection_y(all(Total_intersection_yy==inf,2),:) = [];
Connected_fracture(all(Connected_fracture==inf,2),:) = [];
totalorder(totalorder==inf) = [];
Effective_fracture = Connected_fracture;
%重新对主干裂隙中的裂隙进行编号
totalorder2 = zeros(1,size(Connected_fracture,1));
for i = 1:size(Connected_fracture,1)
    for j = 1:length(all_x)
        if isequal(Connected_fracture(i,1:4),[all_x(j) all_y(j) all_x1(j) all_y1(j)])
            totalorder2(i) = j;
            break
        end
    end
end
% 生成非连通裂隙
Unconnected_fracture = setdiff([all_x all_y all_x1 all_y1],Connected_fracture(:,1:4),'rows');  
Unconnected_fracture = [Unconnected_fracture,zeros(size(Unconnected_fracture,1),1)];
for j = 1:size(Unconnected_fracture,1)
    for i = 1:length(all_x)
        if all_x(i) == Unconnected_fracture(j,1)&&all_y(i) == Unconnected_fracture(j,2)&&...
                all_x1(i) == Unconnected_fracture(j,3)&&all_y1(i) == Unconnected_fracture(j,4)
            Unconnected_fracture(j,5) = i;
            break
        end
    end
end

% save('Effective fracture3.txt','Effective_fracture','-ascii','-double');

% x = Connected_fracture(:,1); y = Connected_fracture(:,2);
% x1 = Connected_fracture(:,3); y1 = Connected_fracture(:,4);
% 
% n = size(x,1);
% 
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)');
% for i = 1:n
%     hold on;
%     plot([x(i),x1(i)],[y(i),y1(i)],'k');
% end

%% 绘制最终裂隙计算图
disp('正在计算最终裂隙网络');
%对同一条裂隙x坐标从小到大排列,
frac_point_2 = [];   %为裂隙的两个端点坐标矩阵
for i = 1:size(Connected_fracture,1)
    if isempty(frac_point_2)
        frac_point_2 = [Connected_fracture(i,1:2);Connected_fracture(i,3:4)];
    else
        frac_point_2 = [frac_point_2 [Connected_fracture(i,1:2);Connected_fracture(i,3:4)]];
    end
end
%   
Intersection_XY2 = [];
for i = 1 : size(Total_intersection_x,2)
    
    Intersection_XY(:,2*i-1:2*i) = [...
        Total_intersection_x(:,i) Total_intersection_y(:,i)];
    
    Intersection_XY(:,2*i-1:2*i) = sortrows(...
                            Intersection_XY(:,2*i-1:2*i));
                        
end
Intersection_XY2 = [Intersection_XY;frac_point_2];   %将裂隙两个端点放在矩阵最下面
for i = 1:size(Total_intersection_x,2)
    Intersection_XY2(:,2*i-1:2*i) = sortrows(Intersection_XY2(:,2*i-1:2*i));
end

%生成死端的线单元
SD_frac = zeros(size(Intersection_XY2,2),6);
for k = 1:2:size(Intersection_XY2,2)
        kkkkk = Intersection_XY2(Intersection_XY2(:,k) ~= 0,k:k+1);
        kkkkk = roundn(kkkkk,-4);
        kkkkk = unique(kkkkk,'rows');
        SD_frac(k,1:4) = [kkkkk(1,:) kkkkk(2,:)];
        SD_frac(k+1,1:4) = [kkkkk(size(kkkkk,1)-1,:) kkkkk(size(kkkkk,1),:)];
end
% 死端裂隙所属裂隙编号
for i = 1:length(totalorder2)
    SD_frac(2*i-1:2*i,6) = totalorder2(i);
end
%将裂隙、节点、线单元编号、隙宽储存到一个矩阵中
e = 0; 
Information = zeros(length(find(Total_intersection_x)),9);
for i = 1 : 2 : size(Intersection_XY,2)
    for j = 1 : size(Intersection_XY,1)
        if Intersection_XY(j,i) ~= 0
            e = e + 1;
            Information(e,1) = Intersection_XY(j,i); 
            Information(e,2) = Intersection_XY(j,i+1);
            Information(e,3) = (i+1)/2;   %裂隙编号
            Information(e,9) = totalorder2((i+1)/2);  %所属裂隙的原始编号
        end
    end
end

N1 = 0;   %内交点个数
for i = 1:size(Information,1)-1
    u = 0;
    if Information(i,4) ~= 0
        continue
    end
    N1 = N1 + 1;
    Information(i,4) = N1;               %节点编号
    for j = i+1 : size(Information,1)
        if Information(i,1) == Information(j,1)&&...
           Information(i,2) == Information(j,2)
               if Information(i,1) == a
                   Information(i,4) = -a;
                   Information(j,4) = -a;
                   continue
               else if Information(i,1) == b
                       Information(i,4) = -b;
                       Information(j,4) = -b;
                       continue
                   end
               end
       
               Information(j,4) = Information(i,4);
               u = 1;
               break
        end
    end
    if u == 0
       Information(i,4) = 0;
       N1 = N1 - 1;
    end
end
intersection_xy = [Information(Information(:,4)>0&Information(:,4)<=N1,1:2)];

N2 = 0; %第二类边界点个数
N3 = 0; %第一类边界点个数
for i = 1:size(Information,1)
    if (Information(i,4) == 0)&&...
       (Information(i,2) == c && (Information(i,1)~=a && Information(i,1)~=b)) ||...
        (Information(i,2) == d &&(Information(i,1)~=a && Information(i,1)~=b))
    
            N2 = N2 + 1;
            Information(i,4) = N1 + N2;
            
    end
end

for i = 1:size(Information,1)
    if (Information(i,4) == 0)&&...
       (Information(i,1) == a)
   
            N3 = N3 + 1;
            Information(i,4) = N1+N2+N3;
            
    end
end

n3 = N3;
for i = 1:size(Information,1)
    if (Information(i,4) == 0)&&...
       (Information(i,1) == b)
   
            N3 = N3 + 1;
            Information(i,4) = N1+N2+N3;
            
    end
end
%        
% figure;
% axis equal;
% axis([a b c d]);
% xlabel('(m)');
% ylabel('(m)');
h = 0;
for i = 1:size(Information,1)
    
    if i == size(Information,1)
        break
    end
    
    if Information(i,3) == Information(i+1,3)
        
        h = h + 1;
        Information(i,5) = h;               %线单元编号
        Information(i,6) = sqrt((Information(i,1) - ...
            Information(i+1,1))^2 + (Information(i,2)...
            - Information(i+1,2))^2);           %线单元长度
        Information(i,7) = Effective_fracture(...
            Information(i,3),5);            %线单元宽度
%         hold on;
%         if isempty(final_network)
%             final_network = [Information(i,1),Information(i,2)];
%         else
%             final_network = [final_network;Information(i,1),Information(i,2)];
%         end
%         final_network = plot([Information(i,1),Information(i+1,1)]...
%             ,[Information(i,2),Information(i+1,2)],'k');
    end
end

%求衔接矩阵
Join = zeros(max(Information(:,5)),N1+N2+N3);
for i = 1:max(Information(:,5))
    for j = 1:size(Information,1)
        if Information(j,5) == i
            Join(i,Information(j,4)) = -1;
            Join(i,Information(j+1,4)) = 1;
        end
    end
end
Join = Join';

A1 = Join(1:N1,:);
A2 = Join(N1+1:N1+N2,:);
A3 = Join(N1+N2+1:N1+N2+N3,:);
disp('衔接矩阵计算结束');
%求对角阵  （20℃）
Part = Information(:,5:7);
Part(all(Part==0,2),:) = [];
Diagonal = (0.809e-3)*Part(:,3).^3 ./ Part(:,2);
T = diag(Diagonal);

%%渗流计算（稳定流）  Aq + Q = 0
disp('开始渗流计算');
Q1 = zeros(N1,1); Q2 = zeros(N2,1);  %节点源汇项为零
H3 = zeros(N3,1); %一类节点边界水头
H3(1:n3) = Value_left;
H3(n3+1:N3) = Value_right;

G = [A1*T*A1' A1*T*A2'; A2*T*A1' A2*T*A2'];
Q = [Q1+A1*T*A3'*H3; Q2+A2*T*A3'*H3];
H = -inv(G) * Q;
H1 = H(1:N1);
H2 = H(N1+1:N1+N2);
Q3 = -A3*T*A1'*H1 - A3*T*A2'*H2 - A3*T*A3'*H3;
flow = sum(abs(Q3))/2;                         %总流量
K = (flow * (b-a)) / ((Value_left - Value_right)*(d - c));  %渗透系数
disp('渗流计算结束');
Level = zeros(size(H1,1),3);
H_total = [H;H3];


for i = 1 : size(H1,1)
    for j = 1 : size(Information,1)
        if i == Information(j,4)
            Level(i,1:2) = Information(j,1:2);
            Level(i,3) = H1(i);
            Information(j,8) = H1(i);
        end
    end
end
for i = 1:size(H2,1)
    for j = 1:size(Information,1)
        if i + N1 == Information(j,4)
            Information(j,8) = H2(i);
        end
    end
end
Information(Information(:,1) == 10,8)= Value_left;
Information(Information(:,1) == 70,8)= Value_right;
Information(:,1:2) = roundn(Information(:,1:2),-4);
for i = 1:size(SD_frac,1)
    Circle_num_SD = i
    for j = 1:size(Information,1)
        if length(intersect(SD_frac(i,1:4),Information(j,1:2))) == 2
            SD_frac(i,5) = Information(j,8);
            if SD_frac(i,1) == a
                SD_frac(i,5) = Value_left;
            elseif SD_frac(i,3) == b
                SD_frac(i,5) = Value_right;
            continue
            end
        end
    end
end

intersection_U_C = zeros(size(neighbor_oo,1),4*size(Unconnected_fracture,1));
ZG_frac = [];
for i = 1:size(Information,1)
    if i == size(Information,1)
        break
    end
    if Information(i,3) == Information(i+1,3)
        ZG_frac = [ZG_frac;Information(i,1),Information(i,2),Information(i,8),...
            Information(i+1,1),Information(i+1,2),Information(i+1,8),Information(i,9),Information(i,5)];
    end
end
ZG_frac = roundn(ZG_frac,-4);
% 将主干裂隙中存在的死端裂隙从死端裂隙矩阵中去除
SD_ZG_frac = [];
n = 1;
for i = 1:size(SD_frac,1)
    for j = 1:size(ZG_frac,1)
        if SD_frac(i,6) == ZG_frac(j,7)
            if isequal(SD_frac(i,1:4),ZG_frac(j,[1 2 4 5]))
                SD_ZG_frac(n,1) = i;
                n = n + 1;
                break
            end
        end
    end
end
SD_frac(SD_ZG_frac,:) = [];
% 计算非连通裂隙与主干裂隙的交点
for i = 1:size(Unconnected_fracture,1)
    for j = 1:size(ZG_frac,1)
        if neighbor_oo(Unconnected_fracture(i,5),ZG_frac(j,7)) == 1
            k1 = (Unconnected_fracture(i,4) - Unconnected_fracture(i,2)) / (Unconnected_fracture(i,3) - Unconnected_fracture(i,1));
            b1 =(Unconnected_fracture(i,1) * Unconnected_fracture(i,4) - Unconnected_fracture(i,3) *...
                Unconnected_fracture(i,2))/(Unconnected_fracture(i,1) - Unconnected_fracture(i,3));
            k2 = (ZG_frac(j,5) - ZG_frac(j,2))/(ZG_frac(j,4) - ZG_frac(j,1));
            b2 = (ZG_frac(j,1) * ZG_frac(j,5) - ZG_frac(j,4) * ZG_frac(j,2))/(ZG_frac(j,1) - ZG_frac(j,4));
            x0 = ( b2 - b1 ) / ( k1 - k2 );
            y0 = ( b1*k2 - b2*k1 ) / (k2 - k1);
            if (x0 >= max([Unconnected_fracture(i,1) ZG_frac(j,1)])) && (x0 <= min([Unconnected_fracture(i,3) ZG_frac(j,4)]))&&...
                (y0 >= min([Unconnected_fracture(i,2) Unconnected_fracture(i,4)])) && (y0 <= max([Unconnected_fracture(i,2) Unconnected_fracture(i,4)]))&&...
                (y0 >= min([ZG_frac(j,2) ZG_frac(j,5)])) && (y0 <= max([ZG_frac(j,2) ZG_frac(j,5)]))
            intersection_U_C(ZG_frac(j,7),4*i-3:4*i-2) = [x0 y0];
            H_yy = ZG_frac(j,3) + (ZG_frac(j,6) - ZG_frac(j,3)) * norm([x0 y0]-[ZG_frac(j,1) ZG_frac(j,2)])/norm([ZG_frac(j,1) ZG_frac(j,2)]-[ZG_frac(j,4) ZG_frac(j,5)]);
            intersection_U_C(ZG_frac(j,7),4*i-1) = H_yy;
            end
        end
    end
end
for i = 1:size(Unconnected_fracture,1)
    for j = 1:size(SD_frac,1)
        if neighbor_oo(Unconnected_fracture(i,5),SD_frac(j,6)) == 1
            k1 = (Unconnected_fracture(i,4) - Unconnected_fracture(i,2)) / (Unconnected_fracture(i,3) - Unconnected_fracture(i,1));
            b1 =(Unconnected_fracture(i,1) * Unconnected_fracture(i,4) - Unconnected_fracture(i,3) *...
                Unconnected_fracture(i,2))/(Unconnected_fracture(i,1) - Unconnected_fracture(i,3));
            k2 = (SD_frac(j,4) - SD_frac(j,2))/(SD_frac(j,3) - SD_frac(j,1));
            b2 = (SD_frac(j,1) * SD_frac(j,4) - SD_frac(j,3) * SD_frac(j,2))/(SD_frac(j,1) - SD_frac(j,3));
            x0 = ( b2 - b1 ) / ( k1 - k2 );
            y0 = ( b1*k2 - b2*k1 ) / (k2 - k1);
            if (x0 >= max([Unconnected_fracture(i,1) SD_frac(j,1)])) && (x0 <= min([Unconnected_fracture(i,3) SD_frac(j,3)]))&&...
                (y0 >= min([Unconnected_fracture(i,2) Unconnected_fracture(i,4)])) && (y0 <= max([Unconnected_fracture(i,2) Unconnected_fracture(i,4)]))&&...
                (y0 >= min([SD_frac(j,2) SD_frac(j,4)])) && (y0 <= max([SD_frac(j,2) SD_frac(j,4)]))
                        intersection_U_C(SD_frac(j,6),4*i-3:4*i-2) = [x0 y0];
                        H_yy = SD_frac(j,5);
                        intersection_U_C(SD_frac(j,6),4*i-1) = H_yy;
            end
        end
    end
end
for i = 1:size(Unconnected_fracture,1)
    for j = 1:size(Unconnected_fracture,1)
        if i == j
            continue
        end
        if neighbor_oo(Unconnected_fracture(i,5),Unconnected_fracture(j,5)) == 1
            k1 = (Unconnected_fracture(i,4) - Unconnected_fracture(i,2)) / (Unconnected_fracture(i,3) - Unconnected_fracture(i,1));
            b1 = (Unconnected_fracture(i,1) * Unconnected_fracture(i,4) - Unconnected_fracture(i,3) *...
                Unconnected_fracture(i,2))/(Unconnected_fracture(i,1) - Unconnected_fracture(i,3));
            k2 = (Unconnected_fracture(j,4) - Unconnected_fracture(j,2)) / (Unconnected_fracture(j,3) - Unconnected_fracture(j,1));
            b2 = (Unconnected_fracture(j,1) * Unconnected_fracture(j,4) - Unconnected_fracture(j,3) *...
                Unconnected_fracture(j,2))/(Unconnected_fracture(j,1) - Unconnected_fracture(j,3));
            x0 = ( b2 - b1 ) / ( k1 - k2 );
            y0 = ( b1*k2 - b2*k1 ) / (k2 - k1);
            if (x0 >= max([Unconnected_fracture(i,1) Unconnected_fracture(j,1)])) && (x0 <= min([Unconnected_fracture(i,3) Unconnected_fracture(j,3)]))&&...
                (y0 >= min([Unconnected_fracture(i,2) Unconnected_fracture(i,4)])) && (y0 <= max([Unconnected_fracture(i,2) Unconnected_fracture(i,4)]))&&...
                (y0 >= min([Unconnected_fracture(j,2) Unconnected_fracture(j,4)])) && (y0 <= max([Unconnected_fracture(j,2) Unconnected_fracture(j,4)])) 
            intersection_U_C(Unconnected_fracture(j,5),4*i-3:4*i-2) = [x0 y0];
            intersection_U_C(Unconnected_fracture(j,5),4*i) = 1;
            end
        end
    end
end
Unconnected_point = [];
for i = 1:size(Unconnected_fracture,1)
    if isempty(Unconnected_point)
        Unconnected_point = [Unconnected_fracture(i,1:2) 0 1;Unconnected_fracture(i,3:4) 0 1];
    else
        Unconnected_point = [Unconnected_point [Unconnected_fracture(i,1:2) 0 1;Unconnected_fracture(i,3:4) 0 1]];
    end
end
intersection_U_C = [intersection_U_C;Unconnected_point];
for i = 1:size(Unconnected_fracture,1)
    intersection_U_C(:,4*i-3:4*i) = sortrows(intersection_U_C(:,4*i-3:4*i));
end
 %建立非连通裂隙线单元矩阵
UC_Information = [];
for i = 1 : 4 : size(intersection_U_C,2)
    for j = 1 : size(intersection_U_C,1)
        if intersection_U_C(j,i) ~= 0
            UC_Information = [UC_Information;intersection_U_C(j,i:i+3) (i+3)/4 0 ] ; 
        end
    end
end


h = 0;
for i = 1:size(UC_Information,1)
    if i == size(UC_Information,1)
        break
    end
    if UC_Information(i,5) == UC_Information(i+1,5)
        h = h + 1;
        UC_Information(i,6) = h;               %线单元编号
    end
end
%左右边界点
for i = 1:size(UC_Information,1)
    if UC_Information(i,1) == a
        UC_Information(i,3) = Value_left;
        UC_Information(i,4) = 0;
    else if UC_Information(i,1) == b
            UC_Information(i,3) = Value_right;
            UC_Information(i,4) = 0;
        end
    end
end

UC_Information = [UC_Information zeros(size(UC_Information,1),1)];
n1 = 0;   %内交点个数
for i = 1:size(UC_Information,1)-1
    u = 0;
    if UC_Information(i,7) ~= 0
        continue
    end
    n1 = n1 + 1;
    UC_Information(i,7) = n1;               %节点编号
    for j = i+1 : size(UC_Information,1)
        if UC_Information(i,1) == UC_Information(j,1)&&...
           UC_Information(i,2) == UC_Information(j,2)
               if UC_Information(i,1) == a
                   UC_Information(i,7) = -a;
                   UC_Information(j,7) = -a;
                   continue
               else if UC_Information(i,1) == b
                       UC_Information(i,7) = -b;
                       UC_Information(j,7) = -b;
                       continue
                   end
               end
       
               UC_Information(j,7) = UC_Information(i,7);
               u = 1;
               break
        end
    end
    if u == 0
       UC_Information(i,7) = 0;
       n1 = n1 - 1;
    end
end
disp('内交点编号结束');
n2 = 0; %第二类边界点个数
n3 = 0; %第一类边界点个数
for i = 1:size(UC_Information,1)
    if (UC_Information(i,7) == 0)&&...
       (UC_Information(i,2) == c && (UC_Information(i,1)~=a && UC_Information(i,1)~=b)) ||...
        (UC_Information(i,2) == d &&(UC_Information(i,1)~=a && UC_Information(i,1)~=b))
    
            n2 = n2 + 1;
            UC_Information(i,7) = n1 + n2;
            
    end
end

for i = 1:size(UC_Information,1)
    if (UC_Information(i,7) == 0)&&...
       (UC_Information(i,1) == a)
   
            n3 = n3 + 1;
            UC_Information(i,7) = n1+n2+n3;
            
    end
end

for i = 1:size(UC_Information,1)
    if (UC_Information(i,7) == 0)&&...
       (UC_Information(i,1) == b)
   
            n3 = n3 + 1;
            UC_Information(i,7) = n1+n2+n3;
            
    end
end


%计算与主干裂隙相连的节点之间的其他节点的水头
for i = 1:size(UC_Information,1)
    for j = i+1:size(UC_Information,1)
        if UC_Information(i,5) == UC_Information(j,5)&&UC_Information(i,4) == 0&&UC_Information(j,4) == 0
            if j - i > 1
                x1 = UC_Information(i,1);
                y1 = UC_Information(i,2);
                x2 = UC_Information(j,1);
                y2 = UC_Information(j,2);
                for k = i+1:j-1
                    if UC_Information(k,4) == 1
                        x3 = UC_Information(k,1);
                        y3 = UC_Information(k,2);
                        UC_Information(k,3) = UC_Information(i,3) + (UC_Information(j,3) - UC_Information(i,3)) * norm([x1 y1]-[x3 y3])/norm([x1 y1]-[x2 y2]);
                        UC_Information(k,4) = 0;
                    end
                end
            end
        end
    end
end
%内节点水头计算
for i = 1:size(UC_Information,1)
    if UC_Information(i,4) == 1
        continue
    end
    for j = 1:size(UC_Information,1)
        if UC_Information(i,7) <= n1&&UC_Information(j,7) <= n1&&UC_Information(i,7) > 0&&UC_Information(i,7) > 0&&...
                UC_Information(i,7) == UC_Information(j,7)&&UC_Information(j,4) == 1
            UC_Information(j,3) = UC_Information(i,3);
            UC_Information(j,4) = 0;
        end
    end
end
%内节点和与主干裂隙及边界点相连节点之间其他节点水头
while 1
    Change_num = 0;
    while 1
        change_num = 0;
        for i = 1:size(UC_Information,1)
            for j = i+1:size(UC_Information,1)
                if UC_Information(i,5) == UC_Information(j,5)&&UC_Information(i,4) == 0&&UC_Information(j,4) == 0
                    if j - i > 1
                        x1 = UC_Information(i,1);
                        y1 = UC_Information(i,2);
                        x2 = UC_Information(j,1);
                        y2 = UC_Information(j,2);
                        for k = i+1:j-1
                            if UC_Information(k,4) == 1
                                x3 = UC_Information(k,1);
                                y3 = UC_Information(k,2);
                                UC_Information(k,3) = UC_Information(i,3) + (UC_Information(j,3) - UC_Information(i,3)) * norm([x1 y1]-[x3 y3])/norm([x1 y1]-[x2 y2]);
                                UC_Information(k,4) = 0;
                                change_num = change_num + 1;
                                Change_num = Change_num + 1;
                            end
                        end
                    end
                end
            end
        end
        %再次计算内节点水头
        for i = 1:size(UC_Information,1)
            if UC_Information(i,4) == 1
                continue
            end
            for j = 1:size(UC_Information,1)
                if UC_Information(i,7) <= n1&&UC_Information(j,7) <= n1&&UC_Information(i,7) > 0&&UC_Information(i,7) > 0&&...
                        UC_Information(i,7) == UC_Information(j,7)&&UC_Information(j,4) == 1
                    UC_Information(j,3) = UC_Information(i,3);
                    change_num = change_num + 1;
                    Change_num = Change_num + 1;
                    UC_Information(j,4) = 0;
                end
            end
        end
        if change_num == 0
            break
        end
    end
    %特殊节点计算，两端节点有水头，中间没有
    while 1
        change_num = 0;
        for i = 2:size(UC_Information,1)-1
            if UC_Information(i,4) == 1&&UC_Information(i,7) > 0&&UC_Information(i,7) <= n1&&...
                    UC_Information(i-1,4)*UC_Information(i+1,4) == 0&&UC_Information(i-1,4)+UC_Information(i+1,4) ~= 0
                for j = 2:size(UC_Information,1)-1
                    if UC_Information(j,7) == UC_Information(i,7)&&UC_Information(j-1,4)*UC_Information(j+1,4) == 0&&...
                            UC_Information(j-1,4)+UC_Information(j+1,4) ~= 0
                        uc_information = [UC_Information(i-1,:);UC_Information(i+1,:);UC_Information(j-1,:);UC_Information(j+1,:)];
                        uc_information = sortrows(uc_information,3);
                        UC_Information(i,3) = uc_information(3,3) + (uc_information(4,3) - uc_information(3,3))*...
                            norm(uc_information(3,1:2)-UC_Information(i,1:2))/(norm(uc_information(3,1:2)-UC_Information(i,1:2))+norm(uc_information(4,1:2)-UC_Information(i,1:2)));
                        UC_Information(j,3) = UC_Information(i,3);
                        UC_Information(j,4) = 0;
                        change_num = change_num + 1;
                        Change_num = Change_num + 1;
                    end
                end
            end
        end
        for i = 1:size(UC_Information,1)
            if UC_Information(i,4) == 1
                continue
            end
            for j = 1:size(UC_Information,1)
                if UC_Information(i,7) > 0&&UC_Information(j,7) > 0&&UC_Information(i,7) <= n1&&UC_Information(j,7) <= n1&&...
                        UC_Information(i,7) == UC_Information(j,7)&&UC_Information(j,4) == 1
                    UC_Information(j,3) = UC_Information(i,3);
                    UC_Information(j,4) = 0;
                    change_num = change_num + 1;
                    Change_num = Change_num + 1;
                end
            end
        end
        if change_num == 0
            break
        end
    end
    if Change_num == 0
        break
    end
end
%最后在判断是否有内节点未赋予水头
% for i = 2:size(UC_Information,1)-1
%     if UC_Information(i,4) == 1&&UC_Information(i,7) > 0&&UC_Information(i,7) <= n1
%         if UC_Information(i-1,4) == 0&&UC_Information(i+1,4) == 0&&UC_Information(i,5) == UC_Information(i-1,5)&&UC_Information(i,5) == UC_Information(i+1,5)
%             x1 = UC_Information(i-1,1);
%             y1 = UC_Information(i-1,2);
%             x2 = UC_Information(i+1,1);
%             y2 = UC_Information(i+1,2);
%             x3 = UC_Information(i,1);
%             y3 = UC_Information(i,2);
%             UC_Information(i,3) = UC_Information(i-1,3) + (UC_Information(i+1,3) - UC_Information(i-1,3)) * ...
%                 norm([x1 y1]-[x3 y3])/norm([x1 y1]-[x2 y2]);
%             UC_Information(i,4) = 0;
%         end
%         for j = 2:size(UC_Information,1)-1
%             if UC_Information(i,7) == UC_Information(j,7)
%                 if UC_Information(i,4) == 1&&UC_Information(j,4) == 1&&UC_Information(j-1,4) == 0&&...
%                         UC_Information(j+1,4) ==0&&UC_Information(j,5) == UC_Information(j-1,5)&&...
%                         UC_Information(j,5) == UC_Information(j+1,5)
%                     x1 = UC_Information(j-1,1);
%                     y1 = UC_Information(j-1,2);
%                     x2 = UC_Information(j+1,1);
%                     y2 = UC_Information(j+1,2);
%                     x3 = UC_Information(j,1);
%                     y3 = UC_Information(j,2);
%                     UC_Information(j,3) = UC_Information(j-1,3) + (UC_Information(j+1,3) - UC_Information(j-1,3)) *...
%                         norm([x1 y1]-[x3 y3])/norm([x1 y1]-[x2 y2]);
%                     UC_Information(i,3) = UC_Information(j,3);
%                     UC_Information(i,4) = 0;
%                     UC_Information(j,4) = 0;
%                 elseif UC_Information(i,4) == 0&&UC_Information(j,4) == 1
%                     UC_Information(j,3) = UC_Information(i,3);
%                     UC_Information(j,4) = 0;
%                 end
%             end
%         end
%     end
% end
  
%死端裂隙
for i = 1:size(UC_Information,1)
    if UC_Information(i,4) == 1
        if i == 1
            UC_Information(i,3) = UC_Information(i+1,3);
            UC_Information(i,4) = 0;
        elseif i ~= 1&&i ~= size(UC_Information,1)
            if UC_Information(i,5) == UC_Information(i-1,5)
                UC_Information(i,3) = UC_Information(i-1,3);
                UC_Information(i,4) = 0;
            elseif UC_Information(i,5) == UC_Information(i+1,5)
                UC_Information(i,3) = UC_Information(i+1,3);
                UC_Information(i,4) = 0;
            end
        elseif i == size(UC_Information,1)
            UC_Information(i,3) = UC_Information(i-1,3);
            UC_Information(i,4) = 0;
        end
    end
end

fid = fopen('trend_example_e.dat','a');
fid2 = fopen('trend_example_main_e.dat','a');
fprintf(fid,'TITLE="head distribution"\n');
fprintf(fid2,'TITLE="head distribution"\n');
for i = 1:size(Information,1)
    if i == size(Information,1)
        break
    end
    if Information(i,3) == Information(i+1,3)
        
        fprintf(fid,'VARIABLES="X","Y","head"\n');
        fprintf(fid,'ZONE T="ZONE %g"\n',i);
        fprintf(fid,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
        fprintf(fid,'%g\t%g\t%g\n%g\t%g\t%g\n',Information(i,1),Information(i,2),Information(i,8),...
            Information(i+1,1),Information(i+1,2),Information(i+1,8));
        fprintf(fid2,'VARIABLES="X","Y","head"\n');
        fprintf(fid2,'ZONE T="ZONE %g"\n',i);
        fprintf(fid2,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
        fprintf(fid2,'%g\t%g\t%g\n%g\t%g\t%g\n',Information(i,1),Information(i,2),Information(i,8),...
            Information(i+1,1),Information(i+1,2),Information(i+1,8));
    end
end

for i = 1:size(SD_frac,1)
    if i == size(Information,1)
        break
    end
    fprintf(fid,'VARIABLES="X","Y","head"\n');
    fprintf(fid,'ZONE T="ZONE %g"\n',i + size(Information,1));
    fprintf(fid,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
    fprintf(fid,'%g\t%g\t%g\n%g\t%g\t%g\n',SD_frac(i,1),SD_frac(i,2),SD_frac(i,5),...
        SD_frac(i,3),SD_frac(i,4),SD_frac(i,5));
    fprintf(fid2,'VARIABLES="X","Y","head"\n');
    fprintf(fid2,'ZONE T="ZONE %g"\n',i + size(Information,1));
    fprintf(fid2,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
    fprintf(fid2,'%g\t%g\t%g\n%g\t%g\t%g\n',SD_frac(i,1),SD_frac(i,2),SD_frac(i,5),...
        SD_frac(i,3),SD_frac(i,4),SD_frac(i,5));
end
for i = 1:size(UC_Information,1)
    if i == size(UC_Information,1)
        break
    end
    if UC_Information(i,5) == UC_Information(i+1,5)
        
        fprintf(fid,'VARIABLES="X","Y","head"\n');
        fprintf(fid,'ZONE T="ZONE %g"\n',i + size(Information,1) + size(SD_frac,1));
        fprintf(fid,'I=2, J=1, K=1, ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE)\n');
        fprintf(fid,'%g\t%g\t%g\n%g\t%g\t%g\n',UC_Information(i,1),UC_Information(i,2),UC_Information(i,3),...
            UC_Information(i+1,1),UC_Information(i+1,2),UC_Information(i+1,3));
    end
end
fclose(fid);
fclose(fid2);
Head_distribution = Information(:,[1 2 8]);
for i = 1:size(SD_frac,1)
    Head_distribution = [Head_distribution;SD_frac(i,[1 2 5]);SD_frac(i,[3 4 5])];
end
for i = 1:size(UC_Information,1)
    Head_distribution = [Head_distribution;UC_Information(i,1:3)];
end
Head_distribution = roundn(Head_distribution,-4);
Head_distribution = unique(Head_distribution,'rows');
figure
X = 10:0.1:30;
Y = 10:0.1:30;
[X1,Y1] = meshgrid(X,Y);
Z1 = griddata(Head_distribution(:,1),Head_distribution(:,2),Head_distribution(:,3),X1,Y1,'v4');
contourf(X1,Y1,Z1);
% 裂隙水流方向修正，衔接矩阵正负
% Head_difference_cc = [];
% nn = 1;
% Head_difference=zeros(size(Part,1),1);
% error = [];
% n = 1;
% for j=1:size(Part,1)
%     for i=1:size(Information,1)-1
%         if Information(i,5)==0
%             continue
%         elseif Information(i,5)==j
%             Point_number_left=Information(i,4);
%             Point_number_right=Information(i+1,4);
%             Head_difference(j,1) = H_total(Point_number_left) - H_total(Point_number_right);
%             if length(Head_difference_cc) ~= 0
%                 Head_difference(Head_difference_cc) = -Head_difference(Head_difference_cc);
%             end
%             if Head_difference(j,1) < 0
%                 Head_difference_cc = [Head_difference_cc;j];
%                 error(n,1) = Point_number_left;
%                 error(n,2) = Point_number_right;
%                 n = n + 1;
%             end
%         end
%     end
% end
% for i = 1:size(error,1)
%     ccc1 = find(Information(:,4) == error(i,1));
%     ccc2 = find(Information(:,4) == error(i,2));
%     ccc = intersect(Information(ccc1,5),Information(ccc2-1,5));
%     Join(error(i,1),ccc) = -Join(error(i,1),ccc);
%     Join(error(i,2),ccc) = -Join(error(i,2),ccc);
% end

% final_network(final_network(:,1)==10,3)=20;
% [X,Y]=meshgrid(min(Level(:,1)):max(Level(:,1)),min(Level(:,2)):max(Level(:,2)));
% Z = griddata(Level(:,1),Level(:,2),Level(:,3),X,Y);
% 数据导出
% X = reshape(X,p * q,1);
% Y = reshape(Y,p * q,1);
% Z = reshape(Z,p * q,1);
% Z1 = Z;
% Maaaa = [X Y Z];
% Maaaa(isnan(Z),:) = [];
% DT = delaunay(Maaaa(:,1),Maaaa(:,2));
% m = size(DT,1);
% p = size(Maaaa,1);
% fid = fopen('test_head2.dat','a');
% fprintf(fid,'TITLE="head distribution"\n');
% fprintf(fid,'VARIABLES="X","Y","head"\n');
% fprintf(fid,'ZONE NODES=%g,ELEMENTS=%g,DATAPACKING=POINT,ZONETYPE=FETRIANGLE\n',p,m);
% dlmwrite('test_head2.dat',Maaaa,'-append','delimiter','\t');
% dlmwrite('test_head2.dat',DT,'-append','delimiter','\t');
% fclose(fid);

% cloud=contour(X,Y,Z);
% contour_h = findobj(gca,'Type','hggroup');
% clabel(cloud);
% uistack(cloud,'bottom');
