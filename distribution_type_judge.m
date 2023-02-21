data = textread('C:\Users\Administrator\Desktop\FracPaQ2DEZROSE.txt');
data(:,2) = [];
data = Angle2;
% m = size(data,1);
% dist = ones(m);
% delt_z = ones(m);
% for i = 1:m
%     for j = i + 1:m
%         dist(i,j) = sqrt((data(i,1) - data(j,1))^2 + (data(i,2) - data(j,2))^2);
%         dist(j,i) = dist(i,j);
%         delt_z(i,j) = data(i,3) - data(j,3);
%         delt_z(j,i) = -delt_z(i,j);
%     end
% end
% gradient = delt_z ./ dist;
% connected = zeros(m);
% connected(dist == 5) = 1;
% gradient_con = gradient .* connected;
N=3;  %Cluster number
[m,n]=size(data);
pattern=zeros(m,n+1);
center=zeros(N,1);
pattern(:,1:n)=data(:,:);
for x=1:N
    center(x)=data( randi(m,1),1);
end
a = 1;
while 1   
    distence=zeros(1,N);
    num=zeros(1,N);
    new_center=zeros(N,1);
    
    for x=1:m
        for y=1:N
            distence(y)=abs(data(x,1)-center(y));
        end
        [~, temp]=min(distence);
        pattern(x,n+1)=temp;
    end
    k=0;
    for y=1:N   
        for x=1:m
            if pattern(x,n+1)==y
                new_center(y)=new_center(y)+pattern(x,1);
                num(y)=num(y)+1;
            end
        end
        new_center(y)=new_center(y)/num(y);  
        if abs(new_center(y)-center(y))<1  
            k=k+1;
        end
    end
    if k>=3   
        break;
    else
        center=new_center;
    end
    a = a + 1;
    if a > 1000
        break
    end
end
[m, n]=size(pattern);


p1 = []; p2 = []; p3 = []; p4 = [];
% figure;
% hold on;
for i=1:m
    if pattern(i,n) == 1
        p1 = [p1;pattern(i,1)];
%         p1=plot(pattern(i,1),pattern(i,2),'r*');
    elseif pattern(i,n) == 2
        p2 = [p2;pattern(i,1)];
%         p2=plot(pattern(i,1),pattern(i,2),'g*');
    elseif pattern(i,n) == 3
        p3 = [p3;pattern(i,1)];
%         p3=plot(pattern(i,1),pattern(i,2),'b*');
    elseif pattern(i,n) == 4
        p4 = [p4;pattern(i,1)];
%         p4=plot(pattern(i,1),pattern(i,2),'y*');
    else
%         p5=plot(pattern(i,1),pattern(i,2),'m*');
    end
end
grid on;
axis([10 70 10 50]);
daspect([1 1 1]);
legend([p1,p2,p3],['Cluster1£º',num2str(center(1))],['Cluster2£º',num2str(center(2))],...
    ['Cluster3£º',num2str(center(3))],'location','NorthOutside')
X = 10:1:70;
Y = 10:1:50;
[X1,Y1] = meshgrid(X,Y);
Z1 = griddata(pattern(:,1),pattern(:,2),pattern(:,3),X1,Y1,'v4');
contour(X1,Y1,Z1);