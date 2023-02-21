%% 
cub_size = 2;
num = [3000;1500];
x = 11:cub_size:69;
y = 11:cub_size:49;
Px = 3 * x - 80;
Py = -0.025*(y-40).^2+10;
Px(Px < 0) = 0;
Nx = num .* Px/sum(Px);
zouxjz = [20 70];
zouxbzc = [5 5];
jicjz = [3 3];
jicbzc = [1 1];
xikuan = 1;
frac_data = [];
for I = 1:2
    for i = 1:size(x,2)
        Ny = Nx(I,i) * Py/sum(Py);
        for j = 1:size(y,2)
            if round(Ny(j)) == 0
                continue
            end
            ny = round(Ny(j));
            X0 = x(1,i) + cub_size/2 * unifrnd(-1,1,ny,1);
            Y0 = y(1,j) + cub_size/2 * unifrnd(-1,1,ny,1);
            zouxiang = normrnd(zouxjz(I),zouxbzc(I),ny,1);
            l = xikuan * ones(ny,1);
            m = jicjz(I);
            v = jicbzc(I) ^ 2;
            mu = log((m^2)/sqrt(v + m^2));
            sigma = sqrt(log(v/(m^2) + 1));
            jichang = lognrnd(mu,sigma,ny,1);
            x1 = X0 + (jichang / 2) .* cos(zouxiang * pi/180);
            y1 = Y0 + (jichang / 2) .* sin(zouxiang * pi/180);
            x2 = X0 + (jichang / 2) .* cos((zouxiang + 180) * pi/180);
            y2 = Y0 + (jichang / 2) .* sin((zouxiang + 180) * pi/180);
            frac_data = [frac_data;X0 Y0 zouxiang jichang l];
            for k = 1:ny
                plot([x1(k),x2(k)],[y1(k),y2(k)],'k');
                hold on
                daspect([1 1 1]);
                axis([10 70 10 50]);
            end
        end
    end
end
dlmwrite('Trend_distribution_fracdata.txt',frac_data,'-append','delimiter','\t');
%% 
cub_size = 2;
x = 11:cub_size:69;
y = 11:cub_size:49;
fx1 = 3.49*x - 15.12; fx2 = 1.61*x - 0.83;
Nx = [fx1;fx2];
fy1 = -0.69*y.^2 + 41.97*y - 349.56; fy2 = -0.34*y.^2 + 20.66*y - 175.83;
Ny = [fy1;fy2];
zouxjz = [20 70];
zouxbzc = [5 5];
jicjz = [3 3];
jicbzc = [1 1];
xikuan = 1;
frac_data = [];
for I = 1:2
    for i = 1:size(x,2)
        Num = Nx(I,i) * Ny(I,:)/sum(Ny(I,:));
        Num(Num<0) = 0;
        for j = 1:size(y,2)
            if round(Num(j)) == 0
                continue
            end
            ny = round(Num(j));
            X0 = x(1,i) + cub_size/2 * unifrnd(-1,1,ny,1);
            Y0 = y(1,j) + cub_size/2 * unifrnd(-1,1,ny,1);
            zouxiang = normrnd(zouxjz(I),zouxbzc(I),ny,1);
            l = xikuan * ones(ny,1);
            m = jicjz(I);
            v = jicbzc(I) ^ 2;
            mu = log((m^2)/sqrt(v + m^2));
            sigma = sqrt(log(v/(m^2) + 1));
            jichang = lognrnd(mu,sigma,ny,1);
            x1 = X0 + (jichang / 2) .* cos(zouxiang * pi/180);
            y1 = Y0 + (jichang / 2) .* sin(zouxiang * pi/180);
            x2 = X0 + (jichang / 2) .* cos((zouxiang + 180) * pi/180);
            y2 = Y0 + (jichang / 2) .* sin((zouxiang + 180) * pi/180);
            frac_data = [frac_data;X0 Y0 zouxiang jichang l];
            for k = 1:ny
                plot([x1(k),x2(k)],[y1(k),y2(k)],'k');
                hold on
                daspect([1 1 1]);
                axis([10 70 10 50]);
            end
        end
    end
end
dlmwrite('Trend_distribution_fracdata.txt',frac_data,'-append','delimiter','\t');