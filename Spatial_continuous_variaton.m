Length = [10 70];
Wide = [10 50];
cube = 2;
Directionmean = [40 115];
DirectionSD = [5 5];
tracelengthmean = [2 2];
tracelengthSD = [1 1];
data = load('Sampledata.txt');
S = data(:,1:2);
theta = [10 10]; lob = [0.1 0.1]; upb = [20 20];
XY = gridsamp([Length(1)+cube/2 Wide(1)+cube/2;Length(2)-cube/2 Wide(2)-cube/2],[30,20]);
X1 = XY(:,1);
Y1 = XY(:,2);
figure
for k = 1:2      
    Z = data(:,k+2);
    [dmodel,perf] = dacefit(S,Z,@regpoly0,@corrgauss,theta,lob,upb);
    [Z1,MSE] = predictor(XY,dmodel);
    Z1(Z1<0) = 0;
    for i = 1:length(Z1)
        num = round(Z1(i));       
        X0 = unifrnd(X1(i)-cube/2,X1(i)+cube/2,num,1);
        Y0 = unifrnd(Y1(i)-cube/2,Y1(i)+cube/2,num,1);
        Direction = normrnd(Directionmean(k),DirectionSD(k),num,1);   
        m = tracelengthmean(k);
        v = tracelengthSD(k)^2;
        mu = log((m^2)/sqrt(v + m^2));
        sigma = sqrt(log(v/(m^2) + 1));
        tracelength = lognrnd(mu,sigma,num,1);
        x1 = X0 + (tracelength / 2) .* cos(Direction * pi/180);
        y1 = Y0 + (tracelength / 2) .* sin(Direction * pi/180);
        x2 = X0 + (tracelength / 2) .* cos((Direction + 180) * pi/180);
        y2 = Y0 + (tracelength / 2) .* sin((Direction + 180) * pi/180);
        aperture = ones(num,1);
        frac_data = [frac_data;X0 Y0 Direction tracelength aperture];
        for j = 1:num
            plot([x1(j),x2(j)],[y1(j),y2(j)],'k');
            hold on
            daspect([1 1 1]);
            axis([Length(1) Length(2) Wide(1) Wide(2)]);
        end
    end
end
