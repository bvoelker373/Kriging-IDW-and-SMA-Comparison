%--------------------------------------------------------------------------
%  Analysis and Comparison of Three Different Spatial Interpolation Methods
%     - Kriging
%     - Inverse Distance Weighted
%     - Spatial Moving Average
%
%  Brandon Voelker
%  05/06/2020
%--------------------------------------------------------------------------

[t]=geotiffread('elev_ned_30m.tif');
figure(5); imagesc(t); colorbar;

%% Create Semivariogram Scatter Plot

clc, clear, close all

elevData = imread('elev_ned_30m.tif', 'tif');

[eRow, eCol] = size(elevData);

%compute bandwidth
pixel = 30; %meters
h = (sqrt(eRow^2+eCol^2))/2;
h_meter = pixel*h;

%randomly select 5000 points
        
rng default
randRow = randi(eRow,1,5000); randRow = randRow';
rng default
randCol = randi(eCol,1,5000); randCol = randCol';

randPoint = [randRow,randCol];

%divide h into number of bins of size 10
numMidpoints = ceil(h/10);

n=-5;
itCheck = 0;
sqDiff = 0;

for i = 1:numMidpoints
%each i represents an increment of 10 cells (300m) 
%so the midpoint of the bins will be 5,15,25, etc. (150m,450m,750m, etc.)
    n=n+10;
    itCheck=itCheck+1;
    row=-3;
    for j = 1:length(randCol)
       %within this loop is the same distance, but changing the point that
       %it is centred around
       j;
        row=row+4;
       %calculate square difference from original point of
       %elevData(randRow(j),randCol(j))
       
       %four points: 
       %elevData( randRow(j)+n,randCol(j) )
       %elevData( randRow(j),randCol(j)+n )
       %elevData( randRow(j)-n,randCol(j) )
       %elevData( randRow(j),randCol(j)-n )

       
       %find squared difference from each point
       
       try
            sqDiff(row,i) = ( ...
                elevData(randRow(j)+n,randCol(j))-...
                elevData(randRow(j),randCol(j))...
                )^2;
       catch
           sqDiff(row,i) = NaN;
       end
     
       try
            sqDiff(row+1,i) = ( ...
                elevData(randRow(j),randCol(j)+n)-...
                elevData(randRow(j),randCol(j))...
                )^2;
       catch
           sqDiff(row+1,i) = NaN;
       end
       
       try
            sqDiff(row+2,i) = ( ...
                elevData(randRow(j)-n,randCol(j))-...
                elevData(randRow(j),randCol(j))...
                )^2;
       catch
           sqDiff(row+2,i) = NaN;
       end
       
       try
            sqDiff(row+3,i) = ( ...
                elevData(randRow(j),randCol(j)-n)-...
                elevData(randRow(j),randCol(j))...
                )^2;
       catch
           sqDiff(row+3,i) = NaN;
       end
        
    end
end

%calculate means for each column of square distances, where each column
%represents a different midpoint, then talk half of it to find gamma
G = (nanmean(sqDiff))/2;

%create vector of distances in meters (150,300,450, etc)
for i = 1:numMidpoints
    d(i) = 30*(10*(i-1)+5);
    d_check(i) = 10*(i-1)+5;
end

figure(1)
scatter(d,G)
hold on
%add labels


%% Linear Least-Squares Regression of Semivariogram

H = d'; H(:,2) = 1;
l = G';
x = (inv(H'*H))*(H'*l);

slope = x(1); b = x(2);
fplot(@(x)slope*x+b, [min(d) max(d)])
%add labels
hold off


%% Select 25,000 Random Locations to Exclude (Without repeats)

excludeRow = randi(eRow,1,25000); excludeRow = excludeRow';
rng(1)
excludeCol = randi(eCol,1,25000); excludeCol = excludeCol';

excludedPoints = [excludeRow,excludeCol];

%scatter(excludedPoints(:,1),excludedPoints(:,2))

repeats = 1;
count = 0;
while repeats == 1
    repeats = 0;
    for i = 1:(25000-1) %rows
        for j = (i+1):25000
            if excludedPoints(i,1) == excludedPoints(j,1)
                if excludedPoints(i,2) == excludedPoints(j,2)
                    repeats =1;
                    excludedPoints(j,1) = randi(eRow);
                    count = count + 1; 
                end
            end
        end
    end
end

%% Simple Kriging


%check = [];
%check2=[]
%k_skip = 0;
iterIJ = 0;
krig_variance =[];

for i = 1:eRow
    for j = 1:eCol
        %(i,j) is the location of the point being interpolated
        k=[];
        K=[];
        coords = [];
        iterNM = 1;
        
        
        for n = i-5:i+5
            if n >= 1 && n <= 450
                for m = j-5:j+5
                    %(n,m) is the location of m(u)
                    if m >= 1 && m <= 500
                        
                        distFromCtr = sqrt(    (n-i)^2+(m-j)^2   );
                        if distFromCtr <= 5 %neighborhood size
                            flag = 0;
                            for ex_Chk = 1:length(excludeRow)
                                if n == excludedPoints(ex_Chk,1) && m == excludedPoints(ex_Chk,2)
                                    flag = 1;
                                    %check(end+1) = n;
                                    %check2(end+1)=m;
                                    %k_skip = k_skip + 1;
                                end
                                if n == i && m == j
                                    flag = 1;
                                end
                        
                            end
                            if flag ==0
                                k(end+1) = (slope*max(d)+b) - (b+slope*distFromCtr);
                                coords(iterNM,1) = n;
                                coords(iterNM,2) = m;
                                iterNM = iterNM + 1;
                            end
                        end
                        
                    
                    end
                end
            end
        end; iterIJ = iterIJ + 1;
        
        
        for p = 1:length(coords)
            for q = 1:length(coords)
                distBetween = sqrt( (coords(p,1) - coords(q,1))^2 + (coords(p,2) - coords(q,2))^2 );
                %distCheck(q) = distBetween;
                K(p,q) = (slope*max(d)+b) - (b+slope*distBetween);
            end
        end
        
        k = k';
        lambda = (inv(K))*k;
        
        for p = 1:length(coords)
            %elevData(coords(p,1),coords(p,2));
            zIntermed(p) = lambda(p) * elevData(coords(p,1),coords(p,2));
            z(i,j) = sum(zIntermed);
        end
        krig_variance(i,j) = (max(G)) -(lambda' * k);
        
        
    end
end

%% Kriging Variance
%kv = (max(G)) -(lambda' * k);
kv_mean = mean(mean(krig_variance));


%% IDW

check = [];
check2=[];
iterIJ = 0;

IDW_Q = 12; %set power for IDW


for i = 1:eRow
    for j = 1:eCol
        %(i,j) is the location of the point being interpolated
        inv_Q_dist = [];
        lambda_IDW=[];
        coords = [];
        iterNM = 1;
        
        
        for n = i-5:i+5
            if n >= 1 && n <= 450
                for m = j-5:j+5
                    %(n,m) is the location of m(u)
                    if m >= 1 && m <= 500
                        
                        distFromCtr = sqrt(    (n-i)^2+(m-j)^2   );
                        if distFromCtr <= 5 %neighborhood size
                            flag = 0;
                            for ex_Chk = 1:length(excludeRow)
                                if n == excludedPoints(ex_Chk,1) && m == excludedPoints(ex_Chk,2)
                                    flag = 1;
                                    check(end+1) = n;
                                    check2(end+1)=m;
                                end
                                if n == i && m == j
                                    flag = 1;
                                end
                        
                            end
                            if flag ==0
                                distances(n,m) = distFromCtr;
                                %calculate distance to -q power
                                
                                coords(iterNM,1) = n;
                                coords(iterNM,2) = m;
                                
                                inv_Q_dist(end+1) = distFromCtr^(-IDW_Q);
                                
                                iterNM = iterNM + 1;
                                
                            end
                        end
                        
                    
                    end
                end
            end
        end; iterIJ = iterIJ + 1;
        
        
        for p = 1:length(inv_Q_dist)
            lambda_IDW(p) = (inv_Q_dist(p)) / (sum(inv_Q_dist));
        end
        
        for q = 1:length(coords)
            zIntermed_IDW(q) = lambda_IDW(q) * elevData(coords(q,1),coords(q,2));
            z_IDW(i,j) = sum(zIntermed_IDW); 
        end
        z_IDW(i,j);
        
    end
end


%% Spatial Moving Average

%check = [];
%check2=[];
iterIJ = 0;

for i = 1:eRow
    for j = 1:eCol
        %(i,j) is the location of the point being interpolated
        iterNM = 1;
        nearVals = [];
        
        
        for n = i-5:i+5
            if n >= 1 && n <= 450
                for m = j-5:j+5
                    %(n,m) is the location of m(u)
                    if m >= 1 && m <= 500
                        
                        distFromCtr = sqrt(    (n-i)^2+(m-j)^2   );
                        if distFromCtr <= 5 %neighborhood size
                            flag = 0;
                            for ex_Chk = 1:length(excludeRow)
                                if n == excludedPoints(ex_Chk,1) && m == excludedPoints(ex_Chk,2)
                                    flag = 1;
                                    %check(end+1) = n;
                                    %check2(end+1)=m;
                                end
                                if n == i && m == j
                                    flag = 1;
                                end
                            end
                            
                            if flag ==0
                                
                                nearVals(end+1) = elevData(n,m);
                                
                                iterNM = iterNM + 1;
                                
                            end
                        end
                        
                    
                    end
                end
            end
        end; iterIJ = iterIJ + 1;
        
        z_SMA(i,j) = mean(nearVals);
        
    end
end


%% Plots

%imwrite2tif(elevData,[],'original.tif','uint8');

%imwrite2tif(z,[],'kriging.tif','uint8');

figure(7)
imagesc(elevData)
colormap jet
caxis manual
caxis([50 160]);
colorbar

figure(8)
imagesc(z)
colormap jet
caxis manual
caxis([50 160]);
colorbar

figure(9)
imagesc(z_IDW)
colormap jet
caxis manual
caxis([50 160]);
colorbar

figure(10)
imagesc(z_SMA)
colormap jet
caxis manual
caxis([50 160]);
colorbar

%% Compare Check Points

for i = 1:length(excludedPoints)
    MSE_Krig_Vector(i) = ( ( elevData(excludedPoints(i,1),excludedPoints(i,2)) ) - ...
        ( z(excludedPoints(i,1),excludedPoints(i,2)) ) )^2;
    
    MSE_IDW_Vector(i) = ( ( elevData(excludedPoints(i,1),excludedPoints(i,2)) ) - ...
        ( z_IDW(excludedPoints(i,1),excludedPoints(i,2)) ) )^2;
    
    MSE_SMA_Vector(i) = ( ( elevData(excludedPoints(i,1),excludedPoints(i,2)) ) - ...
        ( z_SMA(excludedPoints(i,1),excludedPoints(i,2)) ) )^2;
end
MSE_Krig = ( sum(MSE_Krig_Vector) ) / ( length(excludedPoints) );
MSE_IDW = ( sum(MSE_IDW_Vector) ) / ( length(excludedPoints) );
MSE_SMA = ( sum(MSE_SMA_Vector) ) / ( length(excludedPoints) );

%% MSE Comparison
%MSE_IDW_M = [82.043472, 18.700159, 3.4592149, 1.1762522, 0.82894462, ...
%    0.73972589, 0.70540792, 0.69005781, 0.68275946, 0.67916983, 0.67735839, 0.67642498];
%MSE_IDW_M_2 = [3.4592149, 1.1762522, 0.82894462, ...
%    0.73972589, 0.70540792, 0.69005781, 0.68275946, 0.67916983, 0.67735839, 0.67642498];
%MSE_IDW_M_3 = [1.1762522, 0.82894462, ...
%    0.73972589, 0.70540792, 0.69005781, 0.68275946, 0.67916983, 0.67735839, 0.67642498];
%MSE_IDW_M_4 = [ ...
%    0.73972589, 0.70540792, 0.69005781, 0.68275946, 0.67916983, 0.67735839, 0.67642498];
%figure(20)
%hold on

%% Statistical Tests
tTest = ttest2(MSE_Krig_Vector,MSE_IDW_Vector);
tTest2 = ttest2(MSE_SMA_Vector,MSE_IDW_Vector);
tTest3 = ttest2(MSE_Krig_Vector,MSE_SMA_Vector);
