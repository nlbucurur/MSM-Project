clear all;
diary log.txt;

% PROBLEM: From an extensive list of mobile phone models with some of their
% properties and respective prices we want to carry out PCA in order to
% compress the number of the phone properties into a set of new P variables
% which could be associated to their final price. As a secondary objective
% we would like our code to be quite general, so that the input considered
% could be changed without having to do much tuning in the code.

% JUSTIFICATION: Obtaining such a set of reduced Pi variables which are
% able to separate phones according to their price (in the Pi space) could
% be useful for example in the prediction of the ideal price point for a new
% phone releasing, using a lower number of variables would reduce computation
% costs during the training and allow for the training of more complex algorithms.

% the raw data taken from https://www.kaggle.com/datasets/rockyjoseph/mobile-
% price-prediction (see mobile_dataset.csv)

% the dataset was cleaned for PCA: reference removed, display and processor
% columns removed, price column formatted and converted from rupees to euros
% (see mobile_dataset_cleaned.csv)

% extract the cleaned data from the CSV file
data = readmatrix('mobile_dataset_cleaned.csv','Range',[2,2]);

% since the price of the phones is our target variable we will not consider
% it to carry out the PCA process, therefore we create our data matrix M
% without this variable
M = data(:,1:end-1);

% size of M
m = size(M,1);
n = size(M,2);

% create the matrix with standarized data Ms : For a value x_i of the
% variable x its standarized form is given by x_is = (x_i - mean(x))/std(x)
Ms = zeros(m,n);
for j = 1:n
    for i = 1:m
        Ms(i,j) = (M(i,j) - mean(M(:,j)))/std(M(:,j));
    end
end

% we find the covariance matrix V which is given by V = tX*X/(m-1) where 
% X is the data matrix and m the number of data entries
V = ((Ms.')*(Ms))/(m - 1);

% in this case, since we used standarized data the covariance matrix
% corresponds to the correlation matrix, we check that it has a correct 
% form(size should be nxn, elements in diagonal should be 1, other elements
% between -1 and 1 and it should be symmetric)
nxn = true;
diagonal_1 = true;
range_nondiagonal = true;
symmetric = true;

if size(V,1) ~= n || size(V,2) ~= n
    nxn = false;
end

for j = 1:size(V,2)
    for i = 1:size(V,1)
        if i == j
            if int32(V(i,j)) ~= 1
                diagonal_1 = false;
            end
        else
            if abs(V(i,j)) > 1
                range_nondiagonal = false;
            end
            if V(i,j) ~= V(j,i)
                symmetric = false;
            end
        end
    end
end

if nxn & diagonal_1 & range_nondiagonal & symmetric
    disp('Construction of V matrix is OK');
    disp(' ');
end
         
% we now want to obtain the eigenvectors and eigenvalues of the V matrix
% since its eigenvectors will define the Pi axes (new variables) and the
% associated eigenvalues I_Pi will be the information held by the corresponding 
% axis. E is the matrix built with each eigenvector  as a column and D is 
% the diagonal matrix built with the eigenvalues

[E,D] = eig(V);

% since we want to reduce the total number of variables used (compression)
% we now look at the percentage of information held by each axis, given by
% (I_Pi/I_total)*100. 

% vector Ip with I_Pi values
Ip = zeros(n,1);
for i = 1:n
    Ip(i) = D(i,i);
end

% vector I with information percentage held by each axis, we check that the
% sum of the elements of I is 100
I = zeros(n,1);
I_total = sum(Ip);
for i = 1:n
    I(i) = (Ip(i)/I_total)*100;
end
if sum(I) == 100
    disp('The percentage of information held by each Pi axis is:');
    disp(I);
end

% we want to keep the minimum amount of Pi axes that allow as to retain at
% least 70% of the original information
% indices vector will contain the indexes associated with the chosen Pi's
indices = zeros([],1);
Isorted = sort(I,'descend');
Iconserved = 0;
for i = 1:n
    Iconserved = Iconserved + Isorted(i);
    for j = 1:n
        if I(j) == Isorted(i)
            indices(i) = j;
        end
    end
    if Iconserved >= 70.
        break
    end
end
disp('By choosing the following Pi axes:');
disp(indices);
disp('The percentage of conserved information is:');
disp(Iconserved);

% we build the Ep (Eprime) and Dp (Dprime) matrices from the E and D matrices
% by choosing only the elements associated with the chosen Pi's axes
Dp = zeros(size(indices,2),size(indices,2));
Ep = zeros(size(E,1),size(indices,2));
for j = 1:size(Dp,2)
    for i = 1:size(Dp,1)
        if i == j
            Dp(i,j) = D(indices(1,i),indices(1,i));
        else
            Dp(i,j) = 0;
        end
    end
end
for j = 1:size(Ep,2)
        Ep(:,j) = E(:,indices(1,j));
end

% we want to get the projection of our original data (from the normalized
% matrix Ms) in the P basis (using only the selected Pi's), so that 
% P = Ms*Ep but in order to normalize the P variables we divide each Pi by
% the sqrt(I_Pi) 

% we need to create a DpNorm containing (1/sqrt(I_Pi)) elements separately
% to avoid dividing by 0 (which would happen if we were to do Dp^(-1/2)
DpNorm = Dp;
for j = 1:size(DpNorm,2)
    for i = 1:size(DpNorm,1)
        if DpNorm(i,j) ~= 0
            DpNorm(i,j) = DpNorm(i,j)^(-1/2);
        else
            continue
        end
    end
end

% now we can multiply Ep*DpNorm without problems in order to obtain our P 
% matrix         
P = Ms*Ep*DpNorm;

% we verify that the size of P is m x the number of selected Pi's
if size(P,1) == m & size(P,2) == size(indices,2)
    disp('The size of P is OK');
end

% we will now plot the Pi values for the phones in our input list to see how
% the phones with different price points are distributed in the Pi space
% (note: this part of the code could still be further generalized by
% plotting inside a for loop according to the total number of Pi's selected
% for the compression)

% plot Pi(1)vsPi(2)
markerSize = 100; 
scatter(P(:,1),P(:,2),markerSize,data(:,end),'filled');
colorbar1v2 = colorbar;
title(colorbar1v2,'Price');
colormap jet;
xlabel(sprintf('P_%i', indices(1,1)));
ylabel(sprintf('P_%i', indices(1,2)));
grid on;
fig1v2 = gcf;
exportgraphics(fig1v2, "Pi(1)vsPi(2).png");

% plot Pi(1)vsPi(3) if more than 2 Pi axes were used
if size(indices,2) > 2
    scatter(P(:,1),P(:,3),markerSize,data(:,end),'filled');
    colorbar1v3 = colorbar;
    title(colorbar1v3,'Price');
    colormap jet;
    xlabel(sprintf('P_%i', indices(1,1)));
    ylabel(sprintf('P_%i', indices(1,3)));
    grid on;
    fig1v3 = gcf;
    exportgraphics(fig1v3, "Pi(1)vsPi(3).png");
end

% plot Pi(2)vsPi(3) if more than 2 Pi axes were used
if size(indices,2) > 2
    scatter(P(:,2),P(:,3),markerSize,data(:,end),'filled');
    colorbar2v3 = colorbar;
    title(colorbar2v3,'Price');
    colormap jet;
    xlabel(sprintf('P_%i', indices(1,2)));
    ylabel(sprintf('P_%i', indices(1,3)));
    grid on;
    fig2v3 = gcf;
    exportgraphics(fig2v3, "Pi(2)vsPi(3).png");
end

% 3D plot of the distribution in the Pi space (if 3 Pi axes were used)
if size(indices,2) == 3
    scatter3(P(:,1),P(:,2),P(:,3),markerSize,data(:,end),'filled');
    colorbar3d = colorbar;
    title(colorbar3d,'Price');
    colormap jet;
    xlabel(sprintf('P_%i', indices(1,1)));
    ylabel(sprintf('P_%i', indices(1,2)));
    zlabel(sprintf('P_%i', indices(1,3)));
    fig3D = gcf;
    exportgraphics(fig3D, "PiSpace.png");
    grid on;
end

% COMMENT ON OUR RESULTS: we can see that there's a relation between
% the Pi values for the different phones with their respective price,
% therefore the reduced set of variables obtained could be in fact used as
% an indicator of a phone's price point and could be used in the training
% of a predictive algorithm.

diary off;
