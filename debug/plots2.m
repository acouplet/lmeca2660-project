
T = csvread('data/P_400_600.csv');
%load('T.mat');
size(T)
T = T(1:265*602,:);
A = reshape(T,602,402,265);
size(A)

contourf(A(:,:,100))
