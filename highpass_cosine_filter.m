function Y=highpass_cosine_filter(tdim,Y,dcut,TR)

%dcut is the inverse of the cutoff frequency, usually dcut=120sec
R=floor(2*tdim*TR/dcut+1); %the highest cosine function
t=1:tdim;
X=zeros(tdim,R);  %design matrix
for r=1:R
    X(:,r)=cos(r*pi*t/tdim);
end
X(:,R+1)=1;  %constant regressor

%regression
beta=inv(X'*X)*X'*Y;
Y=Y-X(:,1:R)*beta(1:R,:);

%disp('done with highpass_cosine_filter');