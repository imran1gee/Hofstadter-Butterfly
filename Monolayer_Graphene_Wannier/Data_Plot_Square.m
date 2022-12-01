clear all
tic
data1 = load('Square.mat','addE','EF1E');
Add = data1.addE;
F1 = data1.EF1E;
Data = load('Square_Data.mat','Q','P');
Q = Data.Q;
P = Data.P;
alpha = P./Q;
figure(2)
for i1 = 1:1:length(Add(:,1))
    count2 = 0;
    add = Add(i1,:);
    g = 10^-16;
    for i = 1:1:length(add)
        EF1 = F1{i1};
        count2 = count2 +1;
        D(count2) = g/(length(EF1))*sum(1./((add(i)-EF1(:)).^2+g^2));
    end
    DE(i1,:)=log(D)/norm(log(D));
    clear add
end
N1 = -1:2/(2000-1):1;
N = N1/max(N1);
alpha = P(1:length(Add(:,1)))./Q(1:length(Add(:,1)));
[c, h] = contourf(alpha,N,DE',unique(DE'));
colormap('hot')
colorbar
set(h,'LineColor','none')
toc