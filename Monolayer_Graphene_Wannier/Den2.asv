function [add,EF1] =Den2(p_input,q_input,pts_input,MN)
p = p_input;
q = q_input;
alpha =p_input/q_input;

pts=pts_input;
% pts_alpha = ceil(pts/q);
Enespec(:,:,:) = MLG_Ham(p,q,pts);
k1max = 2*pi/q;
k2max = 2*pi; 
K1 = linspace(0,k1max,pts);
K2 = linspace(0,k2max,pts);

dim = 2*q;
a = ceil(dim/2);

count = 0;
for i = 1:1:dim
    a1 = max(max(Enespec(:,:,i)));
    count = count + 1;
    Ebq(:,:,count) = Enespec(:,:,i);
end
[k11,k22]= meshgrid(K1(1):(K1(end)-K1(1))/30:K1(end),K2(1):(K2(end)-K2(1))/30:K2(end));
k1 = K1(1):(K1(end)-K1(1))/99:K1(end);
k2 = K2(1):(K2(end)-K2(1))/99:K2(end);
for i = 1:1:length(Ebq(1,1,:))
    Eb(:,:,i) = interp2(K1,K2,Ebq(:,:,i),k11,k22,'cubic',0);
end
count1 = 0;
EG = [];
for i = 1:1:length(Eb(1,1,:))-1
    Gap = min(min(Eb(:,:,i+1)))-max(max(Eb(:,:,i)));
        count1 = count1 + 1;
        EG(count1,:) =  linspace(min(min(Eb(:,:,i+1))),max(max(Eb(:,:,i))),length(k1)*length(k2));
end
        EF1 = sort(reshape(Eb,[],1));
        EF2 = sort(reshape(EG,[],1));
        Lmax = length(EF1)+length(EF2);
        add = zeros(1, 2*MN);
        if Lmax>MN
            for i = 2:1:MN-1
                j = ceil(i*length(EF1)/MN);
                add(i) = EF1(j);
            end
            for i = 2:1:MN-1
                j = ceil(i*length(EF2)/MN);
                add(MN+i) = EF2(j);
            end
            add(1) = EF1(1);
            add(MN) = EF1(end);
            add(MN+1) = EF2(1);
            add(2*MN) = EF2(end);
        elseif Lmax<=MN
            for i = 1:1:length(EF1)
                add(i) = EF1(i);
            end
            for i = 1:1:length(EF2)
                add(i+length(EF1)) = EF2(i);
            end
            manage = 2*MN-(length(EF1)+length(EF2));
            EF3 = linspace(min(EF1),max(EF1),manage);
            for i = 1:1:length(EF3)
                add(i+length(EF1)+length(EF2)) = EF3(i);
            end

        end
add = sort(add);
