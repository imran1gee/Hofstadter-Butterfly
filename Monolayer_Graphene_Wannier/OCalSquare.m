clear all
qmax = 30;
for i = 1:1:qmax                                             
    for j = 1:1:qmax                                         
        A(i,j) = i/j;                                        
    end                                                      
end 
B(:) = unique(reshape(A',1,[]));                             
for i = 1:1:length(B)                                        
    if B(i)<1
        C1(i) = B(i);
    end
end                                                          
C2 = nonzeros(C1');                                          
[p,q1] = numden(sym(C2));                                                                           
Q(:) = double(q1(:));                                        
P(:) = double(p(:)); 
save('Square_Data.mat','Q','P');
pts = 50;
MN = 1000;
for i = 1:length(Q)
    tic
    fprintf('Progress %d of %d... \n',i,length(Q))
    [add,EF1] = DenSquare(P(i),Q(i),pts,MN);
    addE(i,:) = add(:);
    EF1E{i} = EF1(:);
    alpha(i) = P(i)/Q(i);
    toc
end 
alpha = sort(alpha);
save('Square.mat','addE','EF1E')