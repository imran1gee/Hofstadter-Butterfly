tic
clear all
qmax = 15;
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
aalpha(:) = C2(:); 

counter = 0;
a = 1;
w = 1;
k2 = -pi:0.1:pi;
for i = 1:1:length(aalpha)
    k1 = -pi/(Q(i)):0.1/Q(i):pi/(Q(i));
    for ja = 1:1:length(k2)
        for l1 = 1:1:length(k1)
            Hf = zeros(Q(i),Q(i));
            for j = 1:1:Q(i)
                Hf(j,j) = 1*cos(k2(ja)-2*pi*j*aalpha(i));%2 is not because we will add them
            end
            for j = 1:1:Q(i)-1
                Hf(j,j+1) = exp(1i*k1(l1));
            end
            Hf(1,Q(i)) = Hf(1,Q(i)) + 1*exp(-1i*k1(l1));
            L = eig(Hf+ctranspose(Hf));
            for ii = 1:1:Q(i)
                counter = counter+1;
                EE(counter) = L(ii);
                MF(counter) = aalpha(i);
            end
        end
     end
end
figure(1);
hold on
plot(MF,EE,'.k','MarkerSize',6)
xticks([0,1/4,1/2,3/4,1])
xticklabels({'0','1/3','1/2','2/3','1'})
ylabel('E[meV]')
xlabel('p/q')
set(gca,'FontSize',20)
box on
hold off
% save('Rammal.mat','EE','MF')
timeelapsed = toc