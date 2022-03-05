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
    if B(i)<=1
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
k2 = 0:0.1:2*pi;
for i = 1:1:length(aalpha)
    k1 = 0:0.1/Q(i):2*pi/(Q(i));
    for ja = 1:1:length(k2)
        for l1 = 1:1:length(k1)
            h = zeros(Q(i),Q(i));
            Hf = zeros(2*Q(i),2*Q(i));
            for j = 1:2:2*Q(i)-1
                Hf(j,j+1) = (exp(1i*(k2(ja)*a))+exp(-1i*(2*pi*(aalpha(i))*(j+1)/2)));
            end
            for j1 = 2:2:2*Q(i)-1
                Hf(j1,j1+1) = exp(-1i*(2*pi*(aalpha(i))*(j1+2)/2))*exp(-1i*(k1(l1)*a));
            end
            Hf(1,2*Q(i)) = Hf(1,2*Q(i))+exp(1i*(k1(l1)*a));
            L = eig(Hf+ctranspose(Hf));
            for ii = 1:1:2*Q(i)
                counter = counter+1;
                EE(counter) = L(ii);
                MF(counter) = aalpha(i);
            end
        end
     end
end
figure(1);
hold on
plot(MF,(50/3)*EE,'.b','MarkerSize',4)
xlabel('p/q')
ylabel('E[meV]');
set(gca,'Fontsize',20);
box on
hold off
% save('Rammal.mat','EE','MF')
timeelapsed = toc