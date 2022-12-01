function Eigcolumn =Square_Ham(p_input,q_input,pts_input)



x = -p_input/q_input:(1/(pts_input-1))*2*pi/q_input:pi/q_input;
y = 0:(1/(pts_input-1))*2*pi:2*pi;
Eigcolumn = zeros(length(x),length(y),q_input);
for i = 1:1:length(x)
    for t = 1:1:length(y)
        Eigcolumn(i,t,:) = eig(Ham(x(i),y(t),p_input,q_input));
    end
end


    function H = Ham(x,y,p,q)
        H = zeros(q,q);
        for j = 1:q:q
            H(j,j) = 1*cos(y-2*pi*j*p/q);
        end
        for j = 1:1:q-1
            H(j,j+1) = exp(1i*x);
        end
        H(1,q) = H(1,q)+exp(-1i*x);
        H = H+ctranspose(H);
    end
end