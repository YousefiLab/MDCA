function [Y_k_d] = down_Y_k(Y_k, down_rate)

m = size(Y_k);

Y_k_d = zeros(floor(m(1)./down_rate), m(2));

count = 0;
for i=1: down_rate: m(1)
    count = count + 1;
    Y_k_d(count, :) = Y_k(i, :);
end

end