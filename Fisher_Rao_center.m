function [w] = Fisher_Rao_center(W)
w = mean(W);
for i=1:10
u = W - (W*w') * w;
v = bsxfun(@times, u, acos( (W*w') ./ sqrt(sum(u.*u, 2))));
v = 0.5* mean(v);
w = cos(abs(v)).*w + sin(abs(v)).*v ./abs(v);
w = w / sum(w);
end
end