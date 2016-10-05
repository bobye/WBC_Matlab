function [c] = centroid_init(stride, supp, w, options)
  dim = size(supp,1);
  n = length(stride);
  m = length(w);
% start guess
if ~isfield(options, 'support_points')
  c=centroid_rand(stride, supp, w, options);
  support_size=length(c.w);
  save(['cstart' num2str(n) '-' num2str(support_size) '.mat'], 'c', 'support_size');
  %fprintf('Options: support points of centroid are supposed to be optimized!');
elseif ~isempty(options.support_points)
  c.supp = options.support_points;
  support_size=size(c.supp,2);
  c.w = ones(1,support_size)/support_size;
  %fprintf('Options: support points of centroid are supposed to be fixed!');
else
  error('Please provide support points at options.support_points (if not empty), otherwise leave it undefined');
end