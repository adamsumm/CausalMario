rand('state', 0);
rand('seed', 0);
randn('state', 0);
randn('seed', 0);


n = 40;
classnum = 5;
z = repmat(1:classnum, n/classnum,1);
z = z(:);

G = zeros(n);

cols = repmat(z, 1, n);
rows = repmat(z', n, 1);

for i = 1:classnum
  for j = 1:classnum
    inds = (cols == i) & (rows == j);
    cellcount = sum(inds(:));
    mu = rand - 0.5;
    sigma = 0.05;
    vals = random('normal', mu, sigma,  sum(inds(:)), 1);
    G(inds) = vals;
  end
end

G = G - mean(G(:));
G = G./std(G(:));

% IRM code doesn't allow self links
G(sub2ind(size(G), 1:n, 1:n)) = nan; 

colormap(gray);

% because Rdisplay uses 1-R
imagesc(1-G);


printrelsystem({G}, 'cont1.graph', 'contflag', 1)
save('cont1.mat', 'G');

