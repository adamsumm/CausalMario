function printrelsystem(data,outfile, varargin)

% convert the relations in DATA into a form that can be used by the IRM (a file
% with one line per weighted edge)

args=varargin;
missingflags = []; 
contflag = 0;

for i=1:2:length(args)
  switch args{i}
    case 'missingflags',    missingflags=args{i+1};	    % unobserved edges?
    case 'contflag', contflag=args{i+1};		    % continuous data?
  end
end

nreln = length(data);
if length(missingflags) == 0
  missingflags = zeros(1, nreln);
end

fid = fopen(outfile,'w');

for i = 1:nreln
  G = data{i};
  if contflag
    nonans = find(~isnan(G));			% missing data marked with nans
    edgendx = ind2subv(size(G), nonans)-1; 
    for j = 1:size(edgendx,1)
      fprintf(fid,'%5d ', i-1);
      fprintf(fid,'%5d ', edgendx(j,:));
      fprintf(fid,'%5g\n', G(nonans(j)));
    end
  else
    edgendx = ind2subv(size(G), find(G==1))-1; 
    for j = 1:size(edgendx,1)
      fprintf(fid,'%5d ', i-1);
      fprintf(fid,'%5d ', edgendx(j,:));
      fprintf(fid,' 1\n');
    end
    if missingflags(i) 
      edgendx = ind2subv(size(G), find(G==0))-1; 
      for j = 1:size(edgendx,1)
        fprintf(fid,'%5d ', i-1);
        fprintf(fid,'%5d ', edgendx(j,:));
        fprintf(fid,' 0\n');
      end
    end
  end
end
fclose(fid);
