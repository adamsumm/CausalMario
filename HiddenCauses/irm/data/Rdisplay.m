function Rdisplay(R, datacc, varargin)

% show R sorted by class assignments in DATACC

lw= 0.05;	% linewidth
args=varargin;
linesflag=1;
ticksflag=0;
labels=0;
featcc=[];
for i=1:2:length(args)
  switch args{i}
   case 'lines', linesflag=args{i+1};
   case 'ticks', ticksflag=args{i+1};
   case 'labels', labels=args{i+1};
   case 'rectangle', featcc=args{i+1}; 
  end
end

n=size(R,1);
colormap(gray)
numdata=length(datacc);
[s,sind]=sort(datacc);
changes=find(s(2:end)-s(1:end-1))+0.5;
cnum=length(changes);

if ~isempty(featcc) % use separate clustering for the columns
  fnumdata=length(featcc);
  [fs,fsind]=sort(featcc);
  fchanges=find(fs(2:end)-fs(1:end-1))+0.5;
  fcnum=length(fchanges);
else
  fnumdata = numdata; fchanges = changes; fcnum = cnum; fsind = sind;
end


imagesc(1-R(sind, fsind));
axis equal
axis tight

if linesflag 
  hold on
  plot([0.5*ones(1, cnum); fnumdata*ones(1, cnum)+0.5],...
     [changes; changes], 'b-', 'linewidth', lw); 
  plot([fchanges; fchanges], ...
     [0.5*ones(1, fcnum); numdata*ones(1, fcnum)+0.5], 'b-', 'linewidth', lw); 
end

if ticksflag
  set(gca, 'xtick', changes, 'ticklength', [0.02 0.02]);
  set(gca, 'ytick', changes)
  set(gca, 'tickdir', 'out')
else
  set(gca, 'xtick', [], 'ytick', []);
end

if labels==0
  set(gca, 'xticklabel', []);
  set(gca, 'yticklabel', []);
end

hold off

set(gca, 'XAxisLocation', 'top');
box off

