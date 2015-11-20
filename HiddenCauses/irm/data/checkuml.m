cd ..;
currdir = pwd;
cd data;
scorefile    =[currdir, '/output/uml_status'];
conceptzsfile=[currdir, '/output/uml_dom1'];
relzsfile    =[currdir, '/output/uml_dom0'];

% original data
load('uml.mat');

conceptzs=load(conceptzsfile);
relzs    =load(relzsfile);
score    =load(scorefile);

% find configuration with the highest probability 
[m,mind]    =max(score(:,7));
bestassgt   =conceptzs(mind,:);
bestrelassgt=relzs(mind,:);
us	    =unique(bestassgt);

% show concept clusters
for i=1:length(us)
  groups{i} = names(bestassgt==us(i));
  groupnames{i} = gnames(bestassgt==us(i))';
  disp(['Group ', num2str(i)]);
  for j = 1:length(groups{i});
    % print concept and expert-supplied label
    disp(sprintf('%45s%30s', groups{i}{j}, groupnames{i}{j})); 
  end
  disp('*******************');
end

% show relation clusters
groups = {};
us = unique(bestrelassgt);
for i=1:length(us)
  groups{i} = relnames(bestrelassgt==us(i));
  disp(['Group ', num2str(i)]);
  for j = 1:length(groups{i});
    disp(sprintf('%45s', groups{i}{j})); 
  end
  disp('*******************');
end

% display sorted matrix for the relation 'affects'
relind = 2;
Rdisplay(Rs(:,:,relind), bestassgt);
box on
t = relnames{relind};
title(t);


