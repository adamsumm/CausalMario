cd ..;
currdir = pwd;
cd data;
scorefile    =[currdir, '/output/cont1_status'];
conceptzsfile=[currdir, '/output/cont1_dom0'];

% original data
load('cont1.mat');

conceptzs=load(conceptzsfile);
score    =load(scorefile);

% find configuration with the highest probability 
[m,mind]    =max(score(:,5));
bestassgt   =conceptzs(mind,:);

Rdisplay(G, bestassgt);
box on


