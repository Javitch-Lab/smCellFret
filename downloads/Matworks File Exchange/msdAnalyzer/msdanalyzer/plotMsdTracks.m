%% Plot msd tracks
figure;
hold(gca, 'on');

%% Plot Donor Tracks
dataCh2    = load('diffTracksCh2.mat');
nTracksCh2 = length(dataCh2.diffTracks);
indicesCh2 = 1:nTracksCh2;
hpsCh2     = NaN(nTracksCh2, 1);
 
for i = 1 : nTracksCh2
   
    index = indicesCh2(i);
    track = dataCh2.diffTracks{index};
    x = track(:,2);
    y = track(:,3);
    hpsCh2(i) =  plot(gca, x, y,'g');

end

%% Plot Acceptor Tracks
dataCh1 = load('diffTracksCh1.mat');
nTracksCh1 = length(dataCh1.diffTracks);
indicesCh1 = 1:nTracksCh1;
hpsCh1     = NaN(nTracksCh1, 1);
 
for i = 2 : 2 %nTracksCh1
   
    index = indicesCh1(i);
    track = dataCh1.diffTracks{index};
    x = track(:,2);
    y = track(:,3);
    hpsCh1(i) =  plot(gca, x, y,'r');

end

