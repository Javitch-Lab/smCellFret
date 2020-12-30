function numSegments = getNumSegmentsPerTrack(tracksFinal)
%GETNUMSEGMENTS gives back the number of segments in each compound track
%
%SYNPOSIS numSegments = getNumSegmentsPerTrack(tracksFinal)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman
%OUTPUT numSegments: Vector of number of segments in each compound track
%
%Khuloud Jaqaman, 8 September 2009

%find number of tracks
numTracks = length(tracksFinal);

%initialize numSegments
numSegments = NaN(numTracks,1);

%go over the tracks and get number of segments
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end

