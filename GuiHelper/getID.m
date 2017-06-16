function id = getID( expInfo )


[id, idxsort] = sort([expInfo.id]');
is5HT = [expInfo(idxsort).is5HT]';
isc2 = [expInfo(idxsort).isc2]';

[id is5HT isc2]


end

