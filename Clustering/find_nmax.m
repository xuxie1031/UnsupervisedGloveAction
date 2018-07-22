%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the length of the longest segment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nmax = find_nmax(labels)
    [M, ~]=regexp(sprintf('%i', [0 diff(labels')==0]), '1+', 'match');
    [nmax, ~]=max(cellfun('length',M));
    nmax=nmax+1;
end