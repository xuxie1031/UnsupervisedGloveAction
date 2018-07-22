%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve the proceeding segment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segment=next_segment(X,labels,idx_b)
	sub_labels=labels(idx_b:end);
	[M, ~]=regexp(sprintf('%i',[0 diff(sub_labels')==0]), '1+', 'match');
	stride=length(cell2mat(M(1)));
	idx_e=idx_b+stride;
	segment=X(idx_b:idx_e,:);
end

