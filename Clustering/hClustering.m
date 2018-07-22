%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hierarchcical clustering on hand feature %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hClustering(minC, maxC)
    if(~exist('hand_feature', 'dir'))
        error('hand feature folder does not exist, exit...');
    end

    data_list=dir('hand_feature');
    for i=1:length(data_list)
        [~,~,ext]=fileparts(data_list(i).name);
        if(strcmp(ext,'.mat'))
            filename=['hand_feature/',data_list(i).name];
            load(filename)
        else continue
        end
 
        features=wrist_vec;
        fprintf('hierarchical clustering, waiting...\n')
        %Tree=clusterdata(features,'distance','euclidean','linkage','ward','maxclust',30);
        Z=linkage(features,'ward','euclidean');

        fprintf('generating cluster ids...\n')
        if(~exist('hcluster', 'dir'))
            mkdir('hcluster');
        end
        for k=minC:maxC
            [~,Cids]=dendrogram(Z,k);
            idname=['hcluster/C_',num2str(k),'_',data_list(i).name];
            save(idname,'Cids');
        end
    end
end