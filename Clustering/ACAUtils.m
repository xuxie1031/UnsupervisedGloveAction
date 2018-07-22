%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACA utility: pre-compute variables/parameters for ACA procedure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ACAUtils(file_names, NCmin, NCmax)
    if(~exist('hand_feature', 'dir') || ~exist('hcluster', 'dir'))
        error('hierarchical clustering result not prepared, exit...');
    end
    
    save_dir='../ACA/ACAbin/';
    if(~exist(save_dir, 'dir'))
        mkdir(save_dir);
    end

    % iterate through different tasks
    for file_id=1:length(file_names)  
        for NC=NCmin:NCmax
            fprintf('ACA dump file name: %s\nNC: %d\n',file_names{file_id},NC);
            data_name=['hand_feature/',file_names{file_id},'.mat'];
            Cids_name=['hcluster/C_',num2str(NC),'_',file_names{file_id},'.mat'];
            data=load(data_name);
            features=data.wrist_vec;
            result=load(Cids_name);
            Cids=result.Cids;
            
            [n,~]=size(features);
%             nmax=find_nmax(Cids);
            k=length(unique(Cids));
            
            save_name_Cids=[save_dir,'Cids_C_',num2str(NC),'_',file_names{file_id}];
            save_name_Kb=[save_dir,'Kb_C_',num2str(NC),'_',file_names{file_id}];
            name_Kb=['hcluster/Kb_C_',file_names{file_id},'.mat'];
            if exist(name_Kb,'file')
                load(name_Kb);
            else
                sigma_=0.8;
                threshold_=1.0;         % params for caluculating kernel
                Kb=kernel_binary(features,features,sigma_,threshold_);
                save(name_Kb,'Kb');
            end
            
%             save_name_Zcs=['../binaries/Zcs_C_',num2str(NC),'_',file_names{file_id}];
%             Zcs=zeros(k,dim);
%             for i=1:k
%                 feature_c=features(find(Cids==i),:);
%                 Zcs(i,:)=mean(feature_c);
%             end
            
            save_name_mcs=[save_dir,'mcs_C_',num2str(NC),'_',file_names{file_id}];
            mcs=zeros(1,k);
            for c=1:k
                idy=1;
                while idy<n
                    segy=next_segment(features,Cids,idy);
                    [len_segy,~]=size(segy);
                    if Cids(idy)==c
                        mcs(c)=mcs(c)+1;
                    end
                    idy=idy+len_segy;
                end
            end
            
            save_name_thau=[save_dir,'thauYYs_C_',num2str(NC),'_',file_names{file_id}];
            thau_YYs=zeros(1,k);
            for c=1:k
                idy1=1;
                while idy1<n
                    segy1=next_segment(features,Cids,idy1);
                    [len_segy1,~]=size(segy1);
                    if Cids(idy1)==c
                        idy2=1;
                        while idy2<n
                            segy2=next_segment(features,Cids,idy2);
                            [len_segy2,~]=size(segy2);
                            if Cids(idy2)==c
                                Gamma=DTAK(Kb(idy1:(idy1+len_segy1-1),idy2:(idy2+len_segy2-1)));
                                thau_YYs(c)=thau_YYs(c)+Gamma;
                            end
                            idy2=idy2+len_segy2;
                        end
                    end
                    idy1=idy1+len_segy1;
                end
                thau_YYs(c)=thau_YYs(c)/mcs(c)^2;
            end
            thau_YYs=thau_YYs/sum(thau_YYs);
            
            save_name_segpos=[save_dir,'segpos_C_',num2str(NC),'_',file_names{file_id}];
            seg_pos=[];
            idx=1;
            count=0;
            while idx<n
                count=count+1;
                segx=next_segment(features,Cids,idx);
                [len_segx,~]=size(segx);
                seg_pos(count,1)=idx;
                seg_pos(count,2)=idx+len_segx-1;
                idx=idx+len_segx;
            end
            
            fprintf('saving clustering segmentations as index poses...\n');
            fid=fopen(save_name_segpos,'wb');
            [len_seg_pos,~]=size(seg_pos);
            for i=1:len_seg_pos
                fwrite(fid,int32(seg_pos(i,:)),'int32');
            end
            fclose(fid);
            
            fprintf('saving clustering kernels...\n');
            fid=fopen(save_name_Kb,'wb');
            for i=1:n
                fwrite(fid,Kb(i,:),'double');
            end
            fclose(fid);
            
            fprintf('saving clustering ids...\n');
            fid=fopen(save_name_Cids,'wb');
            fwrite(fid,int32(Cids),'int32');
            fclose(fid);
            
            fprintf('saving ACA MCs...\n');
            fid=fopen(save_name_mcs,'wb');
            fwrite(fid,int32(mcs),'int32');
            fclose(fid);
            
            fprintf('saving ACA thauYYs...\n');
            fid=fopen(save_name_thau,'wb');
            fwrite(fid,thau_YYs,'double');
            fclose(fid);
            
            fprintf('saving other params for ACA...\n')
            save_name_params=[save_dir,'params_C_',num2str(NC),'_',file_names{file_id}];
            fid=fopen(save_name_params,'wb');
            fwrite(fid,int32(n),'int32');
            fwrite(fid,int32(k),'int32');
            fwrite(fid,int32(len_seg_pos),'int32');
            fclose(fid);
        end
    end
end