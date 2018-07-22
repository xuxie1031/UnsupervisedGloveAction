%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script of dumping hand data into feature vector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function feature_dump_vec(data_dir, data_name, dim)
    data_list=dir([data_dir, data_name]);
    if(~exist('hand_feature', 'dir'))
        mkdir('hand_feature');
    end

    for i=1:length(data_list)
        filename=[data_dir,data_list(i).name];
        fid=fopen(filename,'r');
        data=textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        
        fprintf('dumping %s, waiting...\n',data_list(i).name);
        wrist_vec=zeros(length(data{:}),dim);
        for j=1:length(data{:})
            str=cell2mat(data{:}(j));
            strs=regexp(str,',','split');
            wrist_vec(j,:)=str2double(strs(1:(end-1)));
        end
        
        % preset force magnitude upperbound and normalization
        fmax=1000.0;
        wrist_vec=wrist_vec/fmax;
        wrist_vec=single(wrist_vec);
        
        save_name=['hand_feature/',data_list(i).name(1:end-4),'.mat'];
        save(save_name,'wrist_vec');
    end
end