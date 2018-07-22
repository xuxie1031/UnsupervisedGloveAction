%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate ground truth gaussian parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gaussian_fit(data_dir, data_name, dim, nlabel)
    means=zeros(nlabel,dim);
    covs=zeros(nlabel,dim^2);

    for l=0:(nlabel-1)
        data_list=dir([data_dir,data_name,'*ground_label',num2str(l),'.csv']);    % sample from multiple trials
        feature=[];
        fprintf('dumping l=%d...\n',l);
        for i=1:length(data_list)
            filename=[data_dir,data_list(i).name];
            fid=fopen(filename,'r');
            data=textscan(fid,'%s','delimiter','\n');
            fclose(fid);

            wrist_vec=zeros(length(data{:}),dim);
            for j=1:length(data{:})
                str=cell2mat(data{:}(j));
                strs=regexp(str,',','split');
                wrist_vec(j,:)=str2double(strs(1:(end-1)));
            end

            fmax=1000.0;
            wrist_vec=wrist_vec/fmax;
            wrist_vec=single(wrist_vec);
            feature=[feature;wrist_vec];
        end
        
        means(l+1,:)=mean(feature);
        covf=cov(feature);                      
        covf=covf+eye(size(covf))*(-2*min(eig(covf)));
        covf=reshape(covf,1,dim^2);
        covs(l+1,:)=covf;
    end

    if(~exist('gaussian_params', 'dir'))
        mkdir('gaussian_params');
    end
    save_name=['gaussian_params/',data_name,'_gaussian_params.mat'];
    save(save_name,'means','covs');
end