%% STORDATA
% stores important/selected data in destination folder
% save per session/run

folder_name = 'Stored_Sim_Data';
folder_name2 = fullfile(pwd,folder_name,Case_name);

[status, msg] = mkdir(folder_name);
save(folder_name2,'Case_name','INP','OPT','SET','FLD','TCC','X','FRS','EDG','GRD','HVR','SOL','FLP','BLC','MON','PLT');
