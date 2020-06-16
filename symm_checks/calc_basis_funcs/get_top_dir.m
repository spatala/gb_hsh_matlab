function top_dir = get_top_dir()
curr_pwd = split(pwd,'/');
top_dir = '';
for ct1=1:length(curr_pwd)
    top_dir = strcat(top_dir,curr_pwd{ct1},'/');
    if (strcmp(curr_pwd{ct1},'gb_hsh_matlab'))
        break;
    end
end


