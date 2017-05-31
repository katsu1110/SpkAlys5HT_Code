function matching_idx = findIdxMatchingID(fnames_sorting, id_ex, monkey)
% find the index that matches the unit's ID in exinfo with the ones from the
% spike sorting file
% 
%
% fname_spksort is the cell array containing the file names from the
% sorting table
% 
% id is the unit's id from the corresponding ex file saved in exinfo.fname /  exinfo.fname_drug 
% 
% monkey is a string specifying 'ma' (mango) or 'ka' (kaki)
%
%
% @CL February 2, 2017


%in case of kaki's id, I added 0.5 to differentiate between monkeys
% this is redundant given the monkey information and has to be undone
id_ex = floor(id_ex); 

s = getStringID(id_ex, monkey); % convert the id to a string with 0 going in front of it
id_ex_s = [monkey '_' s];


id_prefix_sorting = getIDX(fnames_sorting);
matching_idx = find( ~cellfun(@isempty, strfind(id_prefix_sorting, id_ex_s)) );

end

%%
function s = getIDX(fnames)
% returns the monkey and ID of each row prefix
% this eventually allows to read out all entries of the same unit

for i = 1:length(fnames)

    s{i} = fnames{i}(1:7);  % monkey and session number

end


end

