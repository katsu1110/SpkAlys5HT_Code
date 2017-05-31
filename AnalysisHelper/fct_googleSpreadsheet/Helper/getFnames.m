function fnames = getFnames(ex)

if isfield(ex.Header, 'Headers')
    fnames = {ex.Header.Headers.fileName};
elseif isfield(ex.Header, 'Sort')
    fnames = {ex.Header.Sort.exFileName};
else
    fnames = ex.Header.fileName;
end


% replace the wrong directory to get to the data via the true, newly formed directory
if contains(fnames{1}, 'ka_')
    monkey = 'kaki';
elseif contains(fnames{1}, 'ma_')
    monkey = 'mango';
end
id_s = findIDhelper(fnames{1}, monkey); % unit ID


% use the ex files on the network in the data folder
% since some files were created from raw files that were stored at a
% different place, the directory mus be changed.
for i = 1:length(fnames)
   fnames{i} = strrep(fnames{i}, 'Lenka\LenkaSort', ['data\' monkey '\' id_s]);
   fnames{i} = strrep(fnames{i}, 'Lenka\Mango', ['data\mango\' id_s]);
   fnames{i} = strrep(fnames{i}, 'D:\Users\hn\Desktop\Mango', 'Z:\data\mango');
   
   if contains(fnames{i}, 'Lenka\Kaki\0242')
       fnames{i} = strrep(fnames{i}, 'Lenka\Kaki\0242', ['data\kaki\' id_s]);
   else
       fnames{i} = strrep(fnames{i}, 'Lenka\Kaki', ['data\kaki\' id_s]);
   end
   

   % Lenka saved this file with a different ID
    if strcmp(id_s, '0272') && strcmp(monkey, 'kaki')
           fnames{i} = strrep(fnames{i}, '0272', '0273');
    end


end

end


function id_s = findIDhelper(fname, monkey)
% searches the file names for the unit ID number

monkey = monkey(1:2);
strt_i  = regexp(fname, [monkey '_0']);
strt_i = strt_i+3;

id_s = fname(strt_i:strt_i+3);

end