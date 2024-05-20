%script to plot the best objective against the number of function calls

dirnum = 5;
dirstr = num2str(dirnum);
bestJhist = fullfile('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4/inflam+minrad+radcap', dirstr, 'results', 'BestJhist.dat');
log = fullfile('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4/inflam+minrad+radcap', dirstr, 'log.txt');

%READING JHIST
% Initialize an empty array to store the extracted numbers
jhist_numbers = zeros(1,countlines(bestJhist));
% 
% % Open and read the file
% fid = fopen(bestJhist, 'r');
% while ~feof(fid)
%     line = fgetl(fid);
%     if ischar(line)
%         % Split the line by whitespace and convert the first element to a number
%         number = str2double(strsplit(line, '\t'));
%         jhist_numbers = [jhist_numbers, number(1)];
%     end
% end
% fclose(fid);
% 
% 
% %READING LOG
% % Initialize an empty array to store the extracted numbers
% numbers = [];
% 
% % Open and read the file
% fid = fopen('log.txt', 'r');
% while ~feof(fid)
%     line = fgetl(fid);
%     if ischar(line)
%         % Define the regular expression pattern
%         pattern = 'Number of function evaluations till now : (\d+)';
%         
%         % Try to match the pattern in the line
%         match = regexp(line, pattern, 'tokens');
%         
%         if ~isempty(match)
%             % Extract the number and convert it to a numeric value
%             number = str2double(match{1}{1});
%             numbers = [numbers, number];
%         end
%     end
% end
% fclose(fid);
% 
% 
% %PLOTTING
