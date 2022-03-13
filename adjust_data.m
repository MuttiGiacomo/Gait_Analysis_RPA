%% Data reading and format adjustment

% Load the data from a .txt dump file.
% Read line by line the samples, assign the corresponding round (number of 
% revolutions) and write the new line in the output file.
% The information on the round is used only offline for managing and
% plotting the data better, as well as a guideline during labelling.
% In a real time context, the information on the round is unnecessary since
% all the operations needed to extract features and detect the phases work
% on the point cloud of a single turn (+ memory of the features of the
% previous one) -> see other files

clc;
clearvars;

% data_lab: folder with original dump data of measurements taken in lab
% data_lab_adj: output folder with adjusted files

% data_lab/test_lab_1.txt
% data_lab/test_lab_2.txt
% data_lab/test_lab_3.txt
% data_lab/test_lab_4.txt
% data_lab/test_lab_5.txt
% data_lab/test_lab_6.txt
% data_lab/test_lab_7.txt
% data_lab/test_lab_8.txt
% data_lab/test_lab_9.txt
% data_lab/test_lab_10.txt

%% Load data

fidIN = fopen('data_lab/test_lab_5.txt', 'r');
fidOUT = fopen('data_lab_adj/test_5_forward.txt', 'wt');

% The original code for the dump write an 's' in front of the line
% that begins each new round. We kept this convention and so we
% detect each round by looking at the first character at the beginning of
% each line and comparing it to 's'

r = 0;
new_round = false;
while 1
  thisline = fgetl(fidIN);
  if ~ischar(thisline); break; end
    s1 = thisline(1);
    if strcmp(s1,'s')
          r = r+1;
    end
  thisline = sprintf('%s %s \n',num2str(r),thisline(2:end));
  fwrite(fidOUT, thisline, 'uchar');
end
r_max = r;
fclose(fidIN);
fclose(fidOUT);
