%% Compare between structures
%% Clean up

clearvars
close all
matlabrc
Paths
%% Load both structures

% % load the main structure
main_path = structure_file_path;
% str = load(main_path);
% str = str.str;

[old_root,~] = fileparts(main_path);
old_path = fullfile(old_root,'old ivivstr','str_iviv.mat');
% old_path =  "R:\Share\Simon\Drago_Volker_Simon\Full_data_structure\180522_0811_dataStruct.mat";

old_str = load(old_path);
old_str = old_str.str;
% old_str = old_str.invitro_struct;


% in vitro with old ExcInh, std 3 both windows
% old_path = fullfile(old_root,'Stage1_invitro_only_structure','201130_1826_dataStruct.mat');
% % older in vitro structure, potentially with std 3
% old_path = "R:\Share\Simon\Drago_Volker_Simon\Full_data_structure\200320_0928_dataStruct.mat";
% % older full iviv structure
% old_path = "R:\Share\Simon\Drago_Volker_Simon\_Post_Simon\200219_str_final.mat";
% % older in vitro structure, match
% old_path = "R:\Share\Simon\Drago_Volker_Simon\Full_data_structure\180522_0811_dataStruct.mat";
% % in vitro with new excInh, std 2 both windows
% old_path = fullfile(old_root,'Stage1_invitro_only_structure','201201_1247_dataStruct.mat');
% % in vitro with old excInh, std 2 both windows
% old_path = fullfile(old_root,'Stage1_invitro_only_structure','201201_1545_dataStruct');
% in vitro, old exc inh, std 2 , commented lines
% old_path = fullfile(old_root,'Stage1_invitro_only_structure','201201_1956_dataStruct');
% % in vitro, working!!!
% old_path = fullfile(old_root,'Stage1_invitro_only_structure','201202_1916_dataStruct');

% iviv rebuilt
old_path = structure_file_path;

str = load(old_path);
% str = str.invitro_struct;
str = str.str;
%% Compare fields

disp(setdiff(fields(str), fields(old_str)))
disp(setdiff(fields(old_str), fields(str)))
%% Compare the maps

close all

% get the number of maps
cell_number = min([length(str) length(old_str)]);
% define the number of plots per figure
max_plots = 25;


% for polarity
for polarity = 1:2
    
    figure
    % initialize a counter
    counter = 1;
    
     switch polarity
         case 1
             polarity_term = 'exc';
         case 2
             polarity_term = 'inh';
     end
    % for all the maps
    for maps = 1:cell_number
        
        % generate a new figure if max plots is reached
        if counter > max_plots
            figure
            % reset the counter
            counter = 1;
        end
        
        subplot(round(sqrt(max_plots)),ceil(sqrt(max_plots)),counter)

        % get the maps and concatenate
        map1 = str(maps).(strcat('subpixel_',polarity_term,'Map'));
        map2 = old_str(maps).(strcat('subpixel_',polarity_term,'Map'));
        
%         map1 = str(maps).(strcat(polarity_term,'Map'));
%         map2 = old_str(maps).(strcat(polarity_term,'Map'));

        catmap = horzcat(map1,map2);
        
        imagesc(catmap)
        set(gca,'TickLength',[0 0],'XTick',[],'YTick',[])
        axis equal
        axis tight

        % update the counter
        counter = counter + 1;
    end
    
end
autoArrangeFigures