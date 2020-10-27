% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% September 30, 2020
% Initial Simulation for Mutualism

%% Simulation Parameters %%%
% Seed the random number generator for testing
rng_set = rng(123456789);

% Time-related variables
dt = 1;             % timestep, increment by days
simLength = 500;    % length of simulation: 1 year
numIterations = 1 + simLength/dt;
animation_fps = 100;  % Speed of visualization

% Grid dimensions
row_count = 30; % width
col_count = 30; % length

%% Constants %%
EMPTY = 1;
GROWING_PLANT = 2;
PLANT = 3;
POLLINATED_PLANT = 4;
ANIMAL = 5;
POLLINATED_ANIMAL = 6;
POLLEN = 7;

prob_init_plant = 0.01; % initial probability a cell is plant
prob_init_animal = 0.00; % initial probability a cell is animal
prob_plant_death = 0.05; % probability plant death at each timestep
prob_pollination = 0.5; % probability animal will pick up/drop pollen
prob_pollen_production = 0.5; % chance a plant will produce pollen
prob_pollen_spread = 0.15; % chance pollen spreads to an adjacent empty cell
prob_animal_reproduce = 0.75; % chance 2 animals reproduce if conditions are good

init_animal_count = 10000; % Number of animals at initial animal spawn points
init_plant_count = 1; % Number of animals at initial animal spawn points
init_pollen_count = 100; % Amount of pollen spawned when it is spawned
animal_diffuse_rate = 0.01; % Amount of animals that will diffuse over to lower pop tile
pollen_diffuse_rate = 0.01; % Amount of pollen that will diffuse to lower conc tile


%% Counters for statistics
plant_counter = zeros(1, numIterations); % Keep track of the number of plants
animal_counter = zeros(1, numIterations); % Keep track # of animals
pollen_conc_counter = zeros(1, numIterations); % Keep track # of animals


%% Set up grids
% Initialize grid to be all empty
grids = ones(row_count, col_count, numIterations) * EMPTY;

animal_conc_grids = zeros(row_count, col_count, numIterations) * EMPTY;
pollen_conc_grids = zeros(row_count, col_count, numIterations) * EMPTY;


for row = 1:row_count
    for col = 1:col_count
        % Grid initialization, plants first
        plant_chance = rand;
        animal_chance = rand;
        
        % Get random number, if smaller than p, add a plant
        if plant_chance < prob_init_plant
            grids(row, col, 1) = PLANT;

        % Now adding animals if smaller than p
        elseif animal_chance < prob_init_animal
            % grids(row, col, 1) = ANIMAL;
            % animal_conc_grids(row, col, 1) = init_animal_count;
            grids(row, col, 1) = ANIMAL;
            animal_conc_grids(row,col,1) = init_animal_count;

        end
    end
end

% % First value for counters
% animal_counter(1) = sum(sum(grids(:,:,1)==ANIMAL)) +...
%                               sum(sum(grids(:,:,1)==POLLINATED_ANIMAL));
% plant_counter(1) = sum(sum(grids(:,:,1)==PLANT)) + ...
%                              sum(sum(grids(:,:,1)==GROWING_PLANT)) + ...
%                              sum(sum(grids(:,:,1)==POLLINATED_PLANT));

pollen_conc_counter(1) = sum(sum(animal_conc_grids(:,:,1)));
disp("Grids Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    
    %% Absorbing boundary condition
    % Create grids that are the size of the grid + 2 on each side
    extended_grid_size = size(grids( : , : , frame-1))+2;
    % Initialize all of the cells as empty
    extended_grid = ones(extended_grid_size) * EMPTY; 
    extended_animal_grid = zeros(extended_grid_size) * EMPTY; 
    extended_pollen_grid = zeros(extended_grid_size) * EMPTY; 


    % Set the inside portion of the grids equal to the corresponding values 
    % from the previous timestep (a.k.a the previous frame)
    extended_grid(2:end-1, 2:end-1) = grids(:,:,frame-1);
    extended_animal_grid(2:end-1, 2:end-1) = animal_conc_grids(:,:,frame-1);
    extended_pollen_grid(2:end-1, 2:end-1) = pollen_conc_grids(:,:,frame-1);

    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to original(non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1

            % Store current cell
            current_cell = extended_grid(row, col);
            current_animal_conc = extended_animal_grid(row,col);
            current_pollen_conc = extended_pollen_grid(row,col);

            % Updated Cell Values
            updated_cell = current_cell;
            updated_animal_conc = current_animal_conc;
            updated_pollen_conc = current_pollen_conc;


            % Getting Moore Neighborhood
            north = extended_grid(row - 1, col);
            east  = extended_grid(row, col - 1);
            south = extended_grid(row + 1, col);
            west  = extended_grid(row, col + 1);
            
            northeast = extended_grid(row - 1, col - 1);
            southeast  = extended_grid(row + 1, col - 1);
            northwest = extended_grid(row - 1, col + 1);
            southwest  = extended_grid(row + 1, col + 1);

            % Put all of the neighbors into a list
            neighbors = [north, east, south, west, ...
             northeast, southeast, northwest, southwest];

            % Same thing for animal population neighbors
            % Getting Moore Neighborhood
            north_pollen = extended_pollen_grid(row - 1, col);
            east_pollen  = extended_pollen_grid(row, col - 1);
            south_pollen = extended_pollen_grid(row + 1, col);
            west_pollen  = extended_pollen_grid(row, col + 1);
            
            northeast_pollen = extended_pollen_grid(row - 1, col - 1);
            southeast_pollen  = extended_pollen_grid(row + 1, col - 1);
            northwest_pollen = extended_pollen_grid(row - 1, col + 1);
            southwest_pollen  = extended_pollen_grid(row + 1, col + 1);

            % Put all of the neighbors into a list
            neighbors_pollen = [north_pollen, east_pollen, south_pollen, west_pollen, ...
             northeast_pollen, southeast_pollen, northwest_pollen, southwest_pollen];

            % List of Neighbor coordinates
            neighbor_coords = [ row-1 col; row col-1; row+1 col; ...
                row col+1; row-1 col-1; row+1 col-1; row-1 col+1; row+1 col+1 ];                    


            %% Update cell

            % Getting different counts of neighbors
            empty_count = sum(neighbors == EMPTY);
            pollen_count = sum(neighbors == POLLEN);

            norm_plant_count = sum(neighbors == PLANT);
            poll_plant_count = sum(neighbors == POLLINATED_PLANT);
            norm_animal_count = sum(neighbors == ANIMAL);
            poll_animal_count = sum(neighbors == POLLINATED_ANIMAL);
    
            % Total animal count
            animal_count = norm_animal_count + poll_animal_count;
            % Total plant count
            plant_count = poll_plant_count + norm_plant_count;

            % Getting total population of all animals in neighboring cells
            neighbor_animal_pop = 0;
            if(animal_count > 0)
                for val = 1:length(neighbors)
                    neighbor_animal_pop = neighbor_animal_pop + extended_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2));                 
                    neighbor_animal_pop = round(neighbor_animal_pop);
                end
            end

            % Getting total number of all pollen in neighboring cells
            neighbor_pollen_conc = 0;
            if (pollen_count > 0)
                for val = 1:length(neighbors)
                    neighbor_pollen_conc = neighbor_pollen_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    neighbor_pollen_conc = round(neighbor_pollen_conc);
                end
            end


            % Empty cell behavior
            if(current_cell == EMPTY)
                % Plant is producing pollen on this cell
                if(plant_count > 0)
                    % All neighboring plants have a chance to produce pollen
                    for val = 1:plant_count
                        if(rand < prob_pollen_production)
                            updated_cell = POLLEN;
                            updated_pollen_conc = updated_pollen_conc + init_pollen_count;
                        end
                    end
                elseif(pollen_count > 0)
                    % Setting pollen concentration based on neighbors
                    % Using sum of all neighbor values so no * pollen_count needed
                    % Becomes pollen if it will have more than 1 pollen within it after timestep, else stays empty
                    new_pollen_count = neighbor_pollen_conc * pollen_diffuse_rate;
                    if(new_pollen_count >= 1)
                        updated_cell = POLLEN;
                        updated_pollen_conc = new_pollen_count;
                    else
                        updated_cell = EMPTY;
                        updated_pollen_conc = 0;
                    end
                else
                    updated_cell = EMPTY;
                    updated_pollen_conc= 0;
                end
            end


            % Pollen Cell Behavior
            if(current_cell == POLLEN)
                % If neighboring no empty cells or other pollen cells, stay same
                if(empty_count <= 0 && pollen_count <= 0)
                    updated_cell = POLLEN;
                end

                % Plants produce pollen and add it to this cells current pollen
                concentration_gain1 = 0;
                if(plant_count > 0)
                    % All neighboring plants have a chance to produce pollen
                    for val = 1:plant_count
                        if(rand < prob_pollen_production)
                            updated_cell = POLLEN;
                            concentration_gain1 = concentration_gain1 + init_pollen_count;
                        end
                    end
                end

                % If there are empty cells than concentration will be decreased by
                % the empty cell count times the diffuse rate since each empty cell
                % will take current concentration * diffuse rate out
                concentration_loss1 = 0;
                if(empty_count > 0)
                    concentration_loss1 = current_pollen_conc * ((pollen_diffuse_rate) * empty_count);
                    % If losing more concentration than is possible, just lose current amount
                    % Similar for negative values
                    if(concentration_loss1 > current_pollen_conc)
                        concentration_loss1 = current_pollen_conc;
                    elseif(concentration_loss1 < 0)
                        concentration_loss1 = 0;
                    end
                end

                % If there are pollen neighbors then some pollen will diffuse to
                % lower population neighbors. So either gain or lose concentration
                % based on neighbor populations
                concentration_loss2 = 0;
                concentration_gain2 = 0;
                if(pollen_count > 0)
                    for val = 1:length(neighbors)
                        % Only do this stuff to pollen neighbors
                        if(neighbors(val) == POLLEN) 
                            cur_neighbor_conc = neighbors_pollen(val);
                            % Adding to this cell if neighbor has more pollen
                            if(cur_neighbor_conc > current_pollen_conc)
                                % Gaining concentration, so gain is increased
                                concentration_gain2 = concentration_gain2 + (cur_neighbor_conc * pollen_diffuse_rate);
                            % If current population is greater than neighbor than lose some animals
                            elseif(cur_neighbor_conc < current_pollen_conc)
                                % Losing concentration so concentration loss is increased
                                concentration_loss2 = concentration_loss2 + (current_pollen_conc * pollen_diffuse_rate);
                            end
                        end
                    end
                    if(concentration_loss2 > current_pollen_conc)
                        concentration_loss2 = current_pollen_conc;
                    elseif(concentration_loss2 < 0)
                        concentration_loss2 = 0;
                    end
                end


                % Setting the new populations based on population loss / gain
                updated_pollen_conc = (current_pollen_conc - concentration_loss1 - concentration_loss2 + concentration_gain1 + concentration_gain2);

                % Becomes empty if no pollen left, otherwise empty
                if(updated_pollen_conc >= 1)
                    updated_cell = POLLEN;
                else
                    updated_cell = EMPTY;
                    updated_pollen_conc = 0;
                end
            end

            % CONCENTRATION GRIDS SHOULD NOT BE UPDATED DURING LOOP BUT CONSTANTLY POTENTIALLY MAYBE
            % OR MAKE AN CONC CHANGE GRID, temp grid that stores changes to the concentrations and then just adds it to the normal grid
            % Is this even needed or does that happen implicitly with the structure we already have?
            % Animal Stuff
            % If animal next to pollen
                % Find the max pollen neighbor
                % X animals will move to that neighbor based on animal_diffuse_rate (with pollen)









            % Plant Behavior
            if(current_cell == PLANT)
                updated_cell = PLANT;
            end

            
            % Updating next grid with the new cell value
            grids(row-1, col-1 , frame) = updated_cell;
            pollen_conc_grids(row-1, col-1, frame) = updated_pollen_conc;
        end
    end
    
    % % Recalculate the number of animals and plants for tracking
    % animal_counter(frame) = sum(sum(grids(:,:,frame)==ANIMAL)) +...
    %                           sum(sum(grids(:,:,frame)==POLLINATED_ANIMAL));
    % plant_counter(frame) = sum(sum(grids(:,:,frame)==PLANT)) + ...
    %                          sum(sum(grids(:,:,frame)==GROWING_PLANT)) + ...
    %                          sum(sum(grids(:,:,frame)==POLLINATED_PLANT));

    pollen_conc_counter(frame) = sum(sum(pollen_conc_grids(:,:,frame)));
end

disp("All grids calculated");

%% Visualize the grid

% Create the window for the animation
viz_fig = figure;
viz_axes = axes(viz_fig);
viz_fig2 = figure;

% % Set the colors
% map = [ 1       1       1;          % Empty Cell: white
%         0  1  0;                    % Growing Plant Cell
%         109/255 188/255 0;          % Plant Cell
%         109/255 120/255 0;          % Pollinated plant
%         125/255 100/255  70/255;          % Animal Cell: grey blue 
%         0       0      0;          % Pollinated Animal Cell: light blue
%         230/255 255/255 0];         % Pollen color        

% colormap(viz_axes, map); 

% % Remove axis labels, make aspect ratio look good, and maintain that state
% axis off;
% axis equal;
% hold on;

disp("Drawing...");
for i = 1:numIterations

    
    % Uncomment following line to allow for frame-by-frame viewing of grid
    %w = waitforbuttonpress;

    % if(i>1) % Refresh the image
    %     delete(time_counter_text);
    %     delete(animal_counter_text);
    %     delete(plant_counter_text);
    % end
    
    % Turn each grid into an image
    heatmap(viz_fig2, grids(:,:,i));
    hmap = heatmap(viz_fig, pollen_conc_grids(:, :, i));
    caxis(hmap, [0 30]);
    
    % Draw text
    % time_counter_text = text(viz_axes,0.1,1.025,"Current Time Step: " + ...
    %     string(i), 'Units', 'Normalized');
    % animal_counter_text = text(viz_axes,0.9,0,"Animals: " + ...
    %         animal_counter(i), 'Units', 'Normalized');
    % plant_counter_text = text(viz_axes,0.9,0.05,"Plants: " + ...
    %     plant_counter(i), 'Units', 'Normalized');
    
    pause(1/animation_fps);
end

% Create the graphs for the counters
animal_fig = figure;
animal_axes = axes(animal_fig);

% plant_fig = figure;
% plant_axes = axes(plant_fig);

% both_fig = figure;
% both_axes = axes(both_fig);

plot(animal_axes, 1:numIterations, pollen_conc_counter, 'color', [1 0.5 0.17]);
% plot(plant_axes, 1:numIterations, plant_counter);
% plot(both_axes, 1:numIterations, plant_counter);
% hold on;
% plot(both_axes, 1:numIterations, animal_counter);
% hold off;


% title(animal_axes, "Animal Population Throughout Simulation");
% ylabel(animal_axes, "Number of animals (incl. pollen-carrying)");
% xlabel(animal_axes, "Time step");
% ylim(animal_axes, [0, max(animal_counter)+5]); % Extend the maximum y limit

% title(plant_axes, "Plant Population Throughout Simulation");
% ylabel(plant_axes, "Number of plants (incl. pollen-carrying)");
% xlabel(plant_axes, "Time step");
% ylim(plant_axes, [0, max(plant_counter)+5]); % Extend the maximum y limit

% title(both_axes, "Population of Plants and Animals");
% ylabel(both_axes, "Number of Plants and Animals");
% xlabel(both_axes, "Time step");
% legend('Plants', 'Animals');
% ylim(both_axes, [0, max(plant_counter)+5]); % Extend the maximum y limit


disp("Simulation complete!");