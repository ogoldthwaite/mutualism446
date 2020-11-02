% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% November 1, 2020
% Generic pollination mutualism cellular automata model

%% Simulation Parameters %%%
% Seed the random number generator for testing
rng_set = rng(123456789);

% Time-related variables
dt = 1;             % timestep, increment by days, need to go by hours
simLength = 100;    % length of simulation: 1 year
numIterations = 1 + simLength/dt;
animation_fps = 50;  % Speed of visualization

% Grid dimensions
row_count = 30; % width
col_count = 30; % length

%% Constants %%
EMPTY = 0;
GROWING_PLANT = 2;
PLANT = 3;
POLLINATED_PLANT = 4;
ANIMAL = 5;
POLLINATED_ANIMAL = 6;
POLLEN = 7;

prob_init_plant = 0.01; % initial probability a cell is plant
prob_init_animal = 0.01; % initial probability a cell is animal
prob_animal_death = 0.05; % probability of animal death at each timestep
prob_plant_death = 0.05; % probability plant death at each timestep
prob_plant_growth = 0.1; % probability plant growth at each timestep
prob_pollination = 0.1; % probability animal will pick up/drop pollen
prob_pollen_production = 0.8; % chance a plant will produce pollen
prob_pollen_spread = 0.15; % chance pollen spreads to an adjacent empty cell
prob_animal_reproduce = .1; % chance 2 animals reproduce if conditions are good
prob_random_move = 0.00; % chance an animal will move randomly with no stimuli

init_animal_count = 1000000; % Number of animals at initial animal spawn points
init_plant_count = 100; % Number of plants at initial plant spawn points
init_pollen_count = 100;     % Initial amount of pollen a pollinated plant has
animal_pollen_diffuse_rate = 0.5; % Amount of animals that will diffuse to pollen tile
animal_empty_diffuse_rate = 0.01; % Amount of animals that will diffuse to empty tile
pollen_diffuse_rate = 0.01; % Amount of pollen that will diffuse

%% Counters for statistics
plant_counter = zeros(1, numIterations); % Keep track # of plants
animal_counter = zeros(1, numIterations); % Keep track # of animals
animal_pop_counter = zeros(1, numIterations); % Keep track # of animals

%% Set up grids
% Initialize grid to be all empty
grids = ones(row_count, col_count, numIterations) * EMPTY;
animal_pop_grids = ones(row_count, col_count, numIterations) * EMPTY;
pollen_conc_grids = ones(row_count, col_count, numIterations) * EMPTY;

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
            grids(row, col, 1) = ANIMAL;
            animal_pop_grids(row, col, 1) = init_animal_count;
        end
    end
end

% First value for counters
animal_counter(1) = sum(sum(grids(:,:,1)==ANIMAL)) + ...
                              sum(sum(grids(:,:,1)==POLLINATED_ANIMAL));
plant_counter(1) = sum(sum(grids(:,:,1)==PLANT)) + ...
                             sum(sum(grids(:,:,1)==GROWING_PLANT)) + ...
                             sum(sum(grids(:,:,1)==POLLINATED_PLANT));
animal_pop_counter(1) = sum(sum(animal_pop_grids(:,:,1)));
pollen_conc_counter(1) = sum(sum(pollen_conc_grids(:,:,1)));

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

    % Grids for calculating change from timestep to timestep
    delta_pollen_grid = zeros(extended_grid_size);
    delta_animal_grid = zeros(extended_grid_size);

    % Set the inside portion of the grids equal to the corresponding values 
    % from the previous timestep (a.k.a the previous frame)
    extended_grid(2:end-1, 2:end-1) = grids(:,:,frame-1);
    extended_animal_grid(2:end-1, 2:end-1) = animal_pop_grids(:,:,frame-1);
    extended_pollen_grid(2:end-1, 2:end-1) = pollen_conc_grids(:,:,frame-1);

    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to original(non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1

            % Store current cell
            current_cell = extended_grid(row, col);
            current_animal_pop = extended_animal_grid(row,col);
            current_pollen_conc = extended_pollen_grid(row,col);

            updated_cell = current_cell;
            updated_animal_pop_cell = current_pollen_conc;
            updated_pollen_cell = current_pollen_conc;
            
            % Getting Moore Neighborhood
            self = extended_grid(row, col);
            north = extended_grid(row - 1, col);
            east  = extended_grid(row, col - 1);
            south = extended_grid(row + 1, col);
            west  = extended_grid(row, col + 1);          
            northeast = extended_grid(row - 1, col - 1);
            southeast  = extended_grid(row + 1, col - 1);
            northwest = extended_grid(row - 1, col + 1);
            southwest  = extended_grid(row + 1, col + 1);

            % Getting Moore Neighborhood for pollen neighbors
            self_p = extended_pollen_grid(row, col);
            north_p = extended_pollen_grid(row - 1, col);
            east_p  = extended_pollen_grid(row, col - 1);
            south_p = extended_pollen_grid(row + 1, col);
            west_p  = extended_pollen_grid(row, col + 1);          
            northeast_p = extended_pollen_grid(row - 1, col - 1);
            southeast_p  = extended_pollen_grid(row + 1, col - 1);
            northwest_p = extended_pollen_grid(row - 1, col + 1);
            southwest_p  = extended_pollen_grid(row + 1, col + 1);

            % Put all of the neighbors into a list
            neighbors = [self, north, east, south, west, ...
             northeast, southeast, northwest, southwest];

            % Put all of the pollen neighbors into a list
            pollen_neighbors = [self_p, north_p, east_p, south_p, west_p, ...
            northeast_p, southeast_p, northwest_p, southwest_p];
    
            % List of Neighbor coordinates
            neighbor_coords = [ row col; row-1 col; row col-1; row+1 col; ...
                row col+1; row-1 col-1; row+1 col-1; row-1 col+1; row+1 col+1 ];   
            
            %% Update cell
            % Getting different counts of neighbors
            empty_count = sum(neighbors == EMPTY);
            pollen_count = sum(neighbors == POLLEN);
            norm_plant_count = sum(neighbors == PLANT);
            poll_plant_count = sum(neighbors == POLLINATED_PLANT);
            grow_plant_count = sum(neighbors == GROWING_PLANT);
            norm_animal_count = sum(neighbors == ANIMAL);
            poll_animal_count = sum(neighbors == POLLINATED_ANIMAL);
    
            % Total animal count
            animal_count = norm_animal_count + poll_animal_count;
            % Total plant count
            % Note: pollen and growing plants not included 
            plant_count = poll_plant_count + norm_plant_count;        
            % Total amount of neighbors
            num_neighbors = animal_count + plant_count + pollen_count;

            % Getting total population of all animals in neighboring cells
            if(animal_count > 0)
                neighbor_animal_pop = 0;
                for val = 2:length(neighbors==ANIMAL)
                    neighbor_animal_pop = neighbor_animal_pop + extended_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2));                 
                    neighbor_animal_pop = round(neighbor_animal_pop);
                end
            end
            
            % Getting total amount of pollen in neighboring plant cells
            if (plant_count > 0)
                plant_n_pollen_conc = 0;
                for val = 2:length(neighbors==PLANT)
                    plant_n_pollen_conc = plant_n_pollen_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    plant_n_pollen_conc = round(plant_n_pollen_conc);
                end
            end
            
            % Getting total amt of pollen in neighboring pollinated animals
            if (poll_animal_count > 0)
                neighbor_poll_animal_conc = 0;
                for val = 2:length(neighbors==POLLINATED_ANIMAL)
                    neighbor_poll_animal_conc = neighbor_poll_animal_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    neighbor_poll_animal_conc = round(neighbor_poll_animal_conc);
                end
            end        
            
            % Getting total number of all pollen in neighboring cells
            neighbor_pollen_conc = 0;
            for val = 2:length(neighbors)
                neighbor_pollen_conc = neighbor_pollen_conc + extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                neighbor_pollen_conc = round(neighbor_pollen_conc);
            end
            
            % Pollen Diffusion Behavior
            if(current_pollen_conc > 0 || plant_count > 0)
                % Plants produce pollen and add it to this cells current pollen
                concentration_gain1 = 0;
                if(plant_count > 0)
                    % All neighboring plants have a chance to produce pollen
                    for val = 1:plant_count
                        if(rand < prob_pollen_production)
                            concentration_gain1 = concentration_gain1 + init_pollen_count;
                        end
                    end
                end

                % Lose pollen to neighboring cells, pollen. 
                % Give the pollen to the cells, any gain will be given back when it is the other cells turn
                concentration_loss2 = 0;
                for val = 2:length(neighbors)
                % Adding and subtracting the pollen for this current location
                % Only do this stuff to pollen neighbors
                    cur_neighbor_conc = extended_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2));
                    % Adding to this cell if neighbor has more pollen
                    if(cur_neighbor_conc < current_pollen_conc)
                        % Adding pollen to neighbors :)
                        delta_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2)) = delta_pollen_grid(neighbor_coords(val,1), neighbor_coords(val,2)) + (current_pollen_conc * pollen_diffuse_rate);
                        % Losing concentration so concentration loss is increased, pollen going to neighbors
                        concentration_loss2 = concentration_loss2 + (current_pollen_conc * pollen_diffuse_rate);
                    end
                end
                % Setting the new populations based on population loss / gain
                delta_pollen_grid(row,col) = delta_pollen_grid(row,col) + (concentration_gain1 - concentration_loss2);
            end
            
            
            % Animal Diffusion Behavior
            if(current_animal_pop > 0)
                
                % Finding max pollen spot
                [max_pollen_conc, max_index] = max(pollen_neighbors);
                
                for val = 1:length(pollen_neighbors)
                    if(neighbor_pollen_conc > 0 && val == max_index)
                        % Moving big chunk O' animals to the max pollen cell
                        delta_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2)) = delta_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2)) + (current_animal_pop * animal_pollen_diffuse_rate);
                        % Losing animals that are moving
                        delta_animal_grid(row, col) = delta_animal_grid(row, col) - (current_animal_pop * animal_pollen_diffuse_rate);
                    else % Randomly move rest of animals to other spots, if they want too
                        if(rand < prob_random_move)
                            % Moving animals to neighboring spots
                            delta_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2)) = delta_animal_grid(neighbor_coords(val,1), neighbor_coords(val,2)) + (current_animal_pop * animal_empty_diffuse_rate);
                            % Losing animals that are moving
                            delta_animal_grid(row, col) = delta_animal_grid(row, col) - (current_animal_pop * animal_empty_diffuse_rate);
                        end
                    end
                    
                end
            end
            
            updated_cell = current_cell;
            
            % Updating next grid with the new cell value
            grids(row-1, col-1 , frame) = updated_cell;
            % animal_pop_grids(row-1, col-1, frame) = updated_animal_pop_cell;
            %pollen_conc_grids(row-2, col-2, frame) = updated_pollen_cell;
        end
        pollen_conc_grids(:,:,frame) = pollen_conc_grids(:,:,frame-1) + delta_pollen_grid(2:end-1, 2:end-1);
        animal_pop_grids(:,:,frame) = animal_pop_grids(:,:,frame-1) + delta_animal_grid(2:end-1, 2:end-1);
    end
    
    % Recalculate the number of animals and plants for tracking
    animal_counter(frame) = sum(sum(grids(:,:,frame)==ANIMAL)) +...
                              sum(sum(grids(:,:,frame)==POLLINATED_ANIMAL));
    plant_counter(frame) = sum(sum(grids(:,:,frame)==PLANT)) + ...
                             sum(sum(grids(:,:,frame)==GROWING_PLANT)) + ...
                             sum(sum(grids(:,:,frame)==POLLINATED_PLANT));
    animal_pop_counter(frame) = sum(sum(animal_pop_grids(:,:,frame)));
    pollen_conc_counter(frame) = sum(sum(pollen_conc_grids(:,:,frame)));
end

disp("All grids calculated");

%% Visualize the grid

% Create the window for the animation
states_fig = figure;
animal_pop_fig = figure;
pollen_fig = figure;

states_axes = axes(states_fig);

% Set the colors
map = [ 1       1       1;          % Empty Cell: white
        0  1  0;                    % Growing Plant Cell: Lime Green
        109/255 188/255 0;          % Plant Cell: Green
        109/255 120/255 0;          % Pollinated plant: Camo Green
        125/255 100/255  70/255;    % Animal Cell: Brown 
        0       0      0;           % Pollinated Animal Cell: Black
        230/255 255/255 0];         % Pollen color: Yellow     

colormap(states_axes, map); 

% Remove axis labels, make aspect ratio look good, and maintain that state
% axis off;
% axis equal;
% hold on;

% animal_pop_fig = figure;

disp("Drawing...");
for i = 1:numIterations

    
    % Uncomment following line to allow for frame-by-frame viewing of grid
    %   w = waitforbuttonpress;
    
    states_heat = heatmap(states_fig, grids(:,:,i));
    animal_pop_heat=heatmap(animal_pop_fig, animal_pop_grids(:,:,i), 'Colormap', copper);
    pollen_heat=heatmap(pollen_fig, pollen_conc_grids(:,:,i));
    
    caxis(animal_pop_heat, [0 30]);
    caxis(pollen_heat, [0 1000]);

    % Turn each grid into an image
%     heatmap(animal_pop_fig, animal_pop_grids(:,:,i));
%     hmap = heatmap(states_fig, animal_pop_grids(:, :, i));
%     caxis(hmap, [0 30]);
    
    if(i>1) % Refresh the image
        % delete(time_counter_text);
        % delete(animal_counter_text);
        % delete(plant_counter_text);
    end
    
    % Turn each grid into an image
    %image(states_axes, grids(:, :, i));
    
%     % Draw text
    % time_counter_text = text(states_axes,0.1,1.025,"Current Time Step: " + ...
    %     string(i), 'Units', 'Normalized');
    % animal_counter_text = text(states_axes,0.9,0,"Animals: " + ...
    %         animal_counter(i), 'Units', 'Normalized');
    % plant_counter_text = text(states_axes,0.9,0.05,"Plants: " + ...
    %     plant_counter(i), 'Units', 'Normalized');
    
    pause(1/animation_fps);
end

% Slows down simulation
% s1 = plot(nan,'s', 'color', map(1,:), 'MarkerFaceColor', map(1,:));
% s2 = plot(nan, 's', 'color', map(2,:), 'MarkerFaceColor', map(2,:));
% s3 = plot(nan, 's', 'color', map(3,:), 'MarkerFaceColor', map(3,:));
% s4 = plot(nan, 's', 'color', map(4,:), 'MarkerFaceColor', map(4,:));
% s5 = plot(nan, 's', 'color', map(5,:), 'MarkerFaceColor', map(5,:));
% s6 = plot(nan, 's', 'color', map(6,:), 'MarkerFaceColor', map(6,:));
% s7 = plot(nan, 's', 'color', map(7,:), 'MarkerFaceColor', map(7,:));
% legend([s1, s2, s3, s4, s5, s6, s7], {'Empty', 'Growing Plant', 'Plant'...
%     'Pollinated Plant', 'Animal', 'Pollinated Animal', 'Pollen'}, ...
%     'Location', 'bestoutside');

% Create the graphs for the counters
animal_fig = figure;
animal_axes = axes(animal_fig);

% plant_fig = figure;
% plant_axes = axes(plant_fig);

plot(animal_axes, 1:numIterations, pollen_conc_counter, 'color', [1 0.5 0.17]);


% plot(animal_axes, 1:numIterations, animal_counter, 'color', [1 0.5 0.17]);
% plot(plant_axes, 1:numIterations, plant_counter);

% title(animal_axes, "Animal Population Throughout Simulation");
% ylabel(animal_axes, "Number of animals (incl. pollen-carrying)");
% xlabel(animal_axes, "Time step");
% ylim(animal_axes, [0, max(animal_counter)+5]); % Extend the maximum y limit

% title(plant_axes, "Plant Population Throughout Simulation");
% ylabel(plant_axes, "Number of plants (incl. growing and pollinated)");
% xlabel(plant_axes, "Time step");
% ylim(plant_axes, [0, max(plant_counter)+5]); % Extend the maximum y limit


disp("Simulation complete!");