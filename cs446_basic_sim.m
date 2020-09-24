% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% September 15, 2020
% Initial Simulation for Mutualism


%% Simulation Parameters %%%

% Seed the random number generator for testing
rng_set = rng(1234567);

% Time-related variables
dt = 1;              % timestep, increment by days
simLength = 360;     % length of simulation: 1 year
numIterations = 1 + simLength/dt;
animation_fps = .5;  % Speed of visuazation

% Grid dimensions
row_count = 100; % width
col_count = 100; % length

%% Constants %%
EMPTY = 1;
PLANT = 2; 
GROWING_PLANT = 3;
ANIMAL = 4;
POLLINATED_ANIMAL = 5;
MOVED_ANIMAL = 6;

prob_init_plant = 0.1;  % initial probability a cell is plant
prob_init_animal = 0.1;  % initial probability a cell is animal
prob_death_plant = 0.01; % probability of plant death at each timestep
prob_pollination = 0.01; % probability animal will pick up/drop pollen

animal_per_plant_value = 2; % How many animals one plant can provide for


%% Set up grids
% Initialize grid to be all empty
grids = ones(row_count, col_count, numIterations) * EMPTY;

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
        end
    end
end

disp("Grids Initialized");

%% Main Simulation Loop
for frame = 2:numIterations
    
    %% Absorbing boundary condition
    % Create grids that are the size of the grid + 2 on each side
    extended_grid_size = size(grids( : , : , frame-1))+2;
    % Initialize all of the cells as empty
    extended_grid = ones(extended_grid_size) * EMPTY; 
    
    % Set the inside portion of the grids equal to the corresponding values 
    % from the previous timestep (a.k.a the previous frame)
    extended_grid(2:end-1, 2:end-1) = grids(:,:,frame-1);
    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to original(non-extended) grid
    for row = 2:row_count + 1
        for col = 2:col_count + 1

            % Store current cell
            current_cell = extended_grid(row, col);
            
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

            %% Update cell

            % Getting plant neighbor count and animal neighbor count
            plant_count = sum(neighbors == PLANT);
            animal_count = sum(neighbors == ANIMAL);
            num_neighbors = plant_count+animal_count;

            % Default setting next cell to empty. However, if 3 of the 
            % neighbor cells are animal cells, then one will move to 
            % this empty cell.          
            if(current_cell == EMPTY)
                % Something can move here
                if (animal_count == 3)
                    updated_cell = ANIMAL;
                else
                    % Nothing surrounding to move here
                    updated_cell = EMPTY;
                end  
                         
                
            % If current cell is a plant then there is a small chance for 
            % it to die at any given timestep
            elseif (current_cell == PLANT)
                death = rand;             
                if (death < prob_death_plant)
                    updated_cell = EMPTY;
                else
                    updated_cell = PLANT;
                end
                
            % If current cell is a growing plant then either use a 
            % probability or multiple growth stages to test if it should 
            % mature. There's a chance to die and it can't pollinate yet.
            elseif (current_cell == GROWING_PLANT) 
                death = rand;
                if (death < prob_death_plant)
                    updated_cell = EMPTY;
                else
                    growth_rate = .20;  % Likelihood plant grows this time
                    probability_growth = rand;
                    if (probability_growth < growth_rate)
                        updated_cell = PLANT;
                    else
                        updated_cell = GROWING_PLANT;
                    end
                end
                
            % If current cell is a pollinated animal,and neighbors a 
            % growing plant, there is a chance the animal pollinates that 
            % plant and then loses the pollen
            elseif (current_cell == POLLINATED_ANIMAL)
                pollinate_prob = rand;
                if ((pollinate_prob<prob_pollination) && (plant_count>=1))
                    updated_cell = ANIMAL;
                    % Pollinate the plant
                else
                    % Just keep the pollinated animal moving
                    updated_cell = POLLINATED_ANIMAL;
                end

            % If current cell is an animal and it is near a plant, then 
            % there is a chance it will pick up pollen. 
            elseif (current_cell == ANIMAL)
                total_animal_count = animal_count + 1;
                pollen_pickup = rand;
                if (num_neighbors==8) || (animal_count>2)
                    % Now that one of the animals moved to the empty 
                    % cell, we must remove them from their current cell.
                    % Also, death by overcrowding.
                    updated_cell = EMPTY;
                else
                    if (plant_count >= 1) && (pollen_pickup < prob_pollination)
                        updated_cell = POLLINATED_ANIMAL;
                    else
                        updated_cell = ANIMAL;
                    end
                end

 
            end
            
            % Updating next grid with the new cell value
            grids(row - 1, col - 1, frame) = updated_cell;
        end
    end
end

disp("All grids calculated");

%% Visualize the grid

% Create the window for the animation
viz_fig = figure;
viz_axes = axes(viz_fig);

% Set the colors
map = [ 1       1       1;          % Empty Cell: white
        40/255  191/255 118/255;    % Plant Cell: green
        1       1       0;          % Growing Plant Cell: yellow
        125/255 100/255 70/255;     % Animal Cell: brown
        0       0       1];          % Pollinated Animal Cell: blue
colormap(viz_axes, map); 

% Remove axis labels, make aspect ratio look good, and maintain that state
axis off;
axis equal;
hold on;

disp("Drawing...");
for i = 1:numIterations
    
    % Following line allows for frame-by-frame viewing of grid
    %w = waitforbuttonpress;

    % Turn each grid into an image
    image(viz_axes,grids(:, :, i));

    pause(1/animation_fps);
end

disp("Simulation complete!");
