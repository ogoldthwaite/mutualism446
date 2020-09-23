% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% September 15, 2020
% Initial Simulation for Mutualism


%% Simulation Parameters %%%

% Seed the random number generator for testing
rng_set = rng(123456789);

% Time-related variables
dt = 1;             % timestep, increment by days
simLength = 360;    % length of simulation: 1 year
numIterations = 1 + simLength/dt;
animation_fps = 5;  % Speed of visuazation

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
prob_init_animal = 0.01; % initial probability a cell is animal
prob_plant_death = 0.05; % probability plant death at each timestep
prob_pollination = 0.5; % probability animal will pick up/drop pollen
prob_pollen_production = 0.2; % chance that a plant will produce pollen
prob_pollen_spread = 0.15; % chance that pollen will spread to an adjacent empty cell
prob_animal_reproduce = .2 % chance that 2 animals will reproduce if conditions are good
prob_random_move = 0.15; % chance that an animal will move randomly with no stimuli

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
            grids(row, col, 1) = GROWING_PLANT;

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

            % List of Neighbor coordinates
            neighbor_coords = [ row-1 col; row col-1; row+1 col; row col+1; row-1 col-1; row+1 col-1; row-1 col+1; row+1 col+1 ];                    


            %% Update cell

            % Getting different counts of neighbors
            empty_count = sum(neighbors == EMPTY);
            pollen_count = sum(neighbors == POLLEN);

            norm_plant_count = sum(neighbors == PLANT);
            poll_plant_count = sum(neighbors == POLLINATED_PLANT);
            norm_animal_count = sum(neighbors == ANIMAL);
            poll_animal_count = sum(neighbors == POLLINATED_ANIMAL);
    
            animal_count = norm_animal_count + poll_animal_count; % Total animal count
            plant_count = poll_plant_count + norm_plant_count; % Total plant count


            % Empty cell behavior
            % If empty cell neighboring a pollinated plant than it will become a growing plant, highest priority
            % If Empty cell is neighboring pollen it will become pollen, this
            % will have to be restricted so it doesnt increase exponentially somehow
            % If empty cell is neighboring a pollen producing plant it will become pollen
            % Else it will remain empty
            if(current_cell == EMPTY)
                if(poll_plant_count > 0)
                    updated_cell = GROWING_PLANT;
                elseif(pollen_count > 0 && rand < prob_pollen_spread)
                    updated_cell = POLLEN;
                elseif(plant_count > 0 && rand < prob_pollen_production)
                    updated_cell = POLLEN;
                elseif(animal_count > 1 && plant_count > 0 && rand < prob_animal_reproduce)
                    updated_cell = ANIMAL;
                else
                    updated_cell = EMPTY;
                end
            
            % Animal Behavior
            % If next to plant has a chance to become pollen carrying animal, this takes a timestep so it cannot move and do this
            % If next to pollen will become an empty cell because the pollen will become an animal (moving towards flower)
            % If totally surrounded by other animals/plants it will die and become empty
            % If no pollen currently nothing will happen (random movement later?)
            elseif (current_cell == ANIMAL)
                if (plant_count > 0 && rand < prob_pollination) % Pollinated plants can pollinate too, hermaphroditic currently
                    updated_cell = POLLINATED_ANIMAL; % Gets pollen from plant
                elseif (pollen_count > 0)
                    updated_cell = EMPTY; % Moving to pollen (done in pollen check)
                elseif (plant_count + animal_count == 8)
                    updated_cell = EMPTY; % Death
                else
                    updated_cell = ANIMAL; % Nothing happens
                end
            
            % Pollinated animal behavior
            % If next to plant has a chance to lose the pollen it's carrying
            % If next to pollen cell will become empty to move like normal animal
            % If totally surrounded will also die like normal animal
            % If no special conditonis are met will remain same like before
            elseif (current_cell == POLLINATED_ANIMAL)
                if (plant_count > 0 && rand < prob_pollination)
                    updated_cell = ANIMAL; % Loses pollen
                elseif (pollen_count > 0)
                    updated_cell = EMPTY; % Moving to pollen (done in pollen check)
                elseif (plant_count + animal_count == 8)
                    updated_cell = EMPTY; % Death
                else
                    updated_cell = POLLINATED_ANIMAL;
                end

            % Pollen cell state behavior
            % If next to an animal or pollinated animal will become one of them, pollinated higher priority
            % Elseif next to empty cell will become empty (empty will then become pollen)
            elseif (current_cell == POLLEN)
                if(animal_count == 1)
                    if (poll_animal_count > 0)
                        updated_cell = POLLINATED_ANIMAL;
                    elseif (norm_animal_count > 0)
                        updated_cell = ANIMAL;
                    end
                elseif(empty_count > 0)
                    updated_cell = EMPTY;
                end
            
            % Growing plant behavior
            % Just becomes a normal plant next timestep
            elseif (current_cell == GROWING_PLANT)
                updated_cell = PLANT;

            % Normal plant state behavior
            % If next to pollinated animal then has chance to become a pollinated plant
            % Else just has small random chance to die
            % Else if nothing happens, just remains a plant
            elseif (current_cell == PLANT)
                if(poll_animal_count > 0 && rand < prob_pollination)
                    updated_cell = POLLINATED_PLANT;
                elseif(rand < prob_plant_death)
                    updated_cell = EMPTY;
                else
                    updated_cell = PLANT; % Just remains a plant
                end

            % Pollinated Plant state behavior
            % If next to an empty cell then turns into a normal plant since new plant will grow in said empty cell
            % Else just remains pollinated, waiting for empty location
            elseif (current_cell == POLLINATED_PLANT)
                if(empty_count > 0)
                    updated_cell = PLANT;
                else
                    updated_cell = POLLINATED_PLANT;
                end
                          
            end
            % Updating next grid with the new cell value
            grids(row-1, col-1 , frame) = updated_cell;
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
        0  1  0;                    % Growing Plant Cell
        109/255 188/255 0;          % Plant Cell
        109/255 120/255 0;          % Pollinated plant
        125/255 100/255  70/255;          % Animal Cell: grey blue 
        0       0      0;          % Pollinated Animal Cell: light blue
        230/255 255/255 0];         % Pollen color        

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
    image(viz_axes, grids(:, :, i));

    pause(1/animation_fps);
end

disp("Simulation complete!");