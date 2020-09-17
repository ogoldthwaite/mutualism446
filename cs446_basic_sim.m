% CS446 -- Computational Modeling and Simulation II
% Owen Goldthwaite, Gautam Mitra, Lolo Niemiec
% September 15, 2020
% Initial Simulation for Mutualism


%% Simulation Parameters %%%

% Seed the random number generator for testing
rng_set = rng(1234567);

% Time-related variables
dt = 1;             % timestep, increment by days
simLength = 360;    % length of simulation: 1 year
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

prob_init_plant = 0.01; % initial probability a cell is plant
prob_init_animal = 0.1; % initial probability a cell is animal
prob_death_plant = 0.01; % probability plant death at each timestep
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

            % Default setting next updated cell to empty if it doesnt fit
            % any of these scenarios
            if(current_cell == EMPTY)
                updated_cell = EMPTY;                
 
% Do we need this?
%             % If current cell empty and have neighboring plants and animals
%             if(current_cell == EMPTY && plant_count~=0 && animal_count~=0)
% 
%                 % Cell becomes plant if more plants 
%                 if(plant_count > animal_count)
%                     updated_cell = PLANT;
%                 % Cell becomes animal if more animals
%                 elseif(animal_count > plant_count)
%                     updated_cell = ANIMAL;
%                 % Cell randomly chooses if equal amt of plants and animals
%                 elseif(animal_count == plant_count)
%                     choice_num = rand;
%                     if(choice_num < 0.5) % Random 50/50
%                         updated_cell = PLANT;
%                     else
%                         updated_cell = ANIMAL;
%                     end
%                 end
%             end  
                                 

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
                    growth_rate = .20;    % When should the plant grow
                    probability_growth = rand;
                    if (probability_growth < growth_rate)
                        updated_cell = PLANT;
                    else
                        updated_cell = GROWING_PLANT;
                    end
                end

            % If current cell is an animal then it requires at least half 
            % as many plant neighbors as animal neighbors + 1. Essentially 
            % meaning one plant can provide sustenance for 2 animals, this 
            % value of 2 is changeable.
%             if(current_cell == ANIMAL)
%                 total_animal_count = animal_count + 1;
%                 if(total_animal_count/animal_per_plant_value > plant_count)
%                     updated_cell = EMPTY;
%                 else
%                     updated_cell = ANIMAL;
%                 end
                        
%             end

            % If current cell is an animal and it is near a plant, then 
            % there is a chance it will pick up pollen. Then, it moves.
            elseif (current_cell == ANIMAL)
                total_animal_count = animal_count + 1;
                if (plant_count >= 1)
                    pollen_chance = rand;
                    if (pollen_chance < prob_pollination)
                        updated_cell = POLLINATED_ANIMAL;
                    end
                end
                
                % Move the animal to a new cell. If there aren't any plants
                % or animals in the surrounding cells, then it can move to 
                % any cell. If the surrounding cells are all occupied, the 
                % animal dies.
                num_neighbors = plant_count+animal_count;
                if (num_neighbors == 0)
                    % Move to any cell
                elseif ( (0 < num_neighbors) && ...
                        (num_neighbors < length(neighbors)) )
                    % Move animal to limited # of cells
                elseif (num_neighbors == length(neighbors))
                    % Animal can't move, dies
                    updated_cell = EMPTY;
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
                end
            
            
            % Updating next grid with the new cell value
            grids(row - 1, col - 1, frame) = updated_cell;

            end
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
        109/255 188/255 0;          % Plant Cell
        237/255 41/255  57/255;     % Growing Plant Cell
        .6       0      0;          % Animal Cell: grey blue
        1       .8      .8];        % Pollinated Animal Cell: light blue

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