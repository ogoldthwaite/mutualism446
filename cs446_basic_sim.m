% CS446 -- Computational Modeling and Simulation II
% September 15, 2020
%
% Initial Simulation for mutualism
%

%% Simulation Parameters %%%

% Seed the random number generator for testing
rng_set = rng(1234567);

% Time-related variables
dt = 1; % timestep
simLength = 100; % length of simulation
numIterations = 1 + simLength/dt;
animation_fps = 5; % Speed of visuazation

% Grid dimensions
row_count = 100; % width
col_count = 100; % length

%% Constants %%
EMPTY = 1;
PLANT = 2; 
ANIMAL = 3;

prob_init_plant = 0.01; % initial probability a cell is plant
prob_init_animal = 0.1; % initial probability a cell is animal

animal_per_plant_value = 2; % Basically how many animals one plant can provide for


%% Set up grids
% Initialize grid to be all empty
grids = ones(row_count, col_count, numIterations) * EMPTY;

for row = 1:row_count
    for col = 1:col_count
        % grid initialization, plants first
        plant_chance = rand;
        animal_chance = rand;
        % Checking for tr33s
        if plant_chance < prob_init_plant
            grids(row, col, 1) = PLANT;

        % Now checking for animals
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
    extended_grid = ones(extended_grid_size) * EMPTY; % initialize all as empty
    
    % Set the inside portion of the grids equal to the corresponding values from
    % the previous timestep (a.k.a the previous frame)
    extended_grid(2:end-1, 2:end-1) = grids(:,:,frame-1);
    
    %% Loop for updating each cell in the grid
    % Loop over the indices corresponding to the original (non-extended) grid
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

            % Default setting of next updated cell to empty if it doesnt fit
            % any of these scenarios
            if(current_cell == EMPTY)
                updated_cell = EMPTY;
            end


            % If current cell is empty and there are some neighboring plants and animals
            if(current_cell == EMPTY && plant_count ~= 0 && animal_count ~= 0)

                % Cell becomes plant if more plants, animal if more animals or a random choice if equal
                if(plant_count > animal_count)
                    updated_cell = PLANT;
                elseif(animal_count > plant_count)
                    updated_cell = ANIMAL;
                elseif(animal_count == plant_count)
                    choice_num = rand;
                    if(choice_num < 0.5) % Random 50/50
                        updated_cell = PLANT;
                    else
                        updated_cell = ANIMAL;
                    end
                end
            end

            % If current cell is an plant then it will survive no matter what
            % because it just needs sun to live
            if(current_cell == PLANT)
                updated_cell = PLANT;
            end

            % If current cell is an animal then it requires at least half as many
            % plant neighbors as animal neighbors + 1. Essentially meaning one plant
            % can provide sustenance for 2 animals, this value of 2 is changeable
            if(current_cell == ANIMAL)
                total_animal_count = animal_count + 1;
                if(total_animal_count/animal_per_plant_value > plant_count)
                    updated_cell = EMPTY;
                else
                    updated_cell = ANIMAL;
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
empty_color = [1, 1, 1];
plant_color = [109/255, 188/255, 0];
animal_color = [237/255 41/255 57/255];

colormap(viz_axes, [empty_color; plant_color; animal_color]); 

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


