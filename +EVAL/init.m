addpath(genpath('helpers'));
addpath(genpath('deblatting'))
addpath(genpath('deblatting_sroubek'))

warning('off','all');
global check_file;
if isempty(check_file)
	check_file = containers.Map;
end

