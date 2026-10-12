source(file.path("R", "utils.R"))
source(file.path("R", "config.R"))
simulation_config <- mppi_simulation_config(B = 1000L, K = 5L)
mppi_validate_config(simulation_config)
