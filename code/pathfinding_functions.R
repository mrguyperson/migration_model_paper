
# function to select appropriate columns and filter depths  --------

prep_habitat_df <- function(df, fish_parm, habitat_parm, min_depth) {
  df %>%
    as.data.table() %>% 
    setnames(old = grep("dista|lat|velo", names(.), value = TRUE), new = c("distance", "lat_dist", "velocity")) %>%
    .[order(distance, lat_dist)] %>%
    # filter cells accessible to the species based on flow
    .[.[[grep(fish_parm$specie, names(.), value = TRUE)]] == 1,] %>%
    # filter columns by depth
    .[.$depth > min_depth, ] %>%
    adjust_velocity_for_wall_factor(fish_parm, habitat_parm) %>%
    # select columns
    .[, .SD, .SDcols = grep("dista|lat|velo|dept|area|temp", names(.), value = TRUE)] %>%
    .[, id := seq(1, nrow(.))] %>% 
    as_tibble()
}

# changes velocity to reflect the depth of benthic fish -------------------

adjust_velocity_for_wall_factor <- function(df, fish_parm, habitat_parm) {
  if (fish_parm$benthic_fish == 1) {
    df %>%
      .[depth >= habitat_parm$ben_vel_height, velocity := velocity * habitat_parm$base_wall_factor]
  } else {
    df
  }
}

# provides minimal columns of prepped df to provide environmental info for destination cells  ---------------------

trim_df_for_joining <- function(prepped_df) {
  prepped_df %>%
    # select columns
    .[, c("distance", "lat_dist", "velocity", "area", "id")] %>%
    setnames(old = c("id", "velocity"), new = c("to", "to_velocity"))
}

# adds coordinates for possible destination cells depending on forward or diagonal movement --------

get_destination_coords <- function(prepped_df, cell_width_m) {
  dist_opt <- c(cell_width_m)
  lat_opt <- c(cell_width_m, -cell_width_m, 0)
  df_w_cell_distances <- expand_grid(prepped_df, node_distance = dist_opt, node_lat_dist = lat_opt)
  df_w_cell_distances %>%
    mutate(
      node_dir = fifelse(node_lat_dist == 0, "f", "s"),
      across(
        node_distance:node_lat_dist,
        .fns = ~ .x + get(str_remove(cur_column(), "node_"))
      )
    )
}

# adds destinations to prepped df and environment values for all cells --------

add_destination_info_to_start_cells <- function(prepped_df, cell_width_m){
  df_trimmed <- trim_df_for_joining(prepped_df)
  
  get_destination_coords(prepped_df, cell_width_m) %>%
    # use inner join of the trimmed df to add data (e.g., velocity) for the destination coordinates and filter
    # out destinations coordinates that don't actually exist
    inner_join(df_trimmed, by = c("node_distance" = "distance", "node_lat_dist" = "lat_dist")) %>%
    # renaming to clarify which cells are origin ("from") and which are destinations ("to")
    setnames(old = c("area.x", "area.y", "id", "velocity"), new = c("from_area", "to_area", "from", "from_velocity"))
  
}

# calculates metabolic rate at a given swim speed and temp ----------------

get_metabolic_rate <- function(fish_parameters, swim_speed_m_per_s, temperature_C) {
  int <- fish_parameters$met_int
  lm <- fish_parameters$met_lm
  lt <- fish_parameters$met_lt
  v <- fish_parameters$met_v
  lm_lt <- fish_parameters$met_lm_lt
  lm_v <- fish_parameters$met_lm_v
  t <- fish_parameters$met_t
  lm_t <- fish_parameters$met_lm_t
  sqv <- fish_parameters$met_sqv
  
  m_per_cm <- 0.01
  hours_per_day <- 24
  minutes_per_hour <- 60
  seconds_per_minute <- 60
  
  body_mass_g <- fish_parameters$fish_mass_g

  log_metabolic_rate <-
    int +
    lm * log(body_mass_g) +
    lt * log(temperature_C) +
    v * swim_speed_m_per_s  +
    lm_lt * log(body_mass_g) * log(temperature_C) +
    lm_v * log(body_mass_g) * swim_speed_m_per_s +
    t * temperature_C +
    lm_t * log(body_mass_g) * temperature_C +
    sqv * sqrt(swim_speed_m_per_s)
  
  metabolic_rate_j_per_day <- exp(log_metabolic_rate)

  metabolic_rate_j_per_day / (hours_per_day * minutes_per_hour * seconds_per_minute)
}

# functions for basic fish parameters -------------------------------------

get_fish_body_mass <- function(params) {
  a <- params$length_mass_a
  b <- params$length_mass_b
  a * params$eg_adult_length^b
}

get_ucrit <- function(params, temperature_C) {
  # ucrit in body lengths per second
  a <- params$ucrit_a
  b <- params$ucrit_b
  c <- params$ucrit_c
  d <- params$ucrit_d

  ucrit_in_body_lengths <- 
    (a / params$eg_adult_length + b) *
    (1 + (c - temperature_C) / (c - d)) * ((temperature_C / c)^(c / (c - d)))

  m_per_cm <- 0.01
  # ucrit in m per second; body length is in cm
  ucrit_in_body_lengths *  params$eg_adult_length * m_per_cm
}

get_max_swim_speed <- function(params) {
  m_per_cm <- 0.01
  body_length_m <- params$eg_adult_length * m_per_cm
  0.4 + 7.4 * body_length_m
}

# minimum depth a cell can have to be usable by a fish --------------------

get_min_depth <- function(params, factor) {
  m_per_cm <- 0.01
  factor * params$eg_adult_length * m_per_cm
}

get_cell_length <- function(df, cell_width_m) {
  df %>%
    mutate(across(contains("area"),
      ~ .x / cell_width_m,
      .names = "{str_replace(.col, 'area', 'cell_length')}"
    ))
}

get_cell_hypotenuse <- function(df, cell_width_m) {
  df %>%
    mutate(
      across(contains("cell_length"),
        ~ sqrt(.x^2 + cell_width_m^2),
        .names = "{str_replace(.col, 'length', 'hypotenuse')}"
      )
    )
}

# ratio of cell length to hypotenuse --------------------------------------

# used for scaling diagonal swim speed
get_ratios <- function(df) {
  df %>%
    mutate(
      across(contains("hypotenuse"),
        ~ fifelse(node_dir == "f", 
                  1, 
                  get(str_replace(cur_column(), "hypotenuse", "length")) / .x),
        .names = "{str_replace(.col, 'cell_hypotenuse', 'ratio')}"
      )
    )
}

# calculates swim speed of all fish ---------------------------------------

# based on a numerical estimate of the optimal speed for the Martin et al. model 
# substituted with the FHAST metabolic equations

get_cell_swim_speeds <- function(df, fish_parm) {
  ucrit <- fish_parm$ucrit_m_per_s
  ucrit_cutoff <- fish_parm$ucrit_cutoff_m_per_s
  swim_speed_max <- fish_parm$swim_speed_max_m_per_s
  min_vel_ucrit <- fish_parm$min_vel_ucrit_m_per_s
  max_vel_ucrit <- fish_parm$max_vel_ucrit_m_per_s
  max_water_vel <- fish_parm$max_water_vel_m_per_s
  
  calculate_slope_and_intercept <- function(x1,x2,y1,y2){
   slope <- ( y2 - y1 )/ (x2 - x1)
   intercept <- y1 - slope * x1
   list(intercept=intercept, slope=slope)
  }
  min_vel_params <- calculate_slope_and_intercept(0, min_vel_ucrit, ucrit_cutoff, ucrit)
  max_vel_params <- calculate_slope_and_intercept(max_vel_ucrit, max_water_vel, ucrit, swim_speed_max)
  
  df %>%
    mutate(
      across(ends_with("velocity"),
        .fns = ~ fcase(
          # # multiply swim speed by the ratio of the cell length and cell hypotenuse
          ucrit > ucrit_cutoff &
            min_vel_ucrit * get(str_replace(cur_column(), "velocity", "ratio")) > .x, .x / get(str_replace(cur_column(), "velocity", "ratio")) * min_vel_params$slope + ucrit_cutoff,
          # set ucrit swim speed
          fish_parm$max_vel_ucrit * get(str_replace(cur_column(), "velocity", "ratio")) >= .x, ucrit,
          # burst swim speeds between ucrit and max speed are +50% higher than water velocity
          swim_speed_max > (max_vel_params$slope * .x) / get(str_replace(cur_column(), "velocity", "ratio")) + max_vel_params$intercept, (max_vel_params$slope * .x) / get(str_replace(cur_column(), "velocity", "ratio")) + max_vel_params$intercept,
          swim_speed_max <= .x, NA_real_,
          default = swim_speed_max
        ),
        .names = "{str_replace(.col, 'velocity', 'swim_speed')}"
      )
    )
}

# Martin et al. model of movement cost vs. distance -----------------------

get_cost_of_travel <- function(swim_speed_m_per_s, 
                               water_velocity_m_per_s, 
                               ucrit_m_per_s,
                               fish_parm, 
                               temperature_C,

                               ratio) {

  anaerobic_fuel_recovery_parameter <- 1.82

  
  base_metabolic_rate <- get_metabolic_rate(fish_parm, 0, temperature_C)
  critical_metabolic_rate <- get_metabolic_rate(fish_parm, ucrit_m_per_s, temperature_C)
  total_metabolic_rate <- get_metabolic_rate(fish_parm, swim_speed_m_per_s, temperature_C)

  (total_metabolic_rate + # aerobic component
      base_metabolic_rate *    pmax(0, anaerobic_fuel_recovery_parameter * (total_metabolic_rate - critical_metabolic_rate) / (critical_metabolic_rate - base_metabolic_rate)) + # energy spent due to waiting to recover from anerobic activity
      pmax(0, (anaerobic_fuel_recovery_parameter - 1) * (total_metabolic_rate - critical_metabolic_rate))) / # cost of recovering fuel after anerobic activity
    (ratio * swim_speed_m_per_s - water_velocity_m_per_s)
}

# calculation of cheapest paths between a given set of start and end points --------
# based on Djikstra's algorithm
# returns the distance and lat distance coordinates of cells and the number of paths that went through that cell

get_paths <- function(prepped_df, graph, from, to, weights) {
  map(from, ~ shortest_paths(graph = graph, weights = weights, mode = "out", from = .x, to = to)$vpath) %>%
    unlist() %>%
    as.data.table() %>%
    .[, .N, by = .] %>%
    .[as.data.table(prepped_df), on = c("." = "id"), nomatch = 0] %>%
    setnames(old = "N", new = "num_paths") %>%
    .[, c("distance", "lat_dist", "num_paths")] %>%
    as_tibble()
}

# calculates all costs and paths for a given species ----------------------

get_paths_and_costs <- function(hab_df, 
                                fish_parm, 
                                habitat_parm, 
                                fish_id) {
  # set params --------------------------------------------------------------
  fish_parm <- map(fish_parm, ~.x[[fish_id]])
  species <- fish_parm$specie
  date <- hab_df$date[[1]]
  reach_temp_C <- hab_df$temp[[1]]
  min_depth_m <- get_min_depth(factor = 0.5, fish_parm)
  cell_width_m <- habitat_parm$resolution
  fish_parm[["ucrit_m_per_s"]] <- get_ucrit(fish_parm, reach_temp_C)
  fish_parm[["min_vel_ucrit_m_per_s"]] <- get_estimate_for_linear_model(
    fish_parm$pars_min_water_vel_ucrit_int, fish_parm$pars_min_water_vel_ucrit_slope, fish_parm$ucrit_m_per_s)
  fish_parm[["max_vel_ucrit_m_per_s"]] <- get_estimate_for_linear_model(
    fish_parm$pars_max_water_vel_ucrit_int, fish_parm$pars_max_water_vel_ucrit_slope, fish_parm$ucrit_m_per_s)

  # determine cell relationships and costs ----------------------------------
  prepped_df <- prep_habitat_df(hab_df, fish_parm, habitat_parm, min_depth_m)
  relations <- add_destination_info_to_start_cells(prepped_df, cell_width_m) %>% 
    get_cell_length(cell_width_m) %>%
    get_cell_hypotenuse(cell_width_m) %>%
    # hypotenuse multiplier for determining forward velocity component during diagonal movement
    get_ratios() %>%
    get_cell_swim_speeds(
      fish_parm
    ) %>%
    # drop any cells that would have negative overground velocity values
    drop_na() %>%
    mutate(
      energy_cost = rowSums(across(ends_with("swim_speed"),
                                   .fns = ~ get_cost_of_travel(.x, 
                                                               get(str_replace(cur_column(), "swim_speed", "velocity")), 
                                                               fish_parm$ucrit_m_per_s,
                                                               fish_parm, 
                                                               reach_temp_C,
                                                               get(str_replace(cur_column(), "swim_speed", "ratio")) 
                                                               ) *
                                     get(str_replace(cur_column(), "swim_speed", "cell_length"))
      )),
    ) %>%
    .[, c("from", "to", "energy_cost")]

  # build graph and calculate distances ------------------------------------------

  # get ids for all cells
  actors <- prepped_df$id

  # create a graph of all cell relationships
  g <- graph_from_data_frame(relations, vertices = actors)

  # # select cells to start at
  from <- prepped_df[prepped_df$distance == min(prepped_df$distance), ]$id
  to <- prepped_df[prepped_df$distance == max(prepped_df$distance), ]$id


  # costs of moving between cells are equal to the energy cost to the fish
  weights <- relations$energy_cost

  # calculate shortest paths based on energy costs
  distances <- distances(g, weights = weights, mode = "out", v = from, to = to)

  # convert matrix to vector
  distance_vector <- unlist(as.list(distances))
  # check for div by 0
  distance_vector <- distance_vector[is.finite(distance_vector)]

  paths <- get_paths(prepped_df, g, from, to, weights)

  
  data.table(
    date = date, 
    species = species, 
    energy_cost = list(distance_vector),
    paths = list(paths)
  )

}

# runs through the main function for all species --------------------------

get_path_min_costs_all_species <- function(hab_df, fish_parm, habitat_parm, fish_schedule) {
  species <- fish_schedule$species
  fish_ids <- match(species, fish_parm$specie)
  pathfinding_table <- map_dfr(fish_ids, ~ get_paths_and_costs(hab_df, fish_parm, habitat_parm, .x))
  pathfinding_table[as.data.table(fish_schedule), on = "species", nomatch = 0]
}

##### functions for getting swim speed parameters -----------------------------

# adds some additional parameters to the fish_parm list -------------------

calculate_adult_parameters <- function(fish_parm){
  num_species <- length(fish_parm$specie)
  fish_parm_lol <- map(1:num_species, ~ map(fish_parm, .x))
  fish_parm[["fish_mass_g"]] <- map_dbl(fish_parm_lol, get_fish_body_mass)
  fish_parm[["swim_speed_max_m_per_s"]] <- map_dbl(fish_parm_lol, get_max_swim_speed)
  return(fish_parm)
}

# makes a simulated environment with different temp and water velo --------

make_environment_dt <- function(swim_speed_max){
  
  CJ(
    temperature = seq(1,25,1),
    velocity = seq(0, swim_speed_max, swim_speed_max/50))
}


# calculates the parameters need to estimate the optimal swim speed curve --------

get_swim_speed_model_params <- function(fish_parm, fish_id){
  fish_parm <- map(fish_parm, ~.x[[fish_id]])
  fish_length <- fish_parm$eg_adult_length
  fish_mass <- fish_parm$fish_mass_g
  swim_speed_max <- fish_parm$swim_speed_max_m_per_s
  
  dt <- make_environment_dt(swim_speed_max) %>% 
    .[, ucrit := get_ucrit(fish_parm, temperature)]
  dt$swim_speed_martin <- pmap_dbl(list(dt$velocity, dt$temperature, dt$ucrit), 
                                   ~ optimize(get_cost_of_travel, 
                                              interval =  c(..1 + 1e-6, swim_speed_max), 
                                              water_velocity_m_per_s = ..1, 
                                              fish_parm = fish_parm, 
                                              temperature_C = ..2, 
                                              ucrit_m_per_s = ..3,
                                              ratio = 1)$minimum)
  
  min_speeds <- dt[, .(min_speed = min(swim_speed_martin)), by = "temperature"]
  
  ucrit_cutoff <- min_speeds[dt, on = "temperature"] %>%
    .[min_speed == max(min_speed)] %>%
    .$min_speed %>%
    .[[1]]
  
  min_vel_ucrit <- dt[signif(swim_speed_martin,4) == signif(ucrit,4), .(min = min(velocity)), by = "temperature"]
  parameters_for_min_water_velocity_at_ucrit <- min_vel_ucrit[dt, on = "temperature", nomatch = 0] %>%
    .[, c("ucrit", "min")] %>%
    .[ucrit > ucrit_cutoff,] %>% 
    lm(min ~ ucrit, data = .) %>%
    .[[1]]

  max_vel_ucrit <- dt[signif(swim_speed_martin,4) == signif(ucrit,4), .(max = max(velocity)), by = "temperature"]
  
  water_velocity_at_max_burst <- dt[signif(swim_speed_martin,4) == signif(swim_speed_max, 4)] %>% 
    .$velocity %>% 
    min()
  
  parameters_for_max_water_velocity_at_ucrit <- max_vel_ucrit[dt, on = "temperature", nomatch = 0] %>%
    .[, c("ucrit", "max")] %>%
    lm(max ~ ucrit, data = .) %>%
    .[[1]]
  
  list(
    max_water_vel_m_per_s = water_velocity_at_max_burst,
    ucrit_cutoff_m_per_s = ucrit_cutoff,
    pars_min_water_vel_ucrit_int = parameters_for_min_water_velocity_at_ucrit[[1]],
    pars_min_water_vel_ucrit_slope = parameters_for_min_water_velocity_at_ucrit[[2]],
    pars_max_water_vel_ucrit_int = parameters_for_max_water_velocity_at_ucrit[[1]],
    pars_max_water_vel_ucrit_slope = parameters_for_max_water_velocity_at_ucrit[[2]]
  )
}

# runs through the main swim speed param function for all species ---------

get_swim_speed_parameters_for_all_species <- function(fish_parm){
  fish_ids <- 1:length(fish_parm$specie)
  future_map(fish_ids, 
             get_swim_speed_model_params, 
             fish_parm = fish_parm) %>% 
    pmap(c)
}