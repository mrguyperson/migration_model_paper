library(tidyverse)
library(here)
library(igraph)
library(sf)
library(microbenchmark)
library(data.table)

source(here("scripts/R/migration/R/pathfinding_functions.R"))

# start make path -----------------------------------------------------

shape <- st_read(here("scripts/R/migration/data/test_big_spp_depths.shp"), quiet = TRUE) 



daily_data <- daily_file[1,]
hab_df <-  get_daily_v_and_d(daily_data, habitat, flows, 10) %>% left_join(grid_file) %>%
  as_tibble() %>%
  select(-geometry)

cols_to_rename <- colnames(shape)[1:length(colnames(shape)) - 1]
new_names <- colnames(hab_df)

shape %>% setnames(old = cols_to_rename, new = new_names)

# # set params --------------------------------------------------------------

fish_id <- 3
# determine cell relationships and costs ----------------------------------
fish_parm <- map(fish_parm, ~.x[[fish_id]])
species <- fish_parm$specie
date <- hab_df$date[[1]]
# body_length_cm <- fish_parm$eg_adult_length
# body_mass_g <- get_fish_body_mass(fish_parm, body_length_cm)
reach_temp_C <- hab_df$temp[[1]]
min_depth_m <- get_min_depth(factor = 0.5, fish_parm)
cell_width_m <- habitat_parm$resolution
fish_parm[["ucrit_m_per_s"]] <- get_ucrit(fish_parm, reach_temp_C)
fish_parm[["max_vel_ucrit_m_per_s"]] <- get_estimate_for_linear_model(
  fish_parm$pars_max_water_vel_ucrit_int, fish_parm$pars_max_water_vel_ucrit_slope, fish_parm$ucrit_m_per_s)
# swim_speed_max_m_per_s <- get_max_swim_speed(body_length_cm)

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
                                                             #   body_mass_g,
                                                             # body_length_cm,
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

map(from, ~ shortest_paths(graph = g, weights = weights, mode = "out", from = .x, to = to)$vpath)

paths <- get_paths(prepped_df, g, from, to, weights)


list(
  species = species,
  date = date,
  energy_cost = distance_vector,
  paths = paths
)

path <- tibble(id = as.numeric(shortest_paths(g, weights = weights, mode = "out", from = 4, to = to[[3]])[[1]][[1]]))

shape_for_join <- prep_habitat_df(shape, fish_parm, habitat_parm, min_depth_m)

left_join(path, shape_for_join) %>% st_write(here("scripts/R/migration/data/path.shp"), append = FALSE)


# end make path -----------------------------------------------------------



get_opt_speed <- function(vel_tib, pcrit, swim_speed_max, p_m, temperature){
  vel_tib %>% 
    mutate(swim_speed = signif(map_dbl(velocity, ~ optimize(get_cot_jesse, interval =  c(.x + swim_speed_opt/1000, swim_speed_max * 2), water_velocity = .x, p_m, p_crit, y_r, fish_parm, body_mass, temperature, body_length, ratio = 1, fish_id)$minimum), 3)) %>% 
    select(swim_speed)
}

get_min_speed <- function(data){
  data %>% 
    summarize(
     # diff = max(velocity) - min(velocity), 
      min_speed = min(swim_speed)
      )
}

get_diff <- function(data){
  data  %>% 
    mutate(across(.fns = ~ signif(.x, 4))) %>% 
    filter(swim_speed == ucrit) %>%
    summarize(
      max = max(velocity),
      # min_speed = min(swim_speed)
    )
}

get_vel_at_max_speed <- function(data, swim_speed_max){
  data %>% 
    mutate(across(.fns = ~ signif(.x, 3))) %>% 
    filter(swim_speed == signif(swim_speed_max, 3)) %>% 
    summarize(vel_at_max_speed = min(velocity))
}

get_estimate_for_linear_model <- function(intercept, slope, x){
  intercept + slope * x
}

get_intercept <- function(ucrit, swim_speed_max, v_max, v_max_ucrit){
  ucrit - ((swim_speed_max - ucrit) / (v_max - v_max_ucrit)) * v_max_ucrit
}

get_min_speed_for_high_ucrit <- function(df){
  df %>% 
    filter(min_speed == max(min_speed)) %>% 
    filter(ucrit == min(ucrit)) %>% 
    select(min_speed) %>% 
    distinct() %>% 
    pull(min_speed)
}


fish_parm <- fish_parm_temp %>% 
  rename(species_temp = species) %>% 
  pivot_longer(cols=c(-species_temp), names_to="specie")%>%
  pivot_wider(names_from=c(species_temp)) %>% 
  as.list() %>% 
  calculate_adult_parameters() %>% 
  append(get_swim_speed_parameters_for_all_species(.))

fish_id <- 1
fish_parm <- map(fish_parm, ~.x[[fish_id]])
fish_length <- fish_parm$eg_adult_length
fish_mass <- fish_parm$fish_mass_g
swim_speed_max <- fish_parm$swim_speed_max_m_per_s

dt <- make_environment_dt(swim_speed_max) %>% 
  .[, ucrit := get_ucrit(fish_parm, temperature)]
dt$swim_speed_martin <- pmap_dbl(list(dt$velocity, dt$temperature, dt$ucrit), 
                                 ~ optimize(get_cost_of_travel, 
                                            interval =  c(..1 + 1e-3, swim_speed_max), 
                                            water_velocity_m_per_s = ..1, 
                                            fish_parm = fish_parm, 
                                            #    body_mass_g = fish_mass, 
                                            temperature_C = ..2, 
                                            #   body_length_cm = fish_length, 
                                            ucrit_m_per_s = ..3,
                                            ratio = 1)$minimum)
temp <- 15
test <- dt[temperature == temp] %>% 
  .[, ratio := 1] 

fish_parm[['ucrit_m_per_s']] <- test[["ucrit"]][[1]]
fish_parm[["min_vel_ucrit_m_per_s"]] <- get_estimate_for_linear_model(
  fish_parm$pars_min_water_vel_ucrit_int, fish_parm$pars_min_water_vel_ucrit_slope, fish_parm$ucrit_m_per_s)
fish_parm[["max_vel_ucrit_m_per_s"]] <- get_estimate_for_linear_model(
  fish_parm$pars_max_water_vel_ucrit_int, fish_parm$pars_max_water_vel_ucrit_slope, fish_parm$ucrit_m_per_s)
ucrit <- test[["ucrit"]][[1]]

test %>% 
  get_cell_swim_speeds(fish_parm) %>% 
  drop_na() %>% 
  mutate(swim_type = fcase(swim_speed < ucrit, "sub-Ucrit",
                           swim_speed == ucrit, "Ucrit",
                           default = "burst"),
         swim_type = factor(swim_type, levels = c("sub-Ucrit", "Ucrit", "burst"))) %>% 
  ggplot() +
  geom_line(aes(velocity, swim_speed_martin), linewidth = 1) +
 geom_point(aes(velocity, swim_speed, color = swim_type), size = 2) +
  theme_classic(base_size = 25) +
  xlim(0, 3) +
  ylim(1,4) +
  xlab("water velocity (m/s)") +
  ylab("optimal swim speed (m/s)") +
  labs(color = "swimming type") +
  scale_color_manual(labels = c(bquote("sub-U"["crit"]), bquote("U"["crit"]), "burst"),
                       values = c("#D55E00","#009E73","#0072B2")) +
  theme(
    legend.title = element_text(size = 14),
    legend.position = c(0.8,0.28),
    legend.text = element_text(size = 14)
  )


ggsave(here("scripts","R","migration","output","opt_swim_speed_w_estimates.png"),
       device = "png",
       dpi = 300,
       height = 5,
       width = 5
)
# stuff for temp and size independent parameters --------------------------
#TODO functions need to be updated to the latest!
fish_id <- 1

lengths <- tibble(body_length = seq(50,100,5)) %>% 
  mutate(body_mass = get_fish_body_mass(fish_parm, body_length, fish_id))


full3 <- lengths %>% 
  expand_grid(temperature = 1:25) %>% 
  mutate(ucrit = get_ucrit(fish_parm, temperature, body_length, fish_id),
         swim_speed_max = get_max_swim_speed(body_length)
  ) %>% 
  expand_grid(velocity = seq(0,7,0.01))
tic()
full4 <- full3 %>%  
  mutate(swim_speed = future_pmap_dbl(list(velocity, temperature, body_mass, body_length, ucrit), 
                                      ~ optimize(get_cost_of_travel, 
                                                 interval =  c(..1 + 1/1000, 10), 
                                                 water_velocity_m_per_s = ..1, 
                                                 fish_parm = fish_parm, body_mass_g = ..3, reach_temp_C = ..2, body_length_cm = ..4, swim_speed_ucrit = ..5, ratio = 1, fish_id)$minimum),
  ) 
toc()
full5 <- full4 %>% nest(data = -c(temperature, body_mass, body_length)) %>% 
  mutate(future_map_dfr(data, get_diff)) %>% 
  mutate(future_map_dfr(data, get_min_speed)) %>% 
  mutate(future_map_dfr(data, get_vel_at_max_speed)) 



full4  %>% group_by(temperature, body_length) %>% filter(signif(swim_speed,4) == signif(ucrit,4)) %>% 
  summarize(max = max(velocity), .groups = "drop")

full4  %>% 
  group_by(temperature, body_length) %>% 
  summarize(min_speed = min(swim_speed), .groups = "drop") %>% 
  inner_join(tibble_with_swim_speeds %>% ungroup(), by = "temperature") %>% 
  group_by(body_length) %>% 
  filter(min_speed == max(min_speed)) %>% 
  filter(ucrit == min(ucrit)) %>% 
  select(min_speed) %>% 
  distinct() %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(body_length, min_speed))

full4 %>% 
  group_by(temperature, body_length) %>% 
  filter(signif(swim_speed, 4) == signif(swim_speed_max, 4)) %>% 
  summarize(vel_at_max_speed = min(velocity), .groups = "drop") %>% 
  ggplot() +
  geom_line(aes(body_length, vel_at_max_speed))


full6 <- full5 %>% 
  # select(-data) %>% 
  mutate(
    across(.cols = -data, 
           .fns = ~ ifelse(is.infinite(.x), NA_real_, .x))) %>% 
  drop_na() %>% 
  unnest(data) %>% 
  nest(data = -body_length) %>% 
  mutate(ucrit_cutoff = map_dbl(data, get_min_speed_for_high_ucrit)) %>% 
  unnest(data)

# body length vs velocity at max_speed

full6 %>% 
  ggplot() +
  geom_line(aes(body_length, vel_at_max_speed))

bl_vs_max_speed <- full6 %>% 
  nest(data = everything()) %>% 
  mutate(fit = map(data, ~lm(vel_at_max_speed ~ body_length, data = .)),
         tidy = map(fit, broom::tidy)) %>% 
  select(tidy) %>% 
  unnest(tidy) %>% pull(estimate)

# body length vs ucrit cut-off; i.e., when ucrit isn't the lowest speed a fish swims at

full6 %>% 
  nest(data = -c(body_length, ucrit_cutoff)) %>%
  ggplot(aes(body_length, ucrit_cutoff)) +
  geom_line()

bl_vs_ucrit_cutoff <- full6 %>% 
  nest(data = everything()) %>% 
  mutate(fit = map(data, ~lm(ucrit_cutoff ~ body_length, data = .)),
         tidy = map(fit, broom::tidy)) %>% 
  select(tidy) %>% 
  unnest(tidy) %>% pull(estimate)

# min velocity at which ucrit is used when ucrit isn't the lowest speed a fish can swim

full6 %>% 
  filter(ucrit > ucrit_cutoff) %>%
  ggplot() +
  geom_line(aes(ucrit, min)) +
  facet_wrap(~body_length)


full6 %>% 
  filter(ucrit > ucrit_cutoff) %>% 
  nest(data = -body_length) %>% 
  mutate(fit = map(data, ~lm(min ~ ucrit, data = .)),
         tidy = map(fit, broom::tidy)) %>% 
  select(body_length , tidy) %>% 
  unnest(tidy)

# relationship is just ucrit + -ucrit_cutoff!

full6 %>% 
  filter(ucrit > ucrit_cutoff) %>% 
  nest(data = -body_length) %>% 
  mutate(fit = map(data, ~lm(min ~ ucrit, data = .)),
         tidy = map(fit, broom::tidy)) %>% 
  unnest(tidy) %>% select(body_length, estimate, term) %>% 
  pivot_wider(values_from = "estimate", names_from = "term") %>% 
  rename(intercept = `(Intercept)`) %>% 
  left_join(full6 %>% nest(data = -c(body_length, ucrit_cutoff))) %>% 
  select(-data) %>% 
  mutate(pct_diff = (intercept + ucrit_cutoff)/intercept * 100)

# max velocity at which ucrit is used before switching to burst speed

mod <- lm(full6$max ~ poly(full6$ucrit, 3))

pred <- predict(mod, full6)


full6 %>% 
  nest(data = everything()) %>% 
  mutate(fit = map(data, ~ lm(max ~ ucrit, data = .)),
         augment = map(fit, broom::augment)) %>% 
  select(augment) %>% 
  unnest(augment) %>% 
  mutate(pred = predict(mod, full6)) %>% 
  #nest(data = -c(ucrit, max)) %>% 
  ggplot() +
  geom_point(aes(ucrit, max)) +
  geom_line(aes(ucrit, .fitted), color = "red")


max_vel_for_ucrit <- full6 %>% 
  nest(data = everything()) %>% 
  mutate(fit = map(data, ~ lm(max ~ ucrit, data = .)),
         tidy = map(fit, broom::tidy)) %>% 
  select(tidy) %>% 
  unnest(tidy) %>% pull(estimate)

full6 %>% 
  # filter(ucrit > ucrit_cutoff) %>%
  ggplot() +
  geom_line(aes(ucrit, max^(1/3))) +
  facet_wrap(~body_length)

temp <- 20
length <- 100
ucrit_cutoff <- get_estimate_for_linear_model(bl_vs_ucrit_cutoff, length)
ucrit <- get_ucrit(fish_parm, temp, length, fish_id)
min_water_velocity_at_ucrit <- ucrit - ucrit_cutoff
max_water_velocity_at_ucrit <- get_estimate_for_linear_model(max_vel_for_ucrit, ucrit)
water_velocity_at_max_burst <- get_estimate_for_linear_model(bl_vs_max_speed, length)
swim_speed_max <- get_max_swim_speed(length)
cell_width <- 20

full7 <- full6 %>% 
  filter(temperature == temp, body_length == length) %>% 
  select(swim_speed, velocity) %>% 
  rename(swim_speed_martin = swim_speed) %>% 
  mutate(
    ratio = 1,
    #node_dir = "f"
  )

res <- get_cell_swim_speeds(full7, ucrit, swim_speed_max, cell_width) %>% 
  mutate(swim_type = fcase(swim_speed < ucrit, "sub-Ucrit swim speed",
                           swim_speed == ucrit, "Ucrit",
                           default = "burst swimming")) %>% 
  filter(swim_speed_martin < swim_speed_max)

ggplot(res) +
  geom_point(aes(velocity, swim_speed, color = swim_type))+
  geom_line(aes(velocity, swim_speed_martin)) +
  xlab("water velocity (m/s)")+
  ylab("optimal swim speed (m/s)") +
  # xlim(c(0,3)) +
  # ylim(c(ucrit_cutoff,4))+
  theme_classic()

cols <- c("num_paths.x", "num_paths.y")

merge(ch1, ch2, all = TRUE, by = c("distance", "lat_dist"))[, ':='(
    num_paths = num_paths.x + num_paths.y
    )][, cols := NULL][]


merge(ch1, ch2, all = TRUE, by = c("distance", "lat_dist"))[, ':='(
  num_paths = num_paths.x + num_paths.y,
  num_paths.x = NULL,
  num_paths.y = NULL
)]

adjust_velocity_for_wall_factor2 <- function(df, fish_parm, habitat_parm) {
  if (fish_parm$benthic_fish == 1) {
    df %>%
      .[depth >= habitat_parm$ben_vel_height, velocity := velocity * habitat_parm$base_wall_factor]
  } else {
    df
  }
}

full_df <- full %>% as_tibble()

microbenchmark(

full_df %>% select(-energy_cost) %>% unnest(paths)
,
full[, data.table::rbindlist(paths), by = .(species, number)]
)


function(df){
  dt[,energy_cost]
}