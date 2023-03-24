#!/usr/bin/env Rscript
library(magrittr)

# 1. read in data ----------------------------------------------------------

data <- sf::read_sf(here::here("data/ds2890.gdb")) %>% 
    dplyr::rename(geometry = Shape)

# 2. functions ------------------------------------------------------------

# extracts numbers from a character string

extract_nums <- function(string) {
  stringr::str_extract_all(string, "\\(?[0-9,.]+\\)?")
}

# calculates the mean height depending on the numbers extracted from a string;
# e.g., "1 - 5m" would return the mean of 1 and 5

calc_height <- function(nums_as_string) {
  nums <- nums_as_string %>%
    as.numeric()
  if (length(nums) > 1) {
    return(mean(nums))
  } else if (nums <= 1) {
    return(nums / 2)
  } else {
    return(nums)
  }
}

# uses map_dbl to return all height values from a list of strings

calc_all_heights <- function(height_list) {
  height_list <- tidyr::replace_na(height_list, "0")
  nums_as_strings <- extract_nums(height_list)
  purrr::map_dbl(.x = nums_as_strings, .f = calc_height)
}

# joins a vector of strings using the "|" character

make_word_filter <- function(list) {
  paste(list, collapse = "|")
}

save_shape_file <- function(file, file_name, output_path) {
  sf::st_write(file, here::here(output_path, paste0(file_name, ".shp")), delete_layer = TRUE)
}

# 3. analysis -------------------------------------------------------------

# list of words related to forests for filtering relevant polygons

wood_words <- c(
  "wood",
  "forest"
)

# bundles the list into a single character separated by "|" (i.e., "OR")

wood_filter <- make_word_filter(wood_words)

gravel_sand_words <- c(
  "sand",
  "gravel"
)

gravel_sand_filter <- make_word_filter(gravel_sand_words)

rock_words <- c(
  "rock",
  "boulder",
  "cliff"
)

rock_filter <- make_word_filter(rock_words)

new_base_file <- data %>%
  # make columns lowercase so they are more consistently named
  dplyr::rename_with(tolower) %>%
  dplyr::filter(!grepl("mines", cv_group)) %>% 
  dplyr::mutate(
    # make the CV_Group column values lowercase for easier filtering
    cv_group = tolower(cv_group),
    # make percents into proportions
    dplyr::across(
      .cols = c(per_hardwo:per_shrub, per_total),
      .fns = ~ tidyr::replace_na(.x, 0) / 100
    ),
    # some values of cover were above 100% for some reason but should be 0
    dplyr::across(
      .cols = c(per_hardwo:per_shrub, per_total),
      # can use dplyr's vectorized if_else func b/c both T and F are the same data type
      .fns = ~ dplyr::if_else(.x >= 1, 0, .x)
    ),
    # rescaling percent cover so that it's 0 to 1
    dplyr::across(
      .cols = c(per_hardwo:per_shrub, per_total),
      .fns = ~ .x / max(per_total)
    ),
    # height was a string, so this pulls out useful values as floats
    height = calc_all_heights(ht_code),
    veg = per_total - per_tree,
    wood = per_tree,
    # create new class names for the polygons
    class = dplyr::case_when(
      grepl(wood_filter, cv_group) ~ "t_wood",
      grepl(gravel_sand_filter, cv_group) ~ "gravel",
      grepl(rock_filter, cv_group) ~ "rock",
      grepl("urban", cv_group) ~ "urban",
      grepl("wat: water", cv_group) ~ "water",
      TRUE ~ "t_veg"
    ),
    # add substrate values; most will be NA as they are unknown
    # can (and should) be updated by the user
    # some substrate data is included to some degree in the shape file
    # using "NA_real_" so that dplyr's if_else will return a number or an NA
    # note that gravel and sand are pooled in the original dataset, so their coverage gets split in half
    fine = dplyr::if_else(class == "gravel", 0.5, NA_real_),
    gravel = dplyr::if_else(class == "gravel", 0.5, NA_real_),
    cobble = NA_real_,
    rock = dplyr::if_else(class == "rock", 1, NA_real_)
  ) %>%
  dplyr::select(class, height:rock) %>%
  dplyr::mutate(dplyr::across(
    .cols = height:rock,
    ~ dplyr::if_else(class == "water" | class == "urban", NA_real_, .x)
  ))

# separate out the base file into shapefiles for ground cover and canopy cover
ground_cover <- new_base_file %>%
  dplyr::select(-height)

# 4. save output ----------------------------------------------------------

output_path <- here::here("data")
filename <- "cover"
save_shape_file(ground_cover, filename, output_path)