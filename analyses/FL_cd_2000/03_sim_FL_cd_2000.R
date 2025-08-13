###############################################################################
# Simulate plans for `FL_cd_2000`
# © ALARM Project, August 2025
###############################################################################

# Run the simulation -----
cli_process_start("Running simulations for {.pkg FL_cd_2000}")

# TODO any pre-computation (VRA targets, etc.)

# Global settings
cluster_tol <- .005
nsims_south <- 2000
nsims_north <- 2000
nsims <- 2000

map$row_num <- 1:nrow(map)

# South Florida

map_south <- map %>% filter(section == "South")

attr(map_south, "pop_bounds") <- attr(map, "pop_bounds")

map_south <- set_pop_tol(map_south, cluster_tol)

map <- map %>%
  mutate(cluster_edge = ifelse(row_num %in% map_south$row_num, 1, 0))

z <- geomander::seam_geom(map$adj, map, admin = "cluster_edge", seam = c(0, 1))

z <- z[z$cluster_edge == 1, ]

map_south$cluster_edge <- map_south$row_num %in% z$row_num

constraints <- redist_constr(map_south) %>%
  add_constr_grp_hinge(3, vap_black, vap, .45) %>%
  add_constr_grp_hinge(-5, vap_black, vap, .2)  %>%
  add_constr_grp_hinge(3, vap_hisp, vap, .6) %>%
  add_constr_grp_hinge(-5, vap_hisp, vap, .3) %>%
  add_constr_custom(
    strength = 10,
    fn = function(plan, distr) {
      as.numeric(!any(plan[map_south$cluster_edge] == 0))
    }
  )

n_steps <- (sum(map_south$pop)/attr(map_south, "pop_bounds")[2]) %>% floor()

plans_south <- redist_smc(map_south,
                          constraints = constraints,
                          nsims = 1e3,
                          runs = 4,
                          ncores = 32,
                          n_steps = n_steps,
                          sampling_space = "spanning_forest",
                          ms_params = list(ms_frequency = 1L, ms_moves_multiplier = 2),
                          split_params = list(splitting_schedule = "any_valid_sizes"),
                          counties = pseudo_county,
                          verbose = TRUE)

plans_south <- plans_south %>%
  mutate(ndv = group_frac(map_south, ndv, ndv + nrv),
         hvap = group_frac(map_south, vap_hisp, vap),
         bvap = group_frac(map_south, vap_black, vap),
         hcvap = group_frac(map_south, cvap_hisp, cvap),
         bcvap = group_frac(map_south, cvap_black, cvap),
         dem16 = group_frac(map_south, adv_16, arv_16 + adv_16),
         dem18 = group_frac(map_south, adv_18, arv_18 + adv_18),
         dem20 = group_frac(map_south, adv_20, arv_20 + adv_20))

plans_south <- plans_south %>% group_by(draw) %>%
  mutate(first = max(bvap), second = sort(bvap, decreasing = TRUE)[2]) %>%
  ungroup() %>% filter((first > 0.4 & second > .25) | draw == "cd_2010") %>%
  select(-c(first, second))

samp <- sample(seq_len(ncol(get_plans_matrix(plans_south))), 35000)

plans_south <- plans_south %>% group_by(draw) %>%
  mutate(id = cur_group_id()) %>%
  ungroup(draw) %>% filter(id %in% samp) %>% select(-c(id))

summary(plans_south)

###

# North Florida

map_north <- map %>% filter(section == "North")

attr(map_north, "pop_bounds") <- attr(map, "pop_bounds")

map_north <- set_pop_tol(map_north, cluster_tol)

map <- map %>%
  mutate(cluster_edge = ifelse(row_num %in% map_north$row_num, 1, 0))

z <- geomander::seam_geom(map$adj, map, admin = "cluster_edge", seam = c(0, 1))

z <- z[z$cluster_edge == 1, ]

map_north$cluster_edge <- map_north$row_num %in% z$row_num

constraints <- redist_constr(map_north) %>%
  add_constr_grp_hinge(6, cvap_black, cvap, .5) %>%
  add_constr_grp_hinge(-6, cvap_black, cvap, .2) %>%
  add_constr_grp_hinge(3, cvap_hisp, cvap, .7) %>%
  add_constr_grp_hinge(-6, cvap_hisp, cvap, .3) %>%
  add_constr_custom(
    strength = 10,
    fn = function(plan, distr) {
      as.numeric(!any(plan[map_north$cluster_edge] == 0))
    }
  )

n_steps <- (sum(map_north$pop)/attr(map, "pop_bounds")[2]) %>% floor()

plans_north <- redist_smc(map_north,
                          counties = pseudo_county,
                          nsims = nsims_north,
                          runs = 2L, ncores = 31L,
                          n_steps = n_steps,
                          constraints = constraints,
                          verbose = T)

plans_north <- plans_north %>%
  mutate(hvap = group_frac(map_north, vap_hisp, vap),
         bvap = group_frac(map_north, vap_black, vap),
         hcvap = group_frac(map_north, cvap_hisp, cvap),
         bcvap = group_frac(map_north, cvap_black, cvap),
         dem16 = group_frac(map_north, adv_16, arv_16 + adv_16),
         dem18 = group_frac(map_north, adv_18, arv_18 + adv_18),
         dem20 = group_frac(map_north, adv_20, arv_20 + adv_20))

plans_north <- plans_north %>% group_by(draw) %>%
  mutate(count = sum(bvap >= .25 & district != 0)) %>% ungroup() %>% filter(count > 0) %>%
  select(-c(count))

samp <- sample(seq_len(ncol(get_plans_matrix(plans_north))), 35000)

plans_north <- plans_north %>% group_by(draw) %>%
  mutate(id = cur_group_id()) %>%
  ungroup(draw) %>% filter(id %in% samp) %>% select(-c(id))

summary(plans_north)

# Merge north and south

plans_north <- plans_north %>% filter(draw != "cd_2010")
plans_south <- plans_south %>% filter(draw != "cd_2010")

plans_south$dist_keep <- ifelse(plans_south$district == 0, FALSE, TRUE)
plans_north$dist_keep <- ifelse(plans_north$district == 0, FALSE, TRUE)

prep_mat <- prep_particles(map = map,
                           map_plan_list = list(
                             list(map = map_south, plans = plans_south),
                             list(map = map_north, plans = plans_north)),
                           uid = row_num,
                           dist_keep = dist_keep,
                           nsims = nsims)

# Central Florida

constraints <- redist_constr(map) %>%
  add_constr_grp_hinge(
    12,
    cvap_hisp,
    total_pop = cvap,
    tgts_group = c(0.55)
  ) %>%
  add_constr_grp_hinge(
    12,
    cvap_black,
    total_pop = cvap,
    tgts_group = c(0.55)
  )

plans <- redist_smc(map, nsims = nsims, runs = 2L, ncores = 31L,
                    counties = pseudo_county,
                    init_particles = prep_mat, verbose = T)

plans <- plans %>% filter(draw != "cd_2010") %>%
  group_by(chain) %>%
  filter(as.integer(draw) < min(as.integer(draw)) + 2500) %>% # thin samples
  ungroup()

plans <- plans %>% add_reference(ref_plan = map$cd_2010)

# IF CORES OR OTHER UNITS HAVE BEEN MERGED:
# make sure to call `pullback()` on this plans object!
plans <- match_numbers(plans, "cd_2010")

cli_process_done()
cli_process_start("Saving {.cls redist_plans} object")

# Output the redist_map object. Do not edit this path.
write_rds(plans, here("data-out/FL_2010/FL_cd_2010_plans.rds"), compress = "xz")
cli_process_done()

# Compute summary statistics -----
cli_process_start("Computing summary statistics for {.pkg FL_cd_2010}")

plans <- add_summary_stats(plans, map) %>%
  mutate(total_cvap = tally_var(map, cvap), .after = total_vap)

cvap_cols <- names(map)[tidyselect::eval_select(starts_with("cvap_"), map)]
for (col in rev(cvap_cols)) {
  plans <- mutate(plans, {{ col }} := tally_var(map, map[[col]]), .after = vap_two)
}

# Output the summary statistics. Do not edit this path.
save_summary_stats(plans, "data-out/FL_2010/FL_cd_2010_stats.csv")

cli_process_done()



#######

### Gibbs constraints
constr <- redist_constr(map) %>%
  add_constr_grp_hinge(
    strength = 5,
    group_pop = vap_hisp,
    total_pop = vap,
  ) %>%
  add_constr_grp_hinge(
    strength = -6,
    group_pop = vap_hisp,
    total_pop = vap,
    tgts_group = .3
  ) %>%
  add_constr_grp_hinge(
    strength = 5,
    group_pop = vap_black,
    total_pop = vap,
  ) %>%
  add_constr_grp_hinge(
    strength = -6,
    group_pop = vap_black,
    total_pop = vap,
    tgts_group = .3
  )

### Hard constraints

map <- map |> mutate(vap_bh = vap_black + vap_hisp)

constr2 <- redist_constr(map) |>
  add_constr_custom_plan(1, function(plan, seats, num_regions) {
    
    mmds <- cbind(map, plan) |>
      filter(plan %in% seats[seats == 1]) |> 
      group_by(plan) |>
      summarise(mvap = sum(vap_bh)/sum(vap)) |>
      filter(mvap > .5) |> nrow()
    
    pot_mmds <- sum((cbind(map, plan) |>
      filter(plan %in% seats[seats > 1]) |> 
      group_by(plan) |>
      summarise(mvap = sum(vap_bh)/sum(vap)) |> pull(mvap)) * seats[seats > 1])
    
    print(glue::glue("MMDs: {mmds}, Potential: {pot_mmds}"))
    
    return(-1*(mmds + pot_mmds))
  }, thresh = -10, only_final_plans = FALSE)

  
# TODO customize as needed. Recommendations:
#  - For many districts / tighter population tolerances, try setting
#  `pop_temper=0.01` and nudging upward from there. Monitor the output for
#  efficiency!
#  - Monitor the output (i.e. leave `verbose=TRUE`) to ensure things aren't breaking
#  - Don't change the number of simulations unless you have a good reason
#  - If the sampler freezes, try turning off the county split constraint to see
#  if that's the problem.
#  - Ask for help!

set.seed(2000)

plans <- redist_smc(map,
                    constraints = constr,
                    nsims = 2e3,
                    runs = 4,
                    ncores = 32,
                    sampling_space = redist:::FOREST_SPACE_SAMPLING,
                    ms_params = list(ms_frequency = 1L, ms_moves_multiplier = 32),
                    split_params = list(splitting_schedule = "any_valid_sizes"),
                    counties = pseudo_county,
                    verbose = TRUE)

# IF CORES OR OTHER UNITS HAVE BEEN MERGED:
# make sure to call `pullback()` on this plans object!

plans <- plans %>%
    group_by(chain) %>%
    filter(as.integer(draw) < min(as.integer(draw)) + 1000) %>% # thin samples
    ungroup()
plans <- match_numbers(plans, "cd_2000")

cli_process_done()
cli_process_start("Saving {.cls redist_plans} object")

# TODO add any reference plans that aren't already included

# Output the redist_map object. Do not edit this path.
write_rds(plans, here("data-out/FL_2000/FL_cd_2000_plans.rds"), compress = "xz")
cli_process_done()

# Compute summary statistics -----
cli_process_start("Computing summary statistics for {.pkg FL_cd_2000}")

plans <- add_summary_stats(plans, map)

# Output the summary statistics. Do not edit this path.
save_summary_stats(plans, "data-out/FL_2000/FL_cd_2000_stats.csv")

cli_process_done()

# Extra validation plots for custom constraints -----
# TODO remove this section if no custom constraints
if (interactive()) {
    library(ggplot2)
    library(patchwork)

    validate_analysis(plans, map |> mutate(state = "FL"))
    summary(plans)
    
    ## VAP charts
    d1 <- redist.plot.distr_qtys(
      plans,
      vap_black/total_vap,
      color_thresh = NULL,
      color = ifelse(
        subset_sampled(plans)$ndv > subset_sampled(plans)$nrv,
        "#3D77BB",
        "#B25D4C"
      ),
      size = 0.5,
      alpha = 0.5
    ) +
      scale_y_continuous("Percent Black by VAP") +
      labs(title = "FL Proposed Plan versus Simulations") +
      scale_color_manual(values = c(cd_2000 = "black"))
    
    d2 <- redist.plot.distr_qtys(
      plans,
      vap_hisp/total_vap,
      color_thresh = NULL,
      color = ifelse(
        subset_sampled(plans)$ndv > subset_sampled(plans)$nrv,
        "#3D77BB",
        "#B25D4C"
      ),
      size = 0.5,
      alpha = 0.5
    ) +
      scale_y_continuous("Percent Hispanic by VAP") +
      labs(title = "FL Proposed Plan versus Simulations") +
      scale_color_manual(values = c(cd_2000 = "black"))
    
    d3 <-
      redist.plot.distr_qtys(
        plans,
        (vap_hisp + vap_black)/total_vap,
        color_thresh = NULL,
        color = ifelse(
          subset_sampled(plans)$ndv > subset_sampled(plans)$nrv,
          "#3D77BB",
          "#B25D4C"
        ),
        size = 0.5,
        alpha = 0.5
      ) +
      scale_y_continuous("HVAP + BVAP / VAP") +
      labs(title = "FL Proposed Plan versus Simulations") +
      scale_color_manual(values = c(cd_2000 = "black"))
    
}
