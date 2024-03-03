library(readr)    # For reading CSV files
library(dplyr)    # For data manipulation
library(ggplot2)  # For plotting
library(MASS)     # For robust linear regression (rlm function)
library(deSolve)  # For solving differntial equations
library(EpiModel)


################################################################################
# Case-to-Hospitalization regression.
################################################################################

# Read in the data
cases_by_day <- read_csv("COVID_pub_data/Bronx_PH_data/cases-by-day.csv")
hosp_by_day <- read_csv("COVID_pub_data/Bronx_PH_data/hosp-by-day.csv")

# Define start and end dates for start of ALpha and Iota overlap
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2021-06-01")

# Subset Alpha & Iota data for the first period (Jan 1, 2021 - Jun 1, 2021) for cases and hospitalizations
cases_by_day_A <- cases_by_day[cases_by_day$date_of_interest >= start_date & cases_by_day$date_of_interest <= end_date, ]
cases_by_day_A <- cases_by_day_A[,c("date_of_interest", "BX_CASE_COUNT")]
hosp_by_day_A <- hosp_by_day[hosp_by_day$date_of_interest >= start_date & hosp_by_day$date_of_interest <= end_date, ]
hosp_by_day_A <- hosp_by_day_A[,c("date_of_interest", "BX_HOSPITALIZED_COUNT")]

# Define start and end dates for the first wave of Delta (Jun 2, 2021 - Oct 31, 2021)
start_date <- as.Date("2021-06-02")
end_date <- as.Date("2021-10-31")

# Subset Delta data for the second period for cases and hospitalizations
cases_by_day_D <- cases_by_day[cases_by_day$date_of_interest >= start_date & cases_by_day$date_of_interest <= end_date, ]
cases_by_day_D <- cases_by_day_D[,c("date_of_interest", "BX_CASE_COUNT")]
hosp_by_day_D <- hosp_by_day[hosp_by_day$date_of_interest >= start_date & hosp_by_day$date_of_interest <= end_date, ]
hosp_by_day_D <- hosp_by_day_D[,c("date_of_interest", "BX_HOSPITALIZED_COUNT")]

# Define start and end dates for the first Omicron wave (Dec 6, 2021 - Feb 28, 2022)
start_date <- as.Date("2021-12-06")
end_date <- as.Date("2022-2-28")

# Subset Omicron data for the third period for cases and hospitalizations
cases_by_day_O <- cases_by_day[cases_by_day$date_of_interest >= start_date & cases_by_day$date_of_interest <= end_date, ]
cases_by_day_O <- cases_by_day_O[,c("date_of_interest", "BX_CASE_COUNT")]
hosp_by_day_O <- hosp_by_day[hosp_by_day$date_of_interest >= start_date & hosp_by_day$date_of_interest <= end_date, ]
hosp_by_day_O <- hosp_by_day_O[,c("date_of_interest", "BX_HOSPITALIZED_COUNT")]


# Create data frames for different virus strains and their respective hospitalizations and cases
Ancestral_epi <- data.frame(virus = "Iota and Alpha", 
                            hospitalizations = hosp_by_day_A$BX_HOSPITALIZED_COUNT,  
                            cases = cases_by_day_A$BX_CASE_COUNT)
Delta_epi <- data.frame(virus = "Delta",
                        hospitalizations = hosp_by_day_D$BX_HOSPITALIZED_COUNT,  
                        cases = cases_by_day_D$BX_CASE_COUNT)
Omicron_epi <- data.frame(virus ="Omicron",
                          hospitalizations = hosp_by_day_O$BX_HOSPITALIZED_COUNT,  
                          cases = cases_by_day_O$BX_CASE_COUNT)

# Combine data frames
epi_df <- bind_rows(Ancestral_epi, Delta_epi)
epi_df2 <- bind_rows(Delta_epi, Omicron_epi)
epi_df3 <- bind_rows(Ancestral_epi, Delta_epi, Omicron_epi)

# Generate histograms for each variable in epi_df except the first one (virus)
for (col in colnames(epi_df)[-1]) {
  hist(epi_df[[col]], main = col, xlab = col, col = "lightblue", border = "black")
}

# Plot hospitalizations against cases for Delta and Omicron periods separately
plot(hosp_by_day_D$BX_HOSPITALIZED_COUNT, cases_by_day_D$BX_CASE_COUNT)
plot(hosp_by_day_A$BX_HOSPITALIZED_COUNT, cases_by_day_A$BX_CASE_COUNT)

# Fit robust linear models for hospitalizations against cases for each period
Delta_model <- MASS::rlm(hosp_by_day_D$BX_HOSPITALIZED_COUNT ~ cases_by_day_D$BX_CASE_COUNT)
Ancestral_model <- MASS::rlm(hosp_by_day_A$BX_HOSPITALIZED_COUNT ~ cases_by_day_A$BX_CASE_COUNT)
Omicron_model <- MASS::rlm(hosp_by_day_O$BX_HOSPITALIZED_COUNT ~ cases_by_day_O$BX_CASE_COUNT)

# Fit combined model using robust linear regression
model_combined1 <- MASS::rlm(hospitalizations ~ cases + virus*cases , data = epi_df)
summary(model_combined1)
model_combined2 <- MASS::rlm(hospitalizations ~ cases + virus*cases , data = epi_df2)
summary(model_combined2)

# Calculate p-values for differences in slopes during different virus dominance periods
p_value_delta_iota <- 2 * (1 - pt(abs(-0.3871), df = 300))
p_value_delta_omicron <- 2 * (1 - pt(abs(-7.2216), df = 233))

# Plot data with regression lines
ggplot(data = epi_df3, aes(x = cases, y = hospitalizations, colour = virus)) +
  geom_point() +
  geom_smooth(aes(x = cases, y = hospitalizations, group = virus), method = MASS::rlm, se = F, fullrange=TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Iota and Alpha" = "gray", "Delta" = "green", "Omicron" = "red", "Linear Fit" = "black")) +
  labs(x = "New Cases Daily", y = "New Hospitalizations Daily", title = "Relationship between hospitalizations and cases") +
  theme_minimal()

################################################################################
# Process Inferred Case Data.
################################################################################
# Function to process virus data and plot incidence over time
process_and_plot_virus_data <- function(virus_data, virus_name, start_date, end_date, zero_fill_indices) {
  # Plot incidence over time for the virus
  plot(as.incidence(virus_data, dates = virus_data$Date))
  
  # Subset data for the virus within specified date range
  virus_data <- virus_data[virus_data$Date >= start_date & virus_data$Date <= end_date, ]
  
  # Calculate days since start and fill in missing values for the virus
  virus_data$days_since_start <- as.numeric(difftime(virus_data$Date, min(virus_data$Date), units = "days"))
  virus_data[virus_name][virus_data[virus_name] == 0] <- NA
  virus_filled <- fill(virus_data, virus_data)
  virus_filled[virus_name][zero_fill_indices] <- 0
  virus_filled$frac <- virus_filled[virus_name] / 1000
  
  return(virus_filled)
}

# Read and rename columns in the data
Estimated_VOC_cases <- read_csv("COVID_pub_data/Estimated_C19_VOC_cases.csv", col_names = c("Date","Others", "Iota", "Omicron", "Delta", "Alpha")) 

# Define date range for Iota, Alpha, and Delta
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2021-10-31")

# Process data for Iota, Alpha, and Delta
viruses <- c("Iota", "Alpha", "Delta")
processed_data <- list()

for (virus in viruses) {
  zero_fill_indices <- ifelse(virus == "Iota", 166:304, 162:304)  # Adjust zero fill indices
  processed_data[[virus]] <- process_and_plot_virus_data(Estimated_VOC_cases, virus, start_date, end_date, zero_fill_indices)
}

# Define date range for Omicron
start_date_omicron <- as.Date("2021-12-06")
end_date_omicron <- as.Date("2022-03-08")

# Process data for Omicron
Omicron_filled <- process_and_plot_virus_data(Estimated_VOC_cases, "Omicron", start_date_omicron, end_date_omicron, NA)

# Create zoo objects for smoothed data
Iota.zoo <- zoo(processed_data$Iota$Iota, order.by = processed_data$Iota$Date)
Alpha.zoo <- zoo(processed_data$Alpha$Alpha, order.by = processed_data$Alpha$Date)

# Calculate rolling averages and exponential moving averages
ema_I <- as.vector(ksmooth(processed_data$Iota$days_since_start, processed_data$Iota$Iota, 'normal', bandwidth = 30))
ema_A <- as.vector(ksmooth(processed_data$Alpha$days_since_start, processed_data$Alpha$Alpha, 'normal' , bandwidth = 30))
ema_D <- as.vector(ksmooth(processed_data$Delta$days_since_start, processed_data$Delta$Delta, 'normal' , bandwidth = 30))

# Create data frames for plotting
Alpha.iota.delta <- data.frame(
  Date = processed_data$Iota$Date,
  Days = processed_data$Iota$days_since_start,
  ema_A = ema_A$y,
  ema_I = ema_I$y,
  ema_D = ema_D$y
)

# Plot the original data and the rolling average
ggplot(Alpha.iota.delta, aes(x = Days)) +
  geom_line(aes(y = ema_I), color = "yellow", size = 1.5, linetype = "solid") +
  geom_line(aes(y = ema_A), color = "purple", size = 1.5, linetype = "solid") +
  geom_line(aes(y = ema_D), color = "green", size = 1.5, linetype = "solid") +
  geom_point(data = processed_data$Iota, aes(x = days_since_start, y = Iota), color = "yellow", alpha = 0.3) +
  geom_point(data = processed_data$Alpha, aes(x = days_since_start, y = Alpha), color = "purple", alpha = 0.3) +
  geom_point(data = processed_data$Delta, aes(x = days_since_start, y = Delta), color = "green", alpha = 0.3) +
  labs(
    title = "Smoothed Alpha, Iota and Delta Cases",
    x = "Time (days)",
    y = "Cases per 100,000 individuals"
  ) +
  theme_minimal()

################################################################################
# Infective Force regression.
################################################################################

# Regression to find infective force
FOI_I <- lm(X_virusI[1:51] ~ R_virusI[1:51])
FOI_A <- lm(X_virusA[8:89] ~ R_virusA[8:89])
FOI_D <- lm(X_virusD[140:229] ~ R_virusD[140:229])

# Summary of regressions
summary(FOI_I)
summary(FOI_A)
summary(FOI_D)

# Create data frames for each model
df_FOI_I <- data.frame(
  model = "FOI_I",
  X = X_virusI[1:51],
  R = R_virusI[1:51]
)

df_FOI_A <- data.frame(
  model = "FOI_A",
  X = X_virusA[8:89],
  R = R_virusA[8:89]
)

df_FOI_D <- data.frame(
  model = "FOI_D",
  X = X_virusD[140:229],
  R = R_virusD[140:229]
)

# Combine data frames into a single data frame
data_df <- bind_rows(df_FOI_I, df_FOI_A, df_FOI_D)

# Plot the data with regression lines
ggplot(data_df, aes(x = R, y = X, color = model)) +
  geom_point() +
  geom_smooth(aes(color = model), method = "lm", se = FALSE) +
  scale_color_manual(
    values = c("FOI_I" = "yellow", "FOI_A" = "purple", "FOI_D" = "green"),
    labels = c("FOI_I" = "Iota", "FOI_A" = "Alpha", "FOI_D" = "Delta")
  ) +
  labs(
    x = "Cumulative Cases Daily",
    y = "New Cases Daily",
    title = "Force of Infection Estimation"
  ) +
  theme_minimal()

################################################################################
# Doubling time regression.
################################################################################

# Define a function to calculate doubling time
calculate_doubling_time <- function(model, start_index, end_index) {
  # Fit the model
  summary_model <- summary(model)
  
  # Extract the growth rate (slope) from the model
  growth_rate <- coef(summary_model)["Days"][2]  # Assuming "Days" is the predictor variable
  
  # Calculate doubling time
  doubling_time <- log(2) / growth_rate
  
  return(doubling_time)
}

# Calculate doubling time for each model
doubling_time_I <- calculate_doubling_time(modelI, 1, 51)
doubling_time_A <- calculate_doubling_time(modelA, 8, 89)
doubling_time_D <- calculate_doubling_time(modelD, 140, 229)

# Print doubling times
doubling_time_I
doubling_time_A
doubling_time_D


################################################################################
# SIR model paramter estimation (dynamic)
################################################################################

mod_SIR_1g_cl <- function(t, initial_conditions, parms) {
  with(as.list(c(initial_conditions, parms)), {
    
    # Derivations
    num <- s.num + i.num + r.num
    
    # Parameters
    lambda <- inf.prob * act.rate * i.num / num
    
    # Flows
    si.flow <- lambda * s.num
    ir.flow <- rec.rate * i.num
    
    # ODEs
    dS <- -si.flow
    dI <- si.flow - ir.flow
    dR <- ir.flow
    
    # Output
    list(c(dS, dI, dR))
  })
}

# Extract data for each virus
True_I1 <- Alpha.iota$ema_I
True_I2 <- Alpha.iota$ema_A
True_I3 <- Alpha.iota$ema_D[155:304]


initial_conditions_iota <- c(
  s.num = 100000 - 120 - 5000,      # Initial susceptible fraction
  i.num = 120,      # Initial infected fraction for virus (estimated from actual data)
  r.num = 500      # Initial recovered fraction (estimated from actual data)
)

initial_conditions_alpha <- c(
  s.num = 100000-6-5000,      # Initial susceptible fraction
  i.num = 6,      # Initial infected fraction for virus (estimated from actual data)
  r.num = 5000      # Initial recovered fraction (estimated from actual data)
)

initial_conditions_delta <- c(
  s.num = 100000-1-25000,      # Initial susceptible fraction
  i.num = 1,      # Initial infected fraction for virus (estimated from actual data)
  r.num = 25000      # Initial recovered fraction (inferred from previous model runs)
)

initial_parameters <- c(
  inf.prob = 0.3082492,
  act.rate = 2.9335948,        # Transmission rate for virus 
  rec.rate = 0.6264471       # Recovery rate for virus 
)

times1 = seq(0, 303, by = 1)
times2 = seq(0, 303, by = 1)
times3 = seq(0, 149, by = 1)

# Define a function to calculate RSS
calculate_RSS <- function(parameters, True_I, times, initial_conditions) {
  parms <- setNames(parameters, c("inf.prob", "act.rate", "rec.rate"))
  solution <- ode(y = initial_conditions, times = times, func = mod_SIR_1g_cl, parms = parms)
  model_I <- solution[, "i.num"]
  RSS_I <- sum((True_I - model_I)^2)
  return(RSS_I)
}

# Perform optimization using SANN method
optimized_params_iota <- optim(par = initial_parameters, fn = calculate_RSS, initial_conditions = initial_conditions_iota,True_I = True_I1, times = times1, method = "Nelder-Mead")

optimized_params_iota  <- optimized_params_iota$par

optimized_params_iota 


optimized_params_alpha <- optim(par = initial_parameters, fn = calculate_RSS, initial_conditions = initial_conditions_alpha,True_I = True_I2, times = times2, method = "Nelder-Mead")

optimized_params_alpha   <- optimized_params_alpha$par

optimized_params_alpha 


optimized_params_delta <- optim(par = initial_parameters, fn = calculate_RSS, initial_conditions = initial_conditions_delta,True_I = True_I3, times = times3, method = "Nelder-Mead")

optimized_params_delta  <- optimized_params_delta$par

optimized_params_delta 

################################################################################
# SIR model visualization (dynamic)
################################################################################

# Define function to run SIR model
run_SIR_model <- function(param, s.num, i.num, r.num, nsteps) {
  init <- init.dcm(s.num = s.num, i.num = i.num, r.num = r.num)
  control <- control.dcm(type = "SIR", nsteps = nsteps, dt = 1)
  mod <- dcm(param, init, control)
  return(mod)
}


# Run SIR model for each virus
mod1 <- run_SIR_model(param = param.dcm(inf.prob = 0.1511724, act.rate = 3.0921119, rec.rate = 0.4173837),
                      s.num = 100000 - 120 - 5000, i.num = 120, r.num = 5000, nsteps = 304)
mod2 <- run_SIR_model(param = param.dcm(inf.prob = 0.3403329, act.rate = 2.8619365, rec.rate = 0.8789471),
                      s.num = 100000 - 6 - 5000, i.num = 6, r.num = 5000, nsteps = 304)
mod3 <- run_SIR_model(param = param.dcm(inf.prob = 0.4975271, act.rate = 2.5700256, rec.rate = 0.8793897),
                      s.num = 100000 - 1 - 25000, i.num = 1, r.num = 25000, nsteps = 150)


# Extract data from the models
mod1_df <- as.data.frame(mod1)
mod2_df <- as.data.frame(mod2)
mod3_df <- as.data.frame(mod3)

# Adjust mod3 to be on same time scale
output3_zero = data.frame(matrix(0, nrow = 154 , ncol = 7))
colnames(output3_zero) = colnames(mod3_df)
mod3_df = rbind(output3_zero, mod3_df)
mod3_df$time = 1:304

# Combine the data frames
combined_output <- rbind(mod1_df, mod2_df, mod3_df)

# Create a scenario column
combined_output$Scenario <- rep(c("Iota", "Alpha", "Delta"), each = nrow(mod1_df))

# Plot the results
ggplot(combined_output) +
  geom_line(aes(x = time, y = i.num, col = Scenario), size = 1.5) +
  labs(title = "Infected Fraction under SIR Model",
       x = "Time (days)",
       y = "Fraction of Infected Population") +
  scale_color_manual(values = c("Iota" = "yellow", "Alpha" = "purple", "Delta" = "green")) +
  geom_point(data = Iota_filled, aes(x = days_since_start, y = Iota), color = "yellow", alpha = 0.3) +
  geom_point(data = Alpha_filled, aes(x = days_since_start, y = Alpha), color = "purple", alpha = 0.3) +
  geom_point(data = Delta_filled, aes(x = days_since_start, y = Delta), color = "green", alpha = 0.3) +
  theme_minimal()


################################################################################
# SIR model paramter estimation (stochastic)
################################################################################

# Define the function to run Euler's method
eulmar <- function(func, xstart, times, dt, ...) {
  out <- array(
    dim = c(length(times), length(xstart)),
    dimnames = list(NULL, names(xstart))
  )
  out[1,] <- x <- xstart
  t <- times[1]
  for (k in seq.int(from = 2, to = length(times))) {
    while (t < times[k]) {
      dx <- func(t, x, dt, ...)
      x <- x + dx
      t <- t + dt
    }
    out[k,] <- x
  }
  as.data.frame(cbind(time = times, out))
}

# Define the SIR model
sir.eulerstep <- function(t, x, dt, Beta, gamma) {
  N <- sum(x)
  means <- c(
    Beta * x[1] * x[2] / N,
    gamma * x[2]
  ) * dt
  dn <- rnorm(n = 2, mean = means, sd(means))
  c(-dn[1], dn[1] - dn[2], dn[2])
}

# Define the objective function
objective_function <- function(params, xstart, times) {
  Beta <- params[1]
  gamma <- params[2]
  
  x <- eulmar(
    func = sir.eulerstep,
    xstart = xstart, Beta = Beta, gamma = gamma,
    times = times, dt = 1
  )
  
  # Assuming data contains the observed values
  observed_values <- True_I3
  predicted_values <- x[, "Y"]
  
  # Use sum of squared differences as the objective to minimize
  sum((observed_values - predicted_values)^2)
}

# Initial guess for parameters
initial_params <- c(Beta = 1.060927, gamma = 0.9501149)

# Initial guesses for each scenario
initial_guesses <- list(
  scenario1 = initial_params,
  scenario2 = initial_params,
  scenario3 = initial_params
)

# Define the scenarios
scenarios <- list(
  scenario1 = list(
    xstart = c(X = 100000 - 1 - 25000, Y = 1, Z = 25000),
    times = seq(0, 149, by = 1)
  ),
  scenario2 = list(
    xstart = c(X = 100000 - 40, Y = 40, Z = 0),
    times = seq(0, 304, by = 1)
  ),
  scenario3 = list(
    xstart = c(X = 100000 - 10000 - 10, Y = 10, Z = 10000),
    times = seq(0, 304, by = 1)
  )
)

# Optimize parameters for each scenario
optimized_params <- list()
for (scenario in names(scenarios)) {
  optimized_params[[scenario]] <- optim(
    par = initial_guesses[[scenario]], 
    fn = objective_function, 
    xstart = scenarios[[scenario]]$xstart, 
    times = scenarios[[scenario]]$times,
    method = "Nelder-Mead"
  )$par
}

# Display the optimized parameters
for (scenario in names(optimized_params)) {
  cat("Optimized parameters for", scenario, ":", optimized_params[[scenario]], "\n")
}


################################################################################
# SIR model visualization (stochastic)
################################################################################

# Parameters and initial conditions
times <- seq(0, 304, by = 1)
times3 <- seq(0, 150, by = 1)

# Define sets of parameters and initial conditions
parameter_sets <- list(
  set1 = c(Beta = 0.7503133, gamma = 0.6816677),
  set2 = c(Beta = 1.187811, gamma = 1.120855),
  set3 = c(Beta = 0.9754699, gamma = 0.8244839)
)

initial_conditions <- list(
  initial_condition1 = c(X = 99960, Y = 40, Z = 0),
  initial_condition2 = c(X = 99999, Y = 1, Z = 0),
  initial_condition3 = c(X = 89990, Y = 10, Z = 10000)
)

# Number of simulations
num_simulations <- 1000

# Create a list to store mean and std. dev data frames for each parameter set
mean_data_frames <- list()

# Outer loop for different parameter sets (n values from 1 to 3)
for (n in 1:3) {
  current_set_data_frames <- list()
  times_used <- if (n <= 2) times else times3
  for (i in 1:num_simulations) {
    common_params <- parameter_sets[[n]]
    xstart <- initial_conditions[[n]]
    simulation_result <- sir_simulation(Beta = common_params[["Beta"]], gamma = common_params[["gamma"]], xstart = xstart, times = times_used)
    df_name <- paste0("simulation_", i)
    simulation_df <- data.frame(time = simulation_result$time, Y = simulation_result$Y)
    current_set_data_frames[[df_name]] <- simulation_df
  }
  mean_data <- data.frame(
    time = unique(unlist(current_set_data_frames[[1]]$time)),
    mean_Y = apply(do.call(cbind, lapply(current_set_data_frames, function(df) df$Y)), 1, mean),
    std_dev_Y = apply(do.call(cbind, lapply(current_set_data_frames, function(df) df$Y)), 1, sd)
  )
  mean_data_frames[[paste0("parameter_set_", n)]] <- mean_data
}

# Extract mean data frames
Iota_Stoch <- mean_data_frames[[1]]
Alpha_Stoch <- mean_data_frames[[2]]
Delta_stoch_short <- mean_data_frames[[3]]

# Regression analysis
rlm_model <- function(data, predictor) {
  MASS::rlm(data[["mean_Y"]] ~ predictor)
}

lm_model <- function(data, predictor) {
  summary(lm(data[["mean_Y"]] ~ predictor))
}

# Perform regression analysis
rlm_Iota <- rlm_model(Iota_Stoch[1:305, ], Alpha.iota.delta$ema_I)
lm_Iota <- lm_model(Iota_Stoch[1:305, ], Alpha.iota$ema_I)

rlm_Alpha <- rlm_model(Alpha_Stoch[1:305, ], Alpha.iota.delta$ema_A)
lm_Alpha <- lm_model(Alpha_Stoch[1:305, ], Alpha.iota.delta$ema_A)

rlm_Delta <- rlm_model(Delta_stoch_short, Alpha.iota.delta$ema_D)
lm_Delta <- lm_model(Delta_stoch_short, Alpha.iota.delta$ema_D)

# Create zero-filled data frame
output3_zero <- data.frame(matrix(0, nrow = 154 , ncol = 3))
colnames(output3_zero) <- colnames(Delta_stoch_short)
Delta_stoch <- rbind(output3_zero, Delta_stoch_short)
Delta_stoch$time <- 0:304

# Combine the output data from both simulations
combined_output_stoch <- rbind(Iota_Stoch, Alpha_Stoch, Delta_stoch)

# Create a data frame from the combined output
combined_output_stoch <- as.data.frame(combined_output_stoch) %>%
  mutate(Scenario = rep(c("Iota", "Alpha", "Delta"), each = nrow(Iota_Stoch)))

# Plot the results using ggplot2
ggplot(combined_output_stoch, aes(x = time, y = mean_Y, color = Scenario, fill = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_Y - std_dev_Y, ymax = mean_Y + std_dev_Y), alpha = 0.2) +
  scale_color_manual(values = c(Iota = "yellow", Alpha = "purple", Delta = "green")) +
  scale_fill_manual(values = c(Iota = "yellow", Alpha = "purple", Delta = "green")) +
  labs(title = "Stochastic SIR Model",
       x = "Time",
       y = "Mean cases per 100,000 individuals") +
  theme_minimal()





