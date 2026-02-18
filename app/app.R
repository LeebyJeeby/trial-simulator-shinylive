library(shiny)
library(dplyr)
library(ggplot2)
library(bslib)
library(tidyr)
library(purrr)

# --- INSTRUCTIONS FOR ANTIGRAVITY ---
# 1. This app simulates 3-level educational trials (School -> Year -> Pupil).
# 2. It supports three randomization types: 
#    - Whole School (Cluster RCT)
#    - Year Group (Block randomized by School)
#    - Individual (Within-Year/School)
# 3. Data generation uses hierarchical sampling with user-defined RANGES for 
#    Year Groups and Class Sizes.
# 4. NOTE: Estimators restricted to OLS only for Shinylive compatibility 
#    (removed lme4/Matrix dependencies).
# ------------------------------------

# --- THE SIMULATION ENGINE ---
run_sim <- function(n_schools, 
                    years_range,    
                    classes_range,  
                    class_size_range, 
                    attrition_sch, 
                    attrition_pup, 
                    effect_size, 
                    icc_school,     
                    icc_year,       
                    estimators, 
                    design_type) {
  
  # --- 1. Variance Components ---
  total_cluster_var <- icc_school + icc_year
  if (total_cluster_var >= 1) {
    stop("Error: Combined ICC (School + Year) must be less than 1.")
  }
  
  sd_school <- sqrt(icc_school)
  sd_year   <- sqrt(icc_year)
  sd_resid  <- sqrt(1 - total_cluster_var)
  
  # --- 2. Data Generation (Hierarchical) ---
  
  # A. Generate Schools
  schools <- data.frame(
    school_id = 1:n_schools,
    u_school = rnorm(n_schools, 0, sd_school)
  )
  
  # Assign Treatment: Whole School Design
  if(design_type == "Whole School") {
    schools$treat_school <- rep(c(0, 1), length.out = n_schools)
  } else {
    schools$treat_school <- NA
  }
  
  # B. Generate Year Groups within Schools
  get_n <- function(range) { 
    if(range[1] == range[2]) return(range[1])
    sample(range[1]:range[2], 1) 
  }
  
  data_list <- list()
  
  for(i in 1:nrow(schools)) {
    s_id <- schools$school_id[i]
    s_u  <- schools$u_school[i]
    s_tr <- schools$treat_school[i]
    
    n_years <- get_n(years_range)
    
    # Create Year Groups
    for(y in 1:n_years) {
      y_id <- paste0(s_id, "_", y) 
      u_year <- rnorm(1, 0, sd_year)
      
      # Assign Treatment: Year Group Design
      if(design_type == "Year Group") {
        treat_val <- rbinom(1, 1, 0.5) 
      } else if (design_type == "Whole School") {
        treat_val <- s_tr
      } else {
        treat_val <- NA # Individual
      }
      
      # Generate Classes/Pupils
      n_classes <- get_n(classes_range)
      
      for(c in 1:n_classes) {
        n_pups <- get_n(class_size_range)
        
        # Treatment: Individual Design
        if(design_type == "Individual") {
          t_vec <- sample(rep(c(0,1), length.out = n_pups))
        } else {
          t_vec <- rep(treat_val, n_pups)
        }
        
        chunk <- data.frame(
          school_id = s_id,
          year_id = y_id,
          class_id = paste0(y_id, "_", c),
          treat = t_vec,
          # Sim Outcome includes random effects regardless of estimator
          score = (effect_size * t_vec) + s_u + u_year + rnorm(n_pups, 0, sd_resid)
        )
        data_list[[length(data_list) + 1]] <- chunk
      }
    }
  }
  
  data <- bind_rows(data_list)
  
  # --- 3. Attrition ---
  if (attrition_sch > 0) {
    drops <- sample(unique(data$school_id), size = round(n_schools * attrition_sch))
    data <- data %>% filter(!school_id %in% drops)
  }
  if (attrition_pup > 0) {
    data <- data %>% sample_frac(1 - attrition_pup)
  }
  
  # --- 4. Estimation ---
  results <- list()
  
  get_res <- function(model, name) {
    coefs <- summary(model)$coefficients
    idx <- grep("treat", rownames(coefs))
    if(length(idx) > 0) {
      est <- coefs[idx, 1]
      se <- coefs[idx, 2]
      t_val <- abs(est / se)
      sig <- as.integer(t_val > 1.96)
      return(data.frame(Estimator = name, Estimate = est, SE = se, Significant = sig))
    }
    return(NULL)
  }
  
  # OLS (Always available)
  if ("OLS" %in% estimators) {
    fit <- lm(score ~ treat, data = data)
    results[[length(results)+1]] <- get_res(fit, "OLS")
  }
  
  # MLM and CRE removed to reduce dependencies for Shinylive
  
  return(bind_rows(results))
}

# --- THE UI ---
ui <- page_sidebar(
  theme = bs_theme(bootswatch = "cosmo"),
  title = "Advanced Trial Simulator (OLS Only)",
  
  sidebar = sidebar(
    title = "Design Parameters",
    
    selectInput("design", "Randomization Level", 
                choices = c("Whole School", "Year Group", "Individual")),
    
    accordion(
      open = FALSE, 
      accordion_panel("Structure & Size",
        numericInput("n_schools", "Number of Schools", 30, min = 10, max = 200),
        sliderInput("years_range", "Year Groups per School (Range)", 
                    min = 1, max = 6, value = c(2, 4), step = 1),
        sliderInput("classes_range", "Classes per Year Group (Range)", 
                    min = 1, max = 5, value = c(1, 2), step = 1),
        sliderInput("size_range", "Pupils per Class (Range)", 
                    min = 10, max = 35, value = c(20, 30), step = 1)
      ),
      accordion_panel("Variance (ICC)",
        sliderInput("icc_sch", "School ICC", 0, 0.4, 0.15, step=0.01),
        sliderInput("icc_yr", "Year Group ICC", 0, 0.4, 0.05, step=0.01),
        helpText("Sum of ICCs must be < 1.0")
      ),
      accordion_panel("Effect & Attrition",
         sliderInput("effect", "True Effect (Cohen's d)", 0, 0.8, 0.2),
         sliderInput("att_sch", "School Attrition", 0, 0.3, 0),
         sliderInput("att_pup", "Pupil Attrition", 0, 0.3, 0)
      )
    ),
    
    # Checkbox disabled/restricted to OLS
    checkboxGroupInput("ests", "Estimators", 
                       choices = c("OLS"), 
                       selected = c("OLS")),
    helpText("MLM/CRE disabled for browser compatibility."),
    
    numericInput("reps", "Simulations", 30, min = 10, max = 200),
    actionButton("go", "Run Simulation", class = "btn-primary w-100")
  ),
  
  layout_column_wrap(
    width = 1,
    card(card_header("Statistical Power & Bias"), tableOutput("perfTable")),
    layout_column_wrap(
      width = 1/2,
      card(card_header("Bias Distribution"), plotOutput("biasPlot")),
      card(card_header("Standard Error Distribution"), plotOutput("sePlot"))
    )
  )
)

# --- THE SERVER ---
server <- function(input, output) {
  
  results <- eventReactive(input$go, {
    req(input$reps > 0)
    
    if((input$icc_sch + input$icc_yr) >= 1) {
      showNotification("Error: Total ICC (School + Year) must be < 1", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Simulating (OLS only)...', value = 0, {
      map_dfr(1:input$reps, function(i) {
        incProgress(1/input$reps)
        
        run_sim(
          n_schools = input$n_schools,
          years_range = input$years_range,
          classes_range = input$classes_range,
          class_size_range = input$size_range,
          attrition_sch = input$att_sch,
          attrition_pup = input$att_pup,
          effect_size = input$effect,
          icc_school = input$icc_sch,
          icc_year = input$icc_yr,
          estimators = input$ests,
          design_type = input$design
        )
      })
    })
  })
  
  output$perfTable <- renderTable({
    req(results())
    results() %>%
      group_by(Estimator) %>%
      summarise(
        `Mean Est` = mean(Estimate),
        Bias = mean(Estimate) - input$effect,
        `Mean SE` = mean(SE),
        `Power/Sig %` = mean(Significant) * 100
      ) %>%
      mutate(across(where(is.numeric), \(x) round(x, 3)))
  }, width = "100%")
  
  output$biasPlot <- renderPlot({
    req(results())
    ggplot(results(), aes(x = Estimate, fill = Estimator)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = input$effect, linetype = "dashed") +
      theme_minimal() +
      labs(title = "Effect Estimates", x = "Estimate (Cohen's d)")
  })
  
  output$sePlot <- renderPlot({
    req(results())
    ggplot(results(), aes(x = SE, fill = Estimator)) +
      geom_density(alpha = 0.5) +
      theme_minimal() +
      labs(title = "Standard Errors", x = "SE")
  })
}

shinyApp(ui, server)