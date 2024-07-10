#' Define a tidymodels recipe
#' @return A recipe
#' @description # We normalised numeric data to have a standard deviation of one and a mean of zero. We created dummy variables from categorical variables (i.e. we converted nominal data into numeric binary model terms for the levels of the original data.). We omitted observations with missing values.
define_recipe <- function(x) {
    recipe(outcome ~ ., data = x) %>%
        update_role(id, new_role = "ID") %>%
        step_zv(all_numeric(), -all_outcomes()) %>%
        step_normalize(all_numeric(), -all_outcomes()) %>%
        step_dummy(all_nominal(), -all_outcomes()) %>% # Create dummy variables, see https://stats.stackexchange.com/questions/136085/can-glmnet-logistic-regression-directly-handle-factor-categorical-variables-wi and https://recipes.tidymodels.org/reference/step_dummy.html
        step_naomit(everything())
}


#' Tune grid with parallel processing
parallel_tune <- function(workflow, train_data, lambda_grid) {
	doParallel::registerDoParallel()
	tune_grid(workflow, resamples = bootstraps(train_data), grid = lambda_grid)
}

