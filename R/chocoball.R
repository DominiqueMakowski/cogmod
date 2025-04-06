# #' @title Force-Ball Simulation for Bounded Data
# #'
# #' @description
# #' Simulates data on the interval `[0, 1]` based on the "force-ball" concept. The ball's position
# #' on the x-axis (`x_ball`) and y-axis (`y_ball`) determines how the data is distributed:
# #' - A lower `y_ball` pushes data toward the edges (0 and 1).
# #' - A higher `y_ball` pulls data toward the center (0.5).
# #' - The horizontal position `x_ball` skews the distribution toward one side.
# #'
# #' @param n Number of data points to simulate.
# #' @param x_ball Horizontal position of the ball (range `[0, 1]`).
# #' @param y_ball Vertical position of the ball (positive, higher values pull data toward 0.5).
# #' @param alpha Sharpness of the force effect (default: 2). Higher values create sharper peaks.
# #'
# #' @return A vector of simulated data points in the range `[0, 1]`.
# #'
# #' @examples
# #' # Simulate data with the ball near the center
# #' x <- rforce_ball(10000, x_ball = 0.5, y_ball = 0.2)
# #' hist(x, breaks = 50, main = "Force-Ball Simulation", xlab = "x")
# #'
# #' # Simulate data with the ball pulling toward 0.5
# #' x <- rforce_ball(10000, x_ball = 0.5, y_ball = 5)
# #' hist(x, breaks = 50, main = "Force-Ball Simulation", xlab = "x")
# #'
# #' @export
# rforce_ball <- function(n, x_ball = 0.5, y_ball = 0.2, alpha = 2) {
#   # Validate inputs
#   if (x_ball < 0 || x_ball > 1)
#     stop("x_ball must be between 0 and 1")
#   if (y_ball <= 0)
#     stop("y_ball must be positive")
#   if (alpha <= 0)
#     stop("alpha must be positive")
#
#   # Define the force function
#   force_function <- function(x) {
#     1 / ((y_ball^2 + (x - x_ball)^2)^alpha)
#   }
#
#   # Generate a grid of x values and compute the force at each point
#   x_grid <- seq(0, 1, length.out = 1000)
#   density <- force_function(x_grid)
#
#   # Normalize the density to make it a valid probability distribution
#   density <- density / sum(density)
#
#   # Sample from the distribution using the computed density
#   sample(x_grid, size = n, replace = TRUE, prob = density)
# }
