## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(irt)

## -----------------------------------------------------------------------------
item1 <- item(a = 1.2, b = -.8, c = .33, model = "3PL")
item1

## -----------------------------------------------------------------------------
item1 <- item(a = 1.2, b = -.8, c = .33, D = 1.7, model = "3PL")
item1

## -----------------------------------------------------------------------------
item1 <- item(a = 1.2, b = -.8, c = .33, D = 1.7, model = "3PL", 
              item_id = "ITM384", content = "Quadratic Equations")
item1

## -----------------------------------------------------------------------------
item1 <- item(a = 1.2, b = -.8, c = .33, D = 1.7, model = "3PL", 
              item_id = "ITM384", content = "Quadratic Equations", 
              misc = list(key = "A", 
                          enemies = c("ITM664", "ITM964"), 
                          seed_year = 2020, 
                          target_grade = "11")
              )
item1

## -----------------------------------------------------------------------------
plot(item1)

## -----------------------------------------------------------------------------
item2 <- item(b = -.8, model = "Rasch")
item2

## -----------------------------------------------------------------------------
item3 <- item(b = -.8, D = 1.7, model = "1PL")
item3

## -----------------------------------------------------------------------------
item4 <- item(a = 1.2, b = -.8, D = 1.702, model = "2PL")
item4

## -----------------------------------------------------------------------------
item5 <- item(a = 1.06, b = 1.76, c = .13, d = .98, model = "4PL", 
              item_id = "itm-5")
item5

## -----------------------------------------------------------------------------
item6 <- item(a = 1.22, b = c(-1.9, -0.37, 0.82, 1.68), model = "GRM", 
              item_id = "itm-6")
item6
plot(item6)

## -----------------------------------------------------------------------------
item7 <- item(a = 1.22, b = c(-1.9, -0.37, 0.82, 1.68), D = 1.7, model = "GPCM", 
              item_id = "itm-7")
item7

## -----------------------------------------------------------------------------
item8 <- item(b = c(-1.9, -0.37, 0.82, 1.68), model = "PCM")
item8

## -----------------------------------------------------------------------------
generate_item("3PL")
generate_item("2PL")
generate_item("Rasch")
generate_item("GRM")
# The number of categories of polytomous items can be specified:
generate_item("GPCM", n_categories = 5)

## -----------------------------------------------------------------------------
item1 <- item(a = 1.2, b = -.8, c = .33, D = 1.7, model = "3PL", 
              item_id = "ITM384", content = "Quadratic Equations")
item2 <- item(a = 0.75, b = 1.8, c = .21, D = 1.7, model = "3PL", 
              item_id = "ITM722", content = "Quadratic Equations")
item3 <- item(a = 1.06, b = 1.76, c = .13, d = .98, model = "4PL", 
              item_id = "itm-5")
t1 <- testlet(c(item1, item2, item3))
t1

## -----------------------------------------------------------------------------
t1 <- testlet(item1, item2, item3, testlet_id = "T1")
t1

## -----------------------------------------------------------------------------
item1 <- generate_item("3PL", item_id = "I1") 
item2 <- generate_item("3PL", item_id = "I2") 
item3 <- generate_item("3PL", item_id = "I3") 
ip1 <- itempool(item1, item2, item3)

## -----------------------------------------------------------------------------
item4 <- generate_item("GRM", item_id = "I4") 
item5 <- generate_item("3PL", item_id = "T1-I1") 
item6 <- generate_item("3PL", item_id = "T1-I2") 
t1 <- testlet(item5, item6, item_id = "T1")
ip2 <- itempool(item1, item2, item3, item4, t1)

## -----------------------------------------------------------------------------
n_item <- 6 # Number of items
ipdf <- data.frame(a = rlnorm(n_item), b = rnorm(n_item), 
                   c = runif(n_item, 0, .3))
ip3 <- itempool(ipdf)
ip3

# Scaling constant can be specified
ip4 <- itempool(ipdf, D = 1.7)
ip4

## -----------------------------------------------------------------------------
ipdf <- data.frame(
  item_id = c("Item_1", "Item_2", "Item_3", "Item_4", "Item_5", "Item_6"), 
  model = c("3PL", "3PL", "3PL", "GPCM", "GPCM", "GPCM"), 
  a = c(1.0253, 1.3609, 1.6617, 1.096, 0.9654, 1.3995), 
  b1 = c(NA, NA, NA, -1.112, -0.1709, -1.1324), 
  b2 = c(NA, NA, NA, -0.4972, 0.2778, -0.5242), 
  b3 = c(NA, NA, NA, -0.0077, 0.9684, NA), 
  D = c(1.7, 1.7, 1.7, 1.7, 1.7, 1.7), 
  b = c(0.7183, -0.4107, -1.5452, NA, NA, NA), 
  c = c(0.0871, 0.0751, 0.0589, NA, NA, NA), 
  content = c("Geometry", "Algebra", "Algebra", "Geometry", "Algebra", 
              "Algebra") 
)

ip5 <- itempool(ipdf)

## -----------------------------------------------------------------------------
as.data.frame(ip2)

## -----------------------------------------------------------------------------
item1 <- generate_item("3PL")
theta <- 0.84
# The probability of correct and incorrect response for `item1` at theta = 0.84
prob(item1, theta)

# Multiple theta values
prob(item1, theta = c(-1, 1))

# Polytomous items:
item2 <- generate_item(model = "GPCM")
prob(item2, theta = 1)
prob(item2, theta = c(-1, 0, 1))

## -----------------------------------------------------------------------------
ip <- generate_ip(model = "3PL", n = 7)
ip
prob(ip, theta = 0)
# When there are multiple theta values, a list where each element corresponds
# to a theta value returned. 
prob(ip, theta = c(-2, 0, 1))

## -----------------------------------------------------------------------------
# Plot ICC of each item in the item pool
plot(ip)

# Plot test characteristic curve
plot(ip, type = "tcc")

## -----------------------------------------------------------------------------
item1 <- generate_item("3PL")
info(item1, theta = -2)

# Multiple theta values
info(item1, theta = c(-1, 1))

# Polytomous items:
item2 <- generate_item(model = "GPCM")
info(item2, theta = 1)
info(item2, theta = c(-1, 0, 1))

## -----------------------------------------------------------------------------
ip <- generate_ip(model = "3PL", n = 7)
ip
info(ip, theta = 0)
info(ip, theta = c(-2, 0, 1))

## -----------------------------------------------------------------------------
# Plot information function of each item
plot_info(ip)
# Plot test information function
plot_info(ip, tif = TRUE)

## -----------------------------------------------------------------------------
# Generate an item pool 
ip <- generate_ip(model = "2PL", n = 10)
true_theta <- rnorm(5)
resp <- sim_resp(ip = ip, theta = true_theta, output = "matrix")

# Calculate raw scores
est_ability(resp = resp, ip = ip, method = "sum_score")
# Estimate ability using maximum likelihood estimation:
est_ability(resp = resp, ip = ip, method = "ml")
# Estimate ability using EAP estimation:
est_ability(resp = resp, ip = ip, method = "eap")
# Estimate ability using EAP estimation with a different prior 
# (prior mean = 0, prior standard deviation = 2):
est_ability(resp = resp, ip = ip, method = "eap", prior_pars = c(0, 2))

