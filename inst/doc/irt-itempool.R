## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(irt)

## -----------------------------------------------------------------------------
# Create an item pool of 2PL items:
ip_dtf <- data.frame(
  a = c(1.1821, 0.6645, 0.8994, 1.0731, 1.0252, 1.2325, 0.9278, 1.0967), 
  b = c(0.4185, -0.5992, 0.2193, 0.8823, 0.4652, 1.4006, -1.1193, -0.3747))

ip <- itempool(ip_dtf, model = "2PL")
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame( 
  a = c(1.1821, 0.6645, 0.8994, 1.0731, 1.0252, 1.2325, 0.9278, 1.0967), 
  b = c(0.4185, -0.5992, 0.2193, 0.8823, 0.4652, 1.4006, -1.1193, -0.3747),
  item_id = c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"), 
  content = c("Geometry", "Geometry", "Algebra", "Algebra", "Algebra", 
              "Geometry", "Algebra", "Algebra")
  )

ip <- itempool(ip_dtf, model = "2PL")
ip

## -----------------------------------------------------------------------------
ip <- itempool(ip_dtf, model = "2PL", D = 1.7)
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame(a = c(0.9303, 1.9423, 0.8417, 1.2622), 
                     b = c(1.3515, -0.5039, -1.7263, 1.3125), 
                     c = c(0.2301, 0.2224, 0.0967, 0.0112))
ip <- itempool(ip_dtf, model = "3PL", D = 1.7)
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame(a = c(1.8619, 1.2458, 1.3213, 0.6174, 1.3625), 
                     b1 = c(-0.3666, -0.9717, -1.1588, 0.1093, 0.0858), 
                     b2 = c(0.3178, 0.2458, -0.4978, 0.6437, 0.5161), 
                     b3 = c(1.0384, 1.2382, 1.2787, 1.3609, 1.2145))

ip <- itempool(ip_dtf, model = "GPCM")
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame(a = c(1.175, 0.981, 1.0625, 0.9545, 0.7763), 
                     b1 = c(-0.9633, -0.4098, -0.298, 0.0576, -0.5342), 
                     b2 = c(-0.6213, NA, 0.4792, 0.538, 0.0363), 
                     b3 = c(0.5938, NA, NA, 0.9815, NA), 
                     b4 = c(NA, NA, NA, 1.3351, NA))
ip <- itempool(ip_dtf, model = "GRM", D = 1.702)
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame(a = c(1.1152, 0.8231, 0.9527, 0.6423), 
                     b = c(0.234, 0.0219,  0.7424, -0.3426), 
                     d1 = c(0.0081, 0.8569, -1.5181, -0.8458), 
                     d2 = c(0.3392, NA, -0.1978, 0.3756), 
                     d3 = c(NA, NA, 0.1677, NA))
ip <- itempool(ip_dtf, model = "GPCM2", D = 1.702)
ip

## -----------------------------------------------------------------------------
ip_dtf <- data.frame(
  model = c("3PL", "3PL", "3PL", "GPCM", "GPCM"), 
  a = c(1.6242, 0.9471, 1.4643, 0.6582, 1.0234), 
  b = c(0.4563, -0.2994, -0.3027, NA, NA), 
  c = c(0.0156, 0.0339, 0.1243, NA, NA), 
  b1 = c(NA, NA, NA, -1.1532, -1.2171), 
  b2 = c(NA, NA, NA, -0.5384, -0.3992), 
  b3 = c(NA, NA, NA, 0.0591, 0.1431), 
  b4 = c(NA, NA, NA, NA, 1.52), 
  D = c(1.7, 1.7, 1.7, 1, 1))
ip <- itempool(ip_dtf)
ip

## -----------------------------------------------------------------------------
ip <- c(item(model = "3PL", a = 2.09, b = 1.17, c = 0.25, item_id = "i1"), 
        item(model = "3PL", a = 0.59, b = 0.77, c = 0.13, item_id = "i2"), 
        item(model = "3PL", a = 1.67, b = 1.05, c = 0.04, item_id = "i3"), 
        item(model = "3PL", a = 0.84, b = -1.8, c = 0.24, item_id = "i4"), 
        item(model = "GPCM", a = 1.96, b = c(-0.94, -0.09, 0.25), item_id = "i5"), 
        item(model = "GPCM", a = 0.59, b = c(0.07, 1.46), item_id = "i6"), 
        item(model = "GPCM", a = 0.73, b = c(-1.2, -0.78, 0.2, 1.8), item_id = "i7"))
ip

## -----------------------------------------------------------------------------
# Create a testlet object with three items. 
t1 <- testlet(c(item(model = "3PL", a = 2.09, b = 1.17, c = 0.25, item_id = "i1"), 
                item(model = "3PL", a = 0.59, b = 0.77, c = 0.13, item_id = "i2"), 
                item(model = "3PL", a = 1.67, b = 1.05, c = 0.04, item_id = "i3")), 
              testlet_id = "Testlet-932")
# Create another testlet object with two items. 
t2 <- testlet(c(item(model = "3PL", a = 0.84, b = -1.8, c = 0.24, item_id = "i4"), 
                item(model = "GPCM", a = 1.96, b = c(-0.94, -0.09, 0.25), 
                     item_id = "i5")), 
              testlet_id = "Testlet-77")
# Standalone items to be added:
i6 <- item(model = "GPCM", a = 0.59, b = c(0.07, 1.46), item_id = "i6")
i7 <- item(model = "GPCM", a = 0.73, b = c(-1.2, -0.78, 0.2, 1.8), item_id = "i7")

# Combine all items and testlets:
ip_testlet <- c(t1, t2, i6, i7)
ip_testlet

## -----------------------------------------------------------------------------
ip1 <- itempool(data.frame( 
  a = c(1.1821, 0.6645, 0.8994, 1.0731, 1.0252, 1.2325, 0.9278, 1.0967), 
  b = c(0.4185, -0.5992, 0.2193, 0.8823, 0.4652, 1.4006, -1.1193, -0.3747),
  item_id = c("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"), 
  content = c("Geometry", "Geometry", "Algebra", "Algebra", "Algebra", 
              "Geometry", "Algebra", "Algebra")
  ))
  
ip_mixed <- itempool(data.frame(
  model = c("3PL", "3PL", "3PL", "GPCM", "GPCM"), 
  a = c(1.6242, 0.9471, 1.4643, 0.6582, 1.0234), 
  b = c(0.4563, -0.2994, -0.3027, NA, NA), 
  c = c(0.0156, 0.0339, 0.1243, NA, NA), 
  b1 = c(NA, NA, NA, -1.1532, -1.2171), 
  b2 = c(NA, NA, NA, -0.5384, -0.3992), 
  b3 = c(NA, NA, NA, 0.0591, 0.1431), 
  b4 = c(NA, NA, NA, NA, 1.52), 
  D = c(1.7, 1.7, 1.7, 1.7, 1.7)))

ip_testlet <- c(
  testlet(c(item(model = "3PL", a = 2.09, b = 1.17, c = 0.25, item_id = "i1"), 
            item(model = "3PL", a = 0.59, b = 0.77, c = 0.13, item_id = "i2"), 
            item(model = "3PL", a = 1.67, b = 1.05, c = 0.04, item_id = "i3")), 
          testlet_id = "Testlet-932"), 
  item(model = "GPCM", a = 0.59, b = c(0.07, 1.46), item_id = "i6"),
  testlet(c(item(model = "3PL", a = 0.84, b = -1.8, c = 0.24, item_id = "i4"), 
            item(model = "GPCM", a = 1.96, b = c(-0.94, -0.09, 0.25), 
                 item_id = "i5")), 
          testlet_id = "Testlet-77"), 
  item(model = "GPCM", a = 0.73, b = c(-1.2, -0.78, 0.2, 1.8), item_id = "i7"))

## -----------------------------------------------------------------------------
ip <- c(ip1, ip_mixed)
ip

## -----------------------------------------------------------------------------
# Subset only the first element of the item pool
ip1[1]

# Create an Itempool using the first and third element:
ip1[c(1, 3)] # Order is important
ip1[c(3, 1)]

# Create an Itempool using all but the second element: 
ip1[-2]

# Subsetting using item ID's:
ip1[c("i2", "i1")]

# Subsetting using logical operators:
ip1[ip1$b < 0]

# Select items with information values larger than 0.2 at theta = 1:
ip1[info(ip1, theta = 1) > 0.2]


## -----------------------------------------------------------------------------
# Extract the second element
ip1[[2]]

# Extract a testlet
ip_testlet[[3]]

## -----------------------------------------------------------------------------
ip_new <- ip1
# Replace the second item with a new item
ip_new[[2]] <- item(a = 1, b = c(-1, 0, 1), model = "GRM", item_id = "NewItm4",
                    D = 1.7, content = "Quadratic Functions")
ip_new

## ----eval=FALSE---------------------------------------------------------------
#  ?`$,Itempool-method`

## -----------------------------------------------------------------------------
# Extract the ID's of the items within an item pool
ip1$item_id

# Extract the contents of the items within an item pool
ip1$content

# Extract the models of the items within an item pool
ip1$model
ip_mixed$model

# Maximum possible score of items
ip1$max_score
ip_mixed$max_score
ip_testlet$max_score

# Maximum scores of each standalone item
ip1$item_max_score
ip_mixed$item_max_score
ip_testlet$item_max_score



## -----------------------------------------------------------------------------
ip1$a
ip1$b
ip1$c
ip1$D

ip_mixed$a
ip_mixed$b
ip_mixed$b1
ip_mixed$b4
ip_mixed$D

## -----------------------------------------------------------------------------
# Extract the number of items within an item pool
ip1$n
# In ip_testlet, there are two testlets and two standalone items. Within those
# two testlets, there are a total of 5 items. At total there are 7 items.
ip_testlet$n

## -----------------------------------------------------------------------------
ip_new <- ip1
ip_new$item_id <- paste0("Question-", 1:length(ip_new))
ip_new$content <- c("M", "M", "R", "M", "E", "R", "E", "M")


## -----------------------------------------------------------------------------
ip_new$a <- 1
ip_new
ip_new$b <- rnorm(length(ip_new))
ip_new

## ----eval=FALSE, include=FALSE------------------------------------------------
#  
#  ip <- generate_ip(model = c("3PL", "3PL", "3PL", "3PL", "GPCM", "GPCM", "GPCM"))
#  dput(ip)
#  dput(as.data.frame(generate_ip(model = c("3PL", "3PL", "3PL", "3PL", "3PL"))))
#  dput(as.data.frame(generate_ip(model = "GPCM", n = 5)))
#  dput(as.data.frame(generate_ip(model = "GPCM2", n = 4, n_categories = c(3, 2, 4, 3))))
#  dput(as.data.frame(generate_ip(model = "GPCM", n = 5, n_categories = c(5, 2, 3, 4, 3))))
#  
#  ip <- generate_ip(model = c("3PL", "3PL", "3PL", "GPCM", "GPCM"))
#  dput(data.frame(ip))
#  # list(key = c("A", "C", "D", "B", "B", "A", "D", "C"),
#  #                            grade = c(11, 9, 9, 10, 10, 10, 11, 9))

