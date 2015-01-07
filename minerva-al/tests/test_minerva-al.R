library("testthat")
source("../minerva-al_learn.R")

test_that(
  "Normalized echo"
  , {
    # echo >= -1 & echo <= 1
    # length(echo) == length(event)
  }
)

test_that(
  "Discrepency encoding"
  , {
    # Numeric values from p. 63, Jamieson, Crump, & Hannah (2012)
    a <- c(1, 1, 0, 0)
    x <- c(0, 0, 1, 1)
    e <- a + x
    c_prime <- c(0.4, 0.1, 0.6, 1.0)
    m <- learn(memory = "test", event = e, normalized_echo = c_prime, p_encode = 1)
    
    expect_that(m, is_a("matrix"))
    expect_that(length(m), equals(length(e)))
    expect_that(m, is_identical_to(matrix(c(0.6, 0.9, 0.4, 0.0), ncol = 4, dimnames = list("e", NULL))))
    
    # Numeric values from pp. 63-64, Jamieson, Crump, & Hannah (2012)
    e <- a
    m <- learn(memory = "test", event = e, normalized_echo = c_prime, p_encode = 1)
    
    expect_that(m, is_a("matrix"))
    expect_that(length(m), equals(length(e)))
    expect_that(m, is_identical_to(matrix(c(0.6, 0.9, -0.6, -1.0), ncol = 4, dimnames = list("e", NULL))))
  }
)


# References
#
# Jamieson, R. K., Crump, M. J. C., & Hannah, S. D. (2012). An instance theory of associative learning. Learning & Behavior, 40(1), 61–82. doi:10.3758/s13420-011-0046-2
