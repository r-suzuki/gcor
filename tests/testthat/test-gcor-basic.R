test_that("gcor example scalar regression", {
  x <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  y <- c(1, 2, 3, 4, 5, 3, 4, 5, 6, 7)

  got <- gcor(x, y)
  expected <- 0.5345224838248488

  expect_type(got, "double")
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("gcor example matrix regression", {
  df <- data.frame(
    x = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
    y = c(1, 2, 3, 4, 5, 3, 4, 5, 6, 7),
    z = c("a", "a", "b", "b", "c", "c", "d", "d", "e", "e")
  )

  got <- gcor(df)

  expect_true(is.matrix(got))
  expect_identical(colnames(got), c("x", "y", "z"))
  expect_identical(rownames(got), c("x", "y", "z"))

  expected <- matrix(
    c(
      1.0,      0.534522, 0.838289,
      0.534522, 1.0,      0.763233,
      0.838289, 0.763233, 1.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("x", "y", "z"), c("x", "y", "z"))
  )

  expect_equal(got, expected, tolerance = 5e-5)
})

test_that("gcor example scalar matches matrix entry", {
  x <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
  y <- c(1, 2, 3, 4, 5, 3, 4, 5, 6, 7)
  df <- data.frame(
    x = x,
    y = y,
    z = c("a", "a", "b", "b", "c", "c", "d", "d", "e", "e")
  )

  scalar <- gcor(x, y)
  mat <- gcor(df)

  expect_equal(as.numeric(scalar), as.numeric(mat["x", "y"]), tolerance = 1e-12)
})

test_that("gcor handles constant columns as in examples", {
  df <- data.frame(
    x = rep(0, 10),
    y = rep(1, 10),
    z = c("a", "a", "b", "b", "c", "c", "d", "d", "e", "e")
  )

  got <- gcor(df)

  expected <- matrix(
    c(
      1.0, 1.0, 0.0,
      1.0, 1.0, 0.0,
      0.0, 0.0, 1.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("x", "y", "z"), c("x", "y", "z"))
  )

  expect_equal(got, expected, tolerance = 5e-5)
})
