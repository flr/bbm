# test-class.R - DESC
# FLPKG/tests/testthat/test-class.R

# Copyright %USER%, %YEAR%
# Author: %USER% <%EMAIL%>
#
# Distributed under the terms of the %LICENSE%

# CONTEXT
context("")

# TEST
test_that("", {

  # EXPECT
  expect_()

})


# ---

# CONTEXT class
context("CLASS")

# TEST prototype
out <- new("CLASS")

test_that(class(out), equals())
test_that(validObject(out), is_true())

# CONTEXT class constructor
context("CLASS constructor")

# TEST constructor
test_constructor <- function(msg, method, object) {

  out <- do.call(method, list(object=object))

  # check validObject
  test_that(validObject(out), is_true())

  # check class
  test_that(class(out), equals())
}

# TEST
test_constructor()

# class accessors
