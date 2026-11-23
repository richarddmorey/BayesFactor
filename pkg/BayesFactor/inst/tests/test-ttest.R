
context('t-test')

data('ToothGrowth')

test_that('ttest works', {
    ttestBF(ToothGrowth$len, ToothGrowth$dose)
    ttestBF(formula=len~supp, data=ToothGrowth)
})

test_that('independent t tests use observed factor levels', {
    data <- data.frame(
        score = c(1, 2, 4, 5, 3, 5, 6, 7, 9, 10),
        group = factor(rep(c('a', 'c'), c(4, 6)), levels = c('c', 'b', 'a'))
    )
    expected <- ttestBF(formula = score ~ group, data = droplevels(data))
    actual <- ttestBF(formula = score ~ group, data = data)
    expect_equal(extractBF(actual, onlybf = TRUE), extractBF(expected, onlybf = TRUE))
})

test_that('rejects bad input', {
    expect_error(
        ttestBF(formula=len~dose, data=ToothGrowth),
        'Indep. groups t test requires a factor with exactly 2 levels.',
        fixed=TRUE
    )
    expect_error(
        ttestBF(formula=len~dose+supp, data=ToothGrowth),
        'Indep. groups t test can only support 1 factor as predictor.',
        fixed=TRUE
    )
    expect_error(
        ttestBF(formula=len~dose:supp, data=ToothGrowth),
        'Interaction terms are not allowed in t test.',
        fixed=TRUE
    )

    missing_group <- data.frame(score = 1:3, group = addNA(c('a', 'b', NA)))
    expect_error(
        ttestBF(formula = score ~ group, data = missing_group),
        'Predictors must not contain missing values.',
        fixed=TRUE
    )
})
