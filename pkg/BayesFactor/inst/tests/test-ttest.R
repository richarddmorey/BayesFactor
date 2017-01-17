
context('t-test')

data('ToothGrowth')

test_that('ttest works', {
    ttestBF(ToothGrowth$len, ToothGrowth$dose)
    ttestBF(formula=len~supp, data=ToothGrowth)
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
})
