
######################################################################################
################## Authors: Bozzolan,Verkerk,Zudin  ##################################
######################################################################################
######################## Task 4.6 HS #################################################
######################################################################################
######################################################################################

import pandas as pd
from pathlib import Path
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import r2_score

# Paths files
input_path = Path(r"C:\Users\nibozzol\OneDrive - European Forest Institute\Projects\HoliSoils\WP4\me4soc\simulations\output_Sergey\Ribera\holisoils_ribera_reduced.csv")
output_path = Path(r"C:\Users\nibozzol\OneDrive - European Forest Institute\Projects\HoliSoils\WP4\me4soc\simulations\output_Sergey\Ribera\model_results_holisoils_ribera_reduced_i.csv")

# Load data
df = pd.read_csv(input_path)

# Features and target
X = df[["harvest_ratio", "mortality"]] # independent variables
y = df["year20"]                        # dependent (target value)

# ---- Linear regression ----
lin_model = LinearRegression().fit(X, y)  # Create a Linear Regression model and fit it with X and y . This finds the best coefficients (slopes) and intercept.
y_lin = lin_model.predict(X)              # Use the fitted model to predict y values (SOC) for the same inputs X.
r2_lin = r2_score(y, y_lin)

linear_formula = (
    f"SOC = {lin_model.intercept_:.4f}"
    f" + {lin_model.coef_[0]:.4f}*HarvestRatio"  # eg. SOC = 84.26 + -0.80*HarvestRatio + 0.12*Mortality. (4f 4 decimal places)
    f" + {lin_model.coef_[1]:.4f}*Mortality"
)


# ---- Linear regression + interaction (no squares) ----
X_int = X.copy()
X_int["harvest_mortality"] = X_int["harvest_ratio"] * X_int["mortality"]

lin_int_model = LinearRegression().fit(X_int, y)
y_lin_int = lin_int_model.predict(X_int)
r2_lin_int = r2_score(y, y_lin_int)

linear_int_formula = (
    f"SOC = {lin_int_model.intercept_:.4f}"
    f" + {lin_int_model.coef_[0]:.6f}*HarvestRatio"
    f" + {lin_int_model.coef_[1]:.6f}*Mortality"
    f" + {lin_int_model.coef_[2]:.6f}*(HarvestRatio*Mortality)"
)

# ---- Polynomial regression (quadratic + interactions) ---- harvest_ratio # mortality #  # harvest_ratio² #  # harvest_ratio * mortality # mortality²
poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform(X)                           # Apply transformation: convert the 2-column input into a bigger matrix with polynomial features.
poly_model = LinearRegression().fit(X_poly, y)           # Fit a new linear regression, but now using polynomial features as predictors
y_poly = poly_model.predict(X_poly)                     # Predict values with the polynomial regression.
r2_poly = r2_score(y, y_poly)

poly_coeffs = ", ".join([f"{coef:.4f}" for coef in poly_model.coef_])
poly_formula = f"Quadratic coefficients: {poly_model.intercept_:.4f}, {poly_coeffs}"


# ---- Save results ----
results = pd.DataFrame({
    "Model": ["Linear","Linear + Interaction", "Polynomial"],
    "Formula": [linear_formula,linear_int_formula, poly_formula],
    "R2 Score": [r2_lin, r2_lin_int, r2_poly]
})

results.to_csv(output_path, index=False)

