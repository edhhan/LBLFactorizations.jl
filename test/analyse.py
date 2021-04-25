import os
import numpy as np
from statistics import mean, stdev
import matplotlib.pyplot as plt
import numpy as np
import scipy 
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import csv


# Read CSV data file
# Place yourself in the working directory ../julia-factorization-bunch-kaufamn
working_dir = os.getcwd()
analysis_path = os.path.join(working_dir, "test")
os.chdir(analysis_path)


# Initialize data arrays
bk_fact_arr = []
bp_fact_arr = []
rk_fact_arr = []

bk_solver_arr = []
bp_solver_arr = []
rk_solver_arr = []

dim_arr = []

with open('data_analysis.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:

            # Factorization indices            
            bk_fact_idx = row.index("bkaufmann")
            bp_fact_idx = row.index("bparlett")
            rk_fact_idx = row.index("rook")

            # Solver indices
            bk_solver_idx = row.index("bkaufmann")
            bp_solver_idx = row.index("bparlett")
            rk_solver_idx = row.index("rook")

            # Dimension
            dim_idx = row.index("dim")

            line_count += 1
        else:
            bk_fact_arr.append(row[bk_fact_idx])
            bp_fact_arr.append(row[bp_fact_idx])
            rk_fact_arr.append(row[rk_fact_idx])

            bk_solver_arr.append(row[bk_fact_idx])
            bp_solver_arr.append(row[bp_solver_idx])
            rk_solver_arr.append(row[rk_solver_idx])
            
            dim_arr.append(row[dim_idx])

            line_count += 1

    print("Number of lines processed {}".format(line_count))


bk_fact_arr = np.array(bk_fact_arr, dtype=np.float64).reshape(-1, 1)
bp_fact_arr = np.array(bp_fact_arr, dtype=np.float64).reshape(-1, 1)
rk_fact_arr = np.array(rk_fact_arr, dtype=np.float64).reshape(-1, 1)

bk_solver_arr = np.array(bk_solver_arr, dtype=np.float64).reshape(-1, 1)
bp_solver_arr = np.array(bp_solver_arr, dtype=np.float64).reshape(-1, 1)
rk_solver_arr = np.array(rk_solver_arr, dtype=np.float64).reshape(-1, 1)

dim_arr = np.array(dim_arr, dtype=np.int64).reshape(-1, 1)


# Data 
fig, ax = plt.subplots()
plt.scatter(dim_arr, bk_fact_arr, label="bk")
plt.scatter(dim_arr, bp_fact_arr, label="bp")
plt.scatter(dim_arr, rk_fact_arr, label="rk")
ax.set_xlabel("Dimension [n]")
ax.set_ylabel("Temps [s]")
#plt.title("Temps de calcul en fonction de la dimension")
plt.legend()
plt.show()


#### Power-test ####
x_log = np.log2(dim_arr)

# 1) Power-test : bk
bk_log = np.log2(bk_fact_arr)

# Linear regression
bk_model_log = LinearRegression()
bk_model_log.fit(x_log, bk_log)

# Parameters
bk_slope_log = float(bk_model_log.coef_)
bk_intercept_log = float(bk_model_log.intercept_)

# Print
fig, ax = plt.subplots()
plt.scatter(x_log, bk_log, label="Data")
plt.plot(x_log, bk_slope_log*x_log + bk_intercept_log, label="Regression $log_2(y)=mlog_2(x)+b$", color='r')
ax.set_xlabel("Dimension $[log_2(n)]$", fontsize=14)
ax.set_ylabel("Temps $[log_2(s)]$", fontsize=14)
#plt.title("Log-Log regression for the factorization time with the bunch-kaufmann strategy")
plt.legend()
plt.show()

# Compute complexity quantities
bk_exposant = round(float(bk_slope_log),3)
bk_constante = round(2**float(bk_intercept_log),3)
print("Bkaufmann polynomial complexity degree : {}".format(bk_exposant))




# 2) Power-test : bp
bp_log = np.log2(bp_fact_arr)

# Linear regression
bp_model_log = LinearRegression()
bp_model_log.fit(x_log, bp_log)

# Parameters
bp_slope_log = float(bp_model_log.coef_)
bp_intercept_log = float(bp_model_log.intercept_)

# Print
fig, ax = plt.subplots()
plt.scatter(x_log, bp_log, label="Data")
plt.plot(x_log, bp_slope_log*x_log + bp_intercept_log, label="Regression $log_2(y)=mlog_2(x)+b$", color='r')
ax.set_xlabel("Dimension $[log_2(n)]$", fontsize=14)
ax.set_ylabel("Temps $[log_2(s)]$", fontsize=14)
#plt.title("Log-Log regression for the factorization time with the bunch-parlett strategy")
plt.legend()
plt.show()

# Compute complexity quantities
bp_exposant = round(float(bp_slope_log),3)
bp_constante = round(2**float(bp_intercept_log),3)
print("Bparlett polynomial complexity degree : {}".format(bp_exposant))




# 3) Power-test : rk
rk_log = np.log2(rk_fact_arr)

# Linear regression
rk_model_log = LinearRegression()
rk_model_log.fit(x_log, rk_log)

# Parameters
rk_slope_log = float(rk_model_log.coef_)
rk_intercept_log = float(rk_model_log.intercept_)

# Print
fig, ax = plt.subplots()
plt.scatter(x_log, rk_log, label="Data")
plt.plot(x_log, rk_slope_log*x_log + rk_intercept_log, label="Regression $log_2(y)=mlog_2(x)+b$", color='r')
ax.set_xlabel("Dimension $[log_2(n)]$", fontsize=14)
ax.set_ylabel("Temps $[log_2(s)]$", fontsize=14)
#plt.title("Log-Log regression for the factorization time with the rook strategy")
plt.legend()
plt.show()

# Compute complexity quantities
rk_exposant = round(float(rk_slope_log),3)
rk_constante = round(2**float(rk_intercept_log),3)
print("Rook polynomial complexity degree : {}".format(rk_exposant))