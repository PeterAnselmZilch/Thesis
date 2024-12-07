ENV["GUROBI_HOME"] = "C:\\gurobi1101\\win64"

# import packages  (only run this once to add packages)
import Pkg
Pkg.add("YAML")
Pkg.add("JuMP")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Tables")
Pkg.add("Statistics")
Pkg.add("StatsPlots")
Pkg.add("Dates")
Pkg.add("Gurobi")
Pkg.add("IterTools")

# activate the packages  (run this every time you start a new session)
using YAML
using JuMP
using Gurobi
using CSV
using DataFrames
using Tables
using Statistics
using StatsPlots
using Dates
using IterTools
using Plots



