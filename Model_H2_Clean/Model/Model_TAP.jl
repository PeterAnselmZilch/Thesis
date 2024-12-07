# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# this model assumes a take-as-produced profile for the hydrogen demand during the whole simulation, through a yearly volume agreement
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Define the model and the solver 
m = Model(Gurobi.Optimizer) 
set_optimizer_attribute(m, "Method", 2)
set_optimizer_attribute(m, "NodeMethod", 2)
set_optimizer_attribute(m, "StartNodeLimit", -2)

# Define the components of the path for the data and the timeseries
dir = "" # the directory where this folder is stored
filename = "Data\\inputs.yaml" # the file with the data 
full_path = joinpath(dir, filename) # the full path to the file

# load the data and timeseries
data = YAML.load_file(full_path) # load the data from the yaml file
ts = CSV.read("Data\\merged_df_all.csv", DataFrame) # load the timeseries data

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# defining economic parameters
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Define the discount rate
r = 0.05  # example discount rate

# Define the electrolyser parameters
capacity = 100 # MW e
investment_cost = 3050000 # €/MW e 
operational_cost =  75320 # €/MW e / year
TSO_cost_power = 143570 # €/MWe / year 
HyNetwork_cost =  21130 # €/MWe / year

# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Model
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# defines the model using the parameters and data for timeseries


function complete_model(n_bundles, volume, capacity_solar, capacity_wind, scenario)

    # defining sets for timeseries, note: each year consists here of 8760 hours for simplicity 
    function define_sets!(m::Model, data::Dict)
        # create dictionary to store sets
        m.ext[:sets] = Dict()

        # add sets to dictionary
        m.ext[:sets][:JY] = 1:10 # years
        m.ext[:sets][:JH] = 1:8760 # hours in each year

        return m
    end

    # process the time series data so the model can use it
    function process_time_series_data(m::Model, data::Dict, ts::DataFrame)
        # extract the relevant sets
        JH = m.ext[:sets][:JH] # Time steps / Hours
        JY = m.ext[:sets][:JY] # Years 

        # create dictionary to store time series 
        m.ext[:timeseries] = Dict() 

        # add time series to dictionary 
        m.ext[:timeseries][:DAP] = [ts[:,"price_day_ahead_$scenario"][jh+8760*(jy-1)] for jh in JH, jy in JY] # day ahead prices
        m.ext[:timeseries][:cf_wind] = [ts[:,"wind_cap_$scenario"][jh+8760*(jy-1)] for jh in JH, jy in JY] # capacity factors wind   
        m.ext[:timeseries][:cf_solar] = [ts[:,"solar_cap_$scenario"][jh+8760*(jy-1)] for jh in JH, jy in JY] # capacity factors solar
        return m
    end

    # extract the parameters from the inputs.yaml file
    function process_parameters!(m::Model, data::Dict)
        # generate a dictonary "parameters"
        m.ext[:parameters] = Dict()

        # Electrolyser parameters
        electrolyser = data["electrolyser"]
        m.ext[:parameters][:capacity_electrolyser] = electrolyser["capacity"] # MW e
        m.ext[:parameters][:power_consumption] = electrolyser["power_consumption"] # kWh e / kg H2
        m.ext[:parameters][:power_consumption_compressor] = electrolyser["power_consumption_compressor"] # kWh e / kg H2
        m.ext[:parameters][:minimum_partial_load] = electrolyser["minimum_partial_load"] # 

        # power supply parameters
        power_supply = data["power_supply"]
        m.ext[:parameters][:Q_W] = capacity_wind # MW e
        m.ext[:parameters][:Q_PV] = capacity_solar # MW e 
        m.ext[:parameters][:λ_W] = power_supply["price_wind"] # €/MWh
        m.ext[:parameters][:λ_PV] = power_supply["price_solar"] # €/MWh

        # H2 Storage parameters
        h2_storage = data["H2_storage"]
        m.ext[:parameters][:storage_capacity] = h2_storage["storage_capacity"] # kg H2 / bundle
        m.ext[:parameters][:storage_capacity_flow_rate] = h2_storage["storage_capacity_flow_rate"] # kg H2 / h / bundle
        m.ext[:parameters][:n_bundles] = n_bundles 
        # m.ext[:parameters][:n_bundles] = h2_storage["n_bundles"] # storage capacity flow rate MWh

        # H2 Demand parameters
        m.ext[:parameters][:H_dem_bl] = (volume/100) * 1960 # H2 demand baseload [kg H2 / h] - as a percentage of the maximum production. 
        return m.ext
    end

    define_sets!(m, data)
    process_time_series_data(m, data, ts)
    process_parameters!(m, data)


    function build_model!(m::Model)
        # Clear m.ext entries "variables", "expressions" and "constraints"
        m.ext[:variables] = Dict()
        m.ext[:expressions] = Dict()
        m.ext[:constraints] = Dict()
    
        # Extract sets
        JH = m.ext[:sets][:JH] # Time steps / Hours
        JY = m.ext[:sets][:JY] # Years

        # Extract time series
        DAP = m.ext[:timeseries][:DAP] # Day-Ahead Price
        cf_wind = m.ext[:timeseries][:cf_wind] # capacity factors wind
        cf_solar = m.ext[:timeseries][:cf_solar] # capacity factors solar

        # Extract parameters with symbols
        Q_E = m.ext[:parameters][:capacity_electrolyser] # MW e
        power_consumption = m.ext[:parameters][:power_consumption] # MWh e / kg H2
        power_consumption_compressor = m.ext[:parameters][:power_consumption_compressor] # MWh e / kg H2
        Q_min = m.ext[:parameters][:minimum_partial_load] # [MW e]
        Q_W = m.ext[:parameters][:Q_W] # wind capacity [MW]
        Q_PV = m.ext[:parameters][:Q_PV] # solar capacity [MW]
        Q_S = m.ext[:parameters][:storage_capacity] # [kg H2]
        Q_I = m.ext[:parameters][:storage_capacity_flow_rate] # [kg H2 / h]
        n_bundles = m.ext[:parameters][:n_bundles] #  []
        H_dem_bl = m.ext[:parameters][:H_dem_bl] #  [kg H2 / h]
        λ_W = m.ext[:parameters][:λ_W] # [eu / MWh]
        λ_PV = m.ext[:parameters][:λ_PV] # [eu / MWh]w

        # create variables 
        P_W = m.ext[:variables][:P_W] = @variable(m, [jh in JH,jy in JY], lower_bound=0, base_name="power_wind") # 
        P_PV = m.ext[:variables][:P_PV] = @variable(m, [jh in JH,jy in JY], lower_bound=0, base_name="power_solar") # 
        P_Grid = m.ext[:variables][:P_Grid] = @variable(m, [jh in JH,jy in JY], base_name="power_grid") # 
        P_el = m.ext[:variables][:P_el] = @variable(m, [jh in JH,jy in JY], base_name="power_electrolyser") #
        # P_Curt = m.ext[:variables][:P_Curt] = @variable(m, [jh in JH,jy in JY], lower_bound=0, base_name="power_curtailment") # 
        C_P = m.ext[:variables][:C_P] = @variable(m, [jh in JH,jy in JY], base_name="cost_power") # 
        H_Pr = m.ext[:variables][:H_Pr] = @variable(m, [jh in JH,jy in JY], lower_bound=0, base_name="hydrogen_production") # 
        H_S = m.ext[:variables][:H_S] = @variable(m, [jh in JH,jy in JY], base_name="hydrogen_storage") # 
        H_S_in = m.ext[:variables][:H_S_in] = @variable(m, [jh in JH,jy in JY], base_name="hydrogen_storage_inflow") # 
        H_del = m.ext[:variables][:H_del] = @variable(m, [jh in JH,jy in JY], base_name="hydrogen_delivery") # 
        z_on = m.ext[:variables][:z_on] = @variable(m, [jh in JH,jy in JY], Bin, base_name="on") # 

        # objective function
        m.ext[:objective] = @objective(m, Min, 
            sum(C_P[jh,jy] for jh in JH, jy in JY)
        ) 

        # constraints 

        # power costs
        m.ext[:constraints][:power_costs] = @constraint(m, [jh=JH,jy=JY], C_P[jh,jy] == DAP[jh,jy] * P_Grid[jh,jy] + λ_W * P_W[jh,jy] + λ_PV * P_PV[jh,jy]) # 

        # power generation from wind and solar 
        m.ext[:constraints][:wind_power] = @constraint(m, [jh=JH,jy=JY], P_W[jh,jy] == Q_W * cf_wind[jh,jy]) # 
        m.ext[:constraints][:solar_power] = @constraint(m, [jh=JH,jy=JY], P_PV[jh,jy] == Q_PV * cf_solar[jh,jy]) # 

        # power balance = equality constraint
        m.ext[:constraints][:power_balance] = @constraint(m, [jh=JH,jy=JY], P_el[jh,jy] == P_W[jh,jy] + P_PV[jh,jy] + P_Grid[jh,jy] ) # 

        # Hydrogen production
        m.ext[:constraints][:hydrogen_production] = @constraint(m, [jh=JH,jy=JY], H_Pr[jh,jy] == P_el[jh,jy] / (power_consumption + power_consumption_compressor) *  z_on[jh,jy]) # hydrogen production = H_pr [kg] = P_el [MW] * power_consumption [MWh/kg] * z_on [binary]        
        
        # electrolyser capacity constraint (lower and upper bound)
        m.ext[:constraints][:electrolyser_capacity_min] = @constraint(m, [jh=JH,jy=JY], P_el[jh,jy] >= Q_min * Q_E * z_on[jh,jy]) # electrolyser minimum partial load   *should have the same unit as capacity
        m.ext[:constraints][:electrolyser_capacity_max] = @constraint(m, [jh=JH,jy=JY], P_el[jh,jy] <= Q_E) # electrolyser capacity

        # hydrogen demand for tap profile
        m.ext[:constraints][:hydrogen_demand_tap] = @constraint(m, [jy=JY], sum(H_del[jh,jy] for jh in JH) == H_dem_bl * 8760) # total hydrogen delivered in a year = hourly volume times number of timesteps
    
        # hydrogen balance constraint
        m.ext[:constraints][:hydrogen_balance] = @constraint(m, [jh=JH,jy=JY], H_del[jh,jy] == H_Pr[jh,jy] - H_S_in[jh,jy]) # hydrogen balance 

        # hydrogen storage capacity constraints
        m.ext[:constraints][:hydrogen_storage_min] = @constraint(m, [jh=JH,jy=JY], H_S[jh,jy] >= 0) # storage level is always above 0 
        m.ext[:constraints][:hydrogen_storage_max] = @constraint(m, [jh=JH,jy=JY], H_S[jh,jy] <= Q_S * n_bundles) # storage level is always below Q_S
        m.ext[:constraints][:hydrogen_storage_min_flow] = @constraint(m, [jh=JH,jy=JY], H_S_in[jh,jy] >= -Q_I * n_bundles ) # storage inflow is always above -Q_I 
        m.ext[:constraints][:hydrogen_storage_max_flow] = @constraint(m, [jh=JH,jy=JY], H_S_in[jh,jy] <= Q_I * n_bundles) # storage inflow is always below Q_I
        
        # hydrogen storage
        m.ext[:constraints][:hydrogen_storage_balance_1] = @constraint(m, [jh=JH[2:end],jy=JY], H_S[jh,jy] == H_S[jh-1,jy] + H_S_in[jh-1,jy]) 

        # the hydrogen flow at the last hour of the year is zero
        m.ext[:constraints][:hydrogen_storage_balance_2] = @constraint(m, [jy=JY], H_S_in[JH[end],jy] == 0) #

        # the storage starts with zero kg of hydrogen
        m.ext[:constraints][:hydrogen_storage_balance_3] = @constraint(m, H_S[1,1] == 0)

        # the storage on the first hour of the next year is equal to the storage on the last hour of the previous year
        m.ext[:constraints][:hydrogen_storage_balance_5] = @constraint(m, [jy=JY[2:end-1]], H_S[JH[end],jy] == H_S[1,jy+1])

    end    

    build_model!(m)

    optimize!(m)

    # timestamp the current run for saving the results
    function timestamp()
        return Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    end

    # other parameters for saving the results and for calculating the LCOH
    n_years = 10
    Timestamp = timestamp()
    Bundles = m.ext[:parameters][:n_bundles] 
    Profile = "TAP" 
    Volume = volume 
    capacity_solar = m.ext[:parameters][:Q_PV] # MW
    capacity_wind = m.ext[:parameters][:Q_W] # MW
    r = 0.05

    # Get the C_P and H_Pr values from the model, C_P is the cost of power and H_Pr is the hydrogen production (every hour)
    C_P = value.(m.ext[:variables][:C_P])
    H_Pr = value.(m.ext[:variables][:H_Pr])

    # Initialize arrays to hold total and weighted values, used for discounting and calculating LCOH
    CP_total = zeros(n_years)
    CP_weighted = zeros(n_years)
    HP_total = zeros(n_years)
    HP_weighted = zeros(n_years)
    CAPEX = zeros(n_years)
    CAPEX_weighted = zeros(n_years)
    OPEX = zeros(n_years)
    OPEX_weighted = zeros(n_years)
    TENNET = zeros(n_years)
    TENNET_weighted = zeros(n_years)
    HyNetwork = zeros(n_years)
    HyNetwork_weighted = zeros(n_years)

    # the values for the Cost components excepts power
    CAPEX[1] = m.ext[:parameters][:capacity_electrolyser] * investment_cost
    OPEX[1:10] .= m.ext[:parameters][:capacity_electrolyser] * operational_cost
    TENNET[1:10] .= m.ext[:parameters][:capacity_electrolyser] * TSO_cost_power
    HyNetwork[1:10] .= m.ext[:parameters][:capacity_electrolyser] * HyNetwork_cost 
    

    # Calculate total and weighted costs and productions for each year
    for year in 1:n_years
        CP_total[year] = sum(C_P[:, year])
        CP_weighted[year] = CP_total[year] / (1 + r)^(year - 1)

        CAPEX_weighted[year] = CAPEX[year] / (1 + r)^(year - 1)

        OPEX_weighted[year] = OPEX[year] / (1 + r)^(year - 1)

        TENNET_weighted[year] = TENNET[year] / (1 + r)^(year - 1)

        HyNetwork_weighted[year] = HyNetwork[year] / (1 + r)^(year - 1)
        
        HP_total[year] = sum(H_Pr[:, year])
        HP_weighted[year] = HP_total[year] / (1 + r)^(year - 1)
    end

    # Calculate the LCOH
    LCOH_P = sum(CP_weighted) / sum(HP_weighted)
    LCOH_CAPEX = sum(CAPEX_weighted) / sum(HP_weighted)
    LCOH_OPEX = sum(OPEX_weighted) / sum(HP_weighted)
    LCOH_TENNET = sum(TENNET_weighted) / sum(HP_weighted)
    LCOH_HyNetwork = sum(HyNetwork_weighted) / sum(HP_weighted)
    LCOH = LCOH_P + LCOH_CAPEX + LCOH_OPEX + LCOH_TENNET + LCOH_HyNetwork

    # read the results in a df, so we can add the new values
    df = CSV.read("Output\\LCOH.csv", DataFrame)

    # add the new values to a new row
    push!(df, (Timestamp, scenario, Profile, Bundles, Volume, capacity_solar, capacity_wind, LCOH_P, LCOH_CAPEX, LCOH_OPEX, LCOH_TENNET, LCOH_HyNetwork, LCOH))

    # save the values into a csv
    CSV.write("Output\\LCOH.csv", df)

end


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////
# defining specific run parameters and solving the model
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Defining specific run parameters - use ranges when studying combinations of inputs
# range(start = 0, stop = 0, step = 0)

n_bundles = 1 # n 
volume = 50 # % of capacity
capacity_solar = 200 # MW
capacity_wind = 200 # MW
scenarios = 2 # there are 6 scenarios in the data

# Create an iterator over all combinations, so we can solve the model for multiple inputs at the same time
combinations = IterTools.product(n_bundles, volume, capacity_solar, capacity_wind, scenarios)

# start process of solving for all combinations of input parameters
for (n, vol, solar, wind, scenario) in combinations
    complete_model(n, vol, solar, wind, scenario)
end

