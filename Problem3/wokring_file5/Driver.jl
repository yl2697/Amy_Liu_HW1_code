# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #

# include -
include("Include.jl")
using PyPlot

# define some colors -
const shaded_color_value = "lightgray"
const mean_color_value = "dimgray"
const experimental_color_value = "black"

const P1_color = "blue"
const P1_shaded_color="skyblue"

const P2_color = "orange"
const P2_shaded_color="navajowhite"

const P3_color = "green"
const P3_shaded_color="lightgreen"

# Script to solve the balance equations -
time_start = 0.0
time_stop = 35.0
time_step_size = 0.01
number_of_timesteps = length(time_start:time_step_size:time_stop)

# Load the data dictionary (default parameter values) -
data_dictionary = DataDictionary(time_start,time_stop,time_step_size)

# # How many samples do we want to explore?
# number_of_samples = 10
# sigma = 0.20

#data_dictionary = deepcopy(data_dictionary)

(T,X) = washout_simulation(time_start,time_stop,time_step_size,data_dictionary)

# What is the size of the system?
number_of_states = data_dictionary["number_of_states"]

state_array = X[:,1:number_of_states]

file_path = "./sensitivity/state_array.dat"
writedlm(file_path,state_array)

plot(T,X[:,7],label="protein_gene_1")
plot(T,X[:,8],label="protein_gene_2")
plot(T,X[:,9],label="protein_gene_3")
xlabel("Time (hr)")
ylabel("Protein Expression (concentration)")
legend()
title("Problem 3- Model 2")
# # main loop -
# parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
# average_scaled_sensitivity_array = zeros(number_of_states,1)
# for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)
#
#   @show parameter_index
#
#   # grab the dictionary -
#   data_dictionary = deepcopy(data_dictionary)
#
#   # Solve the adj simulation -
#   # You need to point this to your specific adj simulation code -
#   #
#   # e.g.,
#   (T,X) = washout_simulation(time_start,time_stop,time_step_size,data_dictionary)
#
#   # split -
#     state_array = X[:,1:number_of_states]
#
#     ## time average -
#     #average_sensitivity_col = time_average_array(T,scaled_sensitivity_array)
#
#     # grab -
#     #average_scaled_sensitivity_array = [average_scaled_sensitivity_array average_sensitivity_col]
#     #data_array = [T X]
#     file_path = "./sensitivity_ss/AdjSimulation-P"*string(parameter_index)*".dat"
#     writedlm(file_path,state_array)
#   end
#
#
# # Load the data dictionary -
# #data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
#
# # Solve the model equations -
# #(T,X) = SolveBalances(time_start,time_stop,time_step_size,data_dictionary)
