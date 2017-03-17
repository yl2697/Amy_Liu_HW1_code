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

function createArr_steady()
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

  # Setup the timescale of the simulation -
  time_start = 0.0
  time_stop = 10.0
  time_step_size = 0.1
  number_of_timesteps = length(time_start:time_step_size:time_stop)

  # Load the data dictionary (default parameter values) -
  data_dictionary = DataDictionary(time_start,time_stop,time_step_size)
  number_of_states = data_dictionary["number_of_states"]
  # How many samples do we want to explore?
  number_of_samples = 10
  sigma = 0.20

  # What is the size of the system?
  #number_of_states = data_dictionary["number_of_states"]

  # main loop -
  parameter_name_mapping_array = data_dictionary["parameter_name_mapping_array"]
  average_scaled_sensitivity_array2 = zeros(number_of_states,1)
  for (parameter_index,parameter_value) in enumerate(parameter_name_mapping_array)

    @show parameter_index

    # grab the dictionary -
    local_data_dictionary = deepcopy(data_dictionary)

    # Solve the adj simulation -
    # You need to point this to your specific adj simulation code -
    #
    # e.g.,
    (T,X) = adj_washout_inducer_steady_state(time_start,time_stop,time_step_size,parameter_index,local_data_dictionary)

    # split -
      state_array = X[:,1:number_of_states]
      sensitivity_array = X[:,number_of_states+1:end]
      scaled_sensitivity_array = scale_sensitivity_array(T,state_array,sensitivity_array,parameter_index,local_data_dictionary)

      # time average -
      average_sensitivity_col = time_average_array(T,scaled_sensitivity_array)

      # grab -
      average_scaled_sensitivity_array2 = [average_scaled_sensitivity_array2 average_sensitivity_col]
      data_array = [T X]
      file_path = "./sensitivity_ss/AdjSimulation-P"*string(parameter_index)*".dat"
      writedlm(file_path,data_array)
    end
      return average_scaled_sensitivity_array2[:,2:end]
end

function makegraph(average_scaled_sensitivity_array)
  parameter_control = svd(average_scaled_sensitivity_array)
  #parameter_control2 = svd(average_scaled_sensitivity_array)
  #file_path2 = "./sensitivity/SVD.dat"
  #file_path3 = ".sensitivity_ss/SVD.dat"

  #writedlm(file_path2,parameter_control)
  #writedlm(file_path3,parameter_control)
  # x=collect(0:size(average_scaled_sensitivity_array,1))
  # y=collect(0:size(average_scaled_sensitivity_array,1))
  abs_avg=abs(average_scaled_sensitivity_array)
  #abs_avg2=abs(average_scaled_sensitivity_array)

  pcolormesh(abs_avg, cmap="OrRd")
  #pcolormesh(average_scaled_sensitivity_array,cmap="OrRd")
  colorbar()
  ax=gca()
  xlabel("Parameters")
  ylabel("Species")
  title("Probelm 2- Sensitivity Array for Steady State")
  ax[:xaxis][:set_ticks](x-.5)
  ax[:xaxis][:set_ticks](labels,rotation=60, fontsize=8)
  ax[:yaxis][:set_ticks](x-.5)
  ax[:yaxis][:set_ticks](labels,rotation=0, fontsize=8)

  # pcolormesh(abs_avg2,cmap="OrRd")
  # colorbar()
  # ax=gca()
  # labels=data_dictionary["parameter_name_mapping_array"]
  # ax[:xaxis][:set_ticks](x-.5)
  # ax[:xaxis][:set_ticks](labels,rotation=60, fontsize=8)
  # ax[:yaxis][:set_ticks](x-.5)
  # ax[:yaxis][:set_ticks](labels,rotation=0, fontsize=8)

end
