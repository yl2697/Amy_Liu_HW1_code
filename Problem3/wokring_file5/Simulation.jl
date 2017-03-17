
function washout_simulation(time_start,time_stop,time_step_size,data_dictionary)

  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS;
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Phase 1: Run the model 1/4 of the final time w/o ATRA
  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_1 = 0.0;
  time_stop_phase_1 = 4.0;

  # Solve the model equations -
  (TP1,XP1) = SolveBalances(time_start_phase_1,time_stop_phase_1,time_step_size,data_dictionary);

  # Phase 2: Express gene_1 and run to end-time (ic = last point phase 1) -
  # Update the IC again -
  initial_condition_array = XP1[end,:];
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Grab the control parameters - turn on gene_1 =
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_gene_1_RNAP"] = 1.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_2 = time_stop_phase_1+time_step_size
  time_stop_phase_2 = time_start_phase_2 + 6.0

  # Solve the model equations -
  (TP2,XP2) = SolveBalances(time_start_phase_2,time_stop_phase_2,time_step_size,data_dictionary);

  # Washout -
  initial_condition_array = XP2[end,:];
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Grab the control parameters - turn on gene_1 =
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_gene_1_RNAP"] = 0.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

  time_start_phase_3 = time_stop_phase_2+time_step_size
  time_stop_phase_3 = time_stop
  (TP3,XP3) = SolveBalances(time_start_phase_3,time_stop_phase_3,time_step_size,data_dictionary);

  # Package the two phases together -
  T = [TP1 ; TP2 ; TP3];
  X = [XP1 ; XP2 ; XP3];

  return (T,X);
end

function adj_washout_inducer(time_start,time_stop,time_step_size,parameter_index,data_dictionary)
  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS;
  initial_condition_array=[initial_condition_array;zeros(number_of_states)]
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Phase 1: Run the model 1/4 of the final time w/o ATRA
  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_1 = 0.0;
  time_stop_phase_1 = 10.0;

  # Solve the model equations -
  (TP1,XP1) = SolveAdjBalances(time_start_phase_1,time_stop_phase_1,time_step_size,parameter_index,data_dictionary);

  # Phase 2: Express gene_1 and run to end-time (ic = last point phase 1) -
  # Update the IC again -
  initial_condition_array = XP1[end,:];
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Grab the control parameters - turn on gene_1 =
  control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  control_parameter_dictionary["W_gene_1_RNAP"] = 1.0;
  data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;

  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_2 = time_stop_phase_1+time_step_size
  time_stop_phase_2 = time_start_phase_2 + 1.0

  # Solve the model equations -
  (TP2,XP2) = SolveAdjBalances(time_start_phase_2,time_stop_phase_2,time_step_size,parameter_index,data_dictionary);

  # # Washout -
  # initial_condition_array = XP2[end,:];
  # data_dictionary["initial_condition_array"] = initial_condition_array;
  #
  # # Grab the control parameters - turn on gene_1 =
  # control_parameter_dictionary = data_dictionary["control_parameter_dictionary"]
  # control_parameter_dictionary["W_gene_1_RNAP"] = 0.0;
  # data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary;
  #
  # time_start_phase_3 = time_stop_phase_2+time_step_size
  # time_stop_phase_3 = time_stop
  # (TP3,XP3) = SolveBalances(time_start_phase_3,time_stop_phase_3,time_step_size,data_dictionary);
  #
  # # Package the two phases together -
  # T = [TP1 ; TP2 ; TP3];
  # X = [XP1 ; XP2 ; XP3];

  return (TP2,XP2);
end

function ss_repressed_adj_simulation(time_start,time_stop,time_step_size,parameter_index,data_dictionary)
  # First - run the model to steady-state w/no ATRA forcing -
  XSS = estimate_steady_state(0.001,data_dictionary)

  # Next, set the IC to the steady-state value -
  initial_condition_array = XSS;
  number_of_states = length(XSS)
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Phase 1: Run the model 1/4 of the final time w/o ATRA
  # Run the model for a section of time w/no ATRA forcing -
  time_start_phase_1 = 0.0;
  time_stop_phase_1 = 10.0;

  # Solve the model equations -
  (TP1,XP1) = SolveBalances(time_start_phase_1,time_stop_phase_1,time_step_size,data_dictionary);

  # Phase 2
  initial_condition_array = XP1[end,:];
  initial_condition_array = [initial_condition_array ; zeros(number_of_states)]
  data_dictionary["initial_condition_array"] = initial_condition_array;

  # Run the model for a section of time -
  time_start_phase_2 = time_stop_phase_1+time_step_size
  time_stop_phase_2 = time_start_phase_2 + 1.0

  # Solve the model equations -
  (TP2,XP2) = SolveAdjBalances(time_start_phase_2,time_stop_phase_2,time_step_size,parameter_index,data_dictionary);

  return (TP2,XP2);
end

function estimate_steady_state(epsilon,data_dictionary)

  initial_condition_vector = data_dictionary["initial_condition_array"];
  ic_array = copy(data_dictionary["initial_condition_array"])
  number_of_states = length(ic_array)

  # Setup loop -
  EPSILON = epsilon;
  TSTART = 0.0;
  Ts = 1.0;
  TSTOP = 1000;
  did_reach_steady_state = false
  while (!did_reach_steady_state)

    # solve the balances -
    (TSIM,X1) = SolveBalances(TSTART,TSTOP,Ts,data_dictionary)

    # Take a few additional steps -
    TNEXT_START = TSTOP+Ts;
    TNEXT_STOP = TNEXT_START+1.0;
    Ts = 0.1;

    # solve the balances again 0
    initial_condition_array = vec(X1[end,:])
    data_dictionary["initial_condition_array"] = initial_condition_array;
    (TSIM,X2) = SolveBalances(TNEXT_START,TNEXT_STOP,Ts,data_dictionary)

    # Find the difference -
    DIFF = norm((X2[end,:] - X1[end,:]));

    # Should we stop -or- go around again?
    if (DIFF<EPSILON)
      did_reach_steady_state = true;
      return (vec(X2[end,:]));
    else

      # No, we did *not* reach steady state ....
      TSTART = TSTOP+Ts
      TSTOP = 1.0 + TSTART;
      Ts = 0.1;

      initial_condition_array = vec(X2[end,:])
      data_dictionary["initial_condition_array"] = initial_condition_array;
    end
  end

  # return
  return XSS;
end

# run -
