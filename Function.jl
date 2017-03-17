using ODE
using PyPlot

function balances(t,x)
# Define parameters when looking at differnce between open and close complex.

mj_close = x[1]
mj_open = x[2]

# # Define paramters when looking at difference between different Ltj.
# mj_open1 = x[1]
# mj_open2 = x[2]
# mj_open3 = x[3]

  # Default: HL-60, we've updated these numbers from bionumbers -
  	# ------------------------------------------------------------------------------------------#
  	# constants (from bionumbers)       units
  	# ------------------------------------------------------------------------------------------#
  	cell_diammax_transcription_rateer = 1.1                 # mum
  	number_of_rnapII = 4600            	# copies/cells
  	mRNA_half_life_TF = 0.083           # hrs
    doubling_time_cell = 1000           # hrs
  	max_transcription_rate = 60.0       # nt/sec
  	average_transcript_length = 1200   	# nt
  	fraction_nucleus = 0.0             	# dimensionless
  	av_number = 6.02e23                 # number/mol
  	avg_gene_number = 2                 # number of copies of a gene
  	polysome_number = 4								# number of ribsomoses per transcript
  	# ------------------------------------------------------------------------------------------#
  	#
  	# ------------------------------------------------------------------------------------------#
  	# Calculate constants using bionumber values
  	# ------------------------------------------------------------------------------------------#
  	# Calculate the volume (convert to L)
  	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diammax_transcription_rateer)^3)*(1e-15)

  	# Calculate the rnapII_concentration and ribosome_concentration
  	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM

  	# degrdation rate constants -
  	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1

  	# kcats for transcription and translation -
  	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1

  	# Maximum specific growth rate (death_rate_constantmaximum_specific_growth_rate) -
  	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1

  	# What is the average gene concentration -
  	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM

  	# How fast do my cells die (kd,j)?
  	death_rate_constant = 0.2*maximum_specific_growth_rate                            # hr^-1

  	# Saturation constants for translation and trascription -
  	saturation_transcription = 600*(1/av_number)*(1/V)*1e9                           # nM
  	# -------------------------------------------------------------------------------------------#

  Lt=100.0; # characteristic gene lenght (nt)
  Ltj=100.0;
  # Ltj_1=10.0;
  # Ltj_2=100.0;
  # Ltj_3=1000.0;
  kt=max_transcription_rate/Lt;
  w_RNAP=49.0
  #Rt_open=5.0*1000; # RNAP aboundance [nM] dRt_open/dt=bo-b'Rt_open

  #Pt=2.4; # promoter concentration (from the paper)
  #R=90.0; # enzyme concentration (from the paper)
  #k_inv=2.0*10.0^(-2); #(made up number)
  #k2=4.0*10.0^(-2); # [sec^-1](obtained from the paper)
  #k1=9.0*10.0^(-3)*k_inv/k2; # use the eq from table 3
  #kobs=k1*R*k2/(k1*R+k_inv+k2);
  #Rt_close=Pt*(1.0-exp(-kobs*t));

  # uj- regulatory term:

  uj_open=w_RNAP/(1.0+w_RNAP);
  uj_close=(1.0-1.0/(1.0+exp(t-2.0)))*uj_open;

  #Transcription rate for Open/Close complexes.
  kinmax_transcription_rateic_rate = kt*number_of_rnapII*Lt/Ltj;

  # # Transcription rate with different Ltj.
  # kinmax_transcription_rateic_rate1 = kt*number_of_rnapII*Lt/Ltj_1;
  # kinmax_transcription_rateic_rate2 = kt*number_of_rnapII*Lt/Ltj_2;
  # kinmax_transcription_rateic_rate3 = kt*number_of_rnapII*Lt/Ltj_3;

  #rT,j in terms of open and close complex
  transcription_rate_open = kinmax_transcription_rateic_rate*uj_open;
  transcription_rate_close = kinmax_transcription_rateic_rate*uj_close;

  # #rT,j in terms of different Lt,j:
  # transcription_rate_open1 = kinmax_transcription_rateic_rate1*uj_open
  # transcription_rate_open2 = kinmax_transcription_rateic_rate2*uj_open
  # transcription_rate_open3 = kinmax_transcription_rateic_rate3*uj_open


#dmjdt for Open and Close Complex
  dmjdt_open = transcription_rate_open-(death_rate_constant+maximum_specific_growth_rate)*mj_open;
  dmjdt_close = transcription_rate_close-(death_rate_constant+maximum_specific_growth_rate)*mj_close;
  dxdt = [dmjdt_close; dmjdt_open]

# # dmjdt for different Ltj:
#   dmjdt_open1 = transcription_rate_open1-(death_rate_constant+maximum_specific_growth_rate)*mj_open1;
#   dmjdt_open2 = transcription_rate_open2-(death_rate_constant+maximum_specific_growth_rate)*mj_open2
#   dmjdt_open3 = transcription_rate_open3-(death_rate_constant+maximum_specific_growth_rate)*mj_open3
#   dxdt=[dmjdt_open1;dmjdt_open2;dmjdt_open3]
  return dxdt;
end

function solver()
  time_start = 0.0
  time_stop = 10.0 # [min]
  time_step_size = 0.5 # [min]
  tspan = collect(time_start:time_step_size:time_stop)
  initial_condition=[0.0,0.0]  # for open/close complex
  #initial_condition=[0.0,0.0,0.0] # for different Ltj

  (t,X)=ode45(balances,initial_condition,tspan)

  # Ploting for Open/Close Complex:
  plot(t,[a[1] for a in X], label="Close Complex")
  plot(t,[a[2] for a in X], label="Open Complex")
  #legend(bbox_to_anchor=(1.05,1) loc=2, borderaxespad=0.)
  axis([0,10,0,25000])
  xlabel("time")
  ylabel("mRNA lavel")
  title("Problem 1")
  legend()

  # # Ploting for differnet Ltj:
  # plot(t,[a[1] for a in X], label="Ltj=10.0")
  # plot(t,[a[2] for a in X], label="Ltj=100.0")
  # plot(t,[a[3] for a in X], label="Ltj=1000.0")
  # #legend(bbox_to_anchor=(1.05,1) loc=2, borderaxespad=0.)
  # axis([0,10,0,250000])
  # xlabel("time")
  # ylabel("mRNA lavel")
  # title("Problem 1")
  # legend()

end
