include("Include.jl")

data_dictionary=DataDictionary(0.0,0.0,0.0)

path_to_sensitivity_files="./sensitivity"
file_pattern="AdjSimulation-P"
time_skip=20;


(T,SA)=calculate_sensitivity_array(path_to_sensitivity_files,file_pattern,time_skip,data_dictionary)
ts=length(T)

SA_mRNA=zeros(ts,24)
for i=1:ts
  SA_mRNA[i,:]=SA[[4+9*(i-1)],:]
end

SA_P3=zeros(ts,24)
for i=1:ts
  SA_P3[i,:]=SA[[9+9*(i-1)],:]
end

SA_simplified=[SA_mRNA;SA_P3]

IP=estimate_identifiable_parameters(SA_simplified,0.1)
