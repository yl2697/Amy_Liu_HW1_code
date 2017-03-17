# This file is to be used as a tool to discriminate between models.

# Time:
time_end=251.0
time_inc=1.0
time_start=1.0


# Extract dat files:
path_to_response1_file="/Users/amyyiwenliu/Desktop/Cornell/Spring 2017/CHEME 7770/HW/Problem 3/working_file3/sensitivity"
x1=readdlm(path_to_response1_file)
path_to_response2_file="/Users/amyyiwenliu/Desktop/Cornell/Spring 2017/CHEME 7770/HW/Problem 3/working_file4/sensitivity"
x2=readdlm(path_to_response2_file)

result=zeros(time_end,1)

for i = time_start:time_inc:time_end
  delta_x=x1[i]-x2[i];
  trans_delta_x=delta_x';
  weight=1.0/((x1[i]+x2[i])/2.0)^2
  func=trans_delta_x*weight*delta_x
  result[i]=func;
end


area=0;

for i =time_start:time_inc:(time_end-1)
  height=i+1-i
  top=result[i]
  bottom=result[i+1]
  area_new=(top+bottom)*height/2.0
  area=area+area_new
  maximum_area=max(area)
end
