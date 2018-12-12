from GM_function import *
import sys

model_code = sys.argv[1]

output = model_code + "/predicted_data.txt"
#print "DEBUG@@@", model_code, output


#param_file_lines = open("/scratch/berbellini/Inversion_Portugal/NA_inversion/NA_algorithm/src/rfi_subs/rfi_param.inc","r").readlines() 
#n_measure = int(param_file_lines[9].split()[3])

T0, E = gm_parallel_smart('modello.d', model_code)

out = open(output, "w")
for i in range(0,len(E)):
	out.write(str(T0[i])+'\t'+str(E[i])+'\n')
out.close()




