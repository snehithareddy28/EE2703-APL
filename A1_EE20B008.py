from sys import argv, exit

# The program will show an error if there aren't precisely two  arguments in the commandline.
if len(argv)!=2:
	print("You shoould provide two arguements only.")
	exit()

# Assigning variables to .circuit and .end 
LINE1 = ".circuit"
LINE2 = ".end"

try:

# Opening the file mentioned in the commandline.
	with open(argv[1]) as f:
		list = f.readlines()

# These are parameters for veryfying file format errors. 		
		beg = -1; beg_check = -1; end = -2; end_check = -1 

# The programme will go through the file and extract only the parts that are required.
		for line in list:
			if LINE1 == line[:len(LINE1)]:
				beg = list.index(line)
				beg_check = 0;

			elif LINE2 == line[:len(LINE2)]:
				end = list.index(line)
				end_check = 0;

# The program will show an error if the circuit definition format is not proper.		
		if beg >= end or beg_check == -1 or end_check == -1:
			print("The circuit definition is not valid.")
			exit()

# These lines of code will reverse and print the required output.
		for list in reversed(list[beg+1:end]):
			a=reversed(list.split('#')[0].split())	
			a = "  ".join(a)
			print(a)

# The program will show this error if 
# the name of the netlist file is not proper. 
# or the netlist file is not found in the same directory as the program.
except FileNotFoundError:
	print("The File is not valid.")
	exit()