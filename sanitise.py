#!/usr/bin/env python

def sanitise(input_file):

    data_items = [] # Create list to store every line in input file as a list of words
    
    def prompt(): # Define a function prompt() that stores output filename
        output_name = input("Enter name for sanitised file\n") # Prompt user for output filename and assigned it to variable output_name
        try:
            f = open("%s.sdf" % output_name, "x") # Check if output file already exists        
        except:
            print("Error: Filename already in use")
            prompt() # Use recursion to keep calling prompt() until a valid output filename is entered
        return output_name
        
    output = prompt()
   

    with open (input_file, "r") as fi: # Read each line from input file
        for line in fi:
            words = line.split() # Splits the current line into a list where each word is a list item
            data_items.append(words) # Add each list of words into the end of data_items list
        
    row = 0
    for i in data_items:
        if (row == 3): # Finds the fourth row   
            k = 0
            for k in range (len(i)):
                if (k == 0): # Extract atom number
                    atom_num = int(i[k])
                    f = open("%s.sdf" % output, 'a')
                    f.write("\n {} ".format(str(i[k])))
                elif (k == 1): # Extract bond number
                    bond_num = int(i[k])
                    f = open("%s.sdf" % output, 'a')
                    f.write("%s  " % str(i[k]))
                elif (k == len(i) - 1):
                    i[k] = 1 # Replace last item with 1
                    f = open("%s.sdf" % output, 'a')
                    f.write("             {}\n".format(str(i[k])))
                elif (k > 1 and k < 6):
                    i[k] = 0 # Replace the third to second last items with zeroes
                    f = open("%s.sdf" % output, 'a')
                    f.write("%s  " % str(i[k]))
                    k = k+1     
            print("There are {} atoms with {} bonds in query".format(atom_num,bond_num))            
            
            
            
        elif (row > 3 and row < (atom_num + 4)): # For every line that contains atom information
            for k in range (len(i)):
                if (k < 3):
                    if (float(i[k]) == 0.0000 or float(i[k]) > 0.0000):
                        f = open("%s.sdf" % output, 'a')
                        f.write("    {}".format(str(i[k])))
                    else:
                        f = open("%s.sdf" % output, 'a')
                        f.write("   {}".format(str(i[k])))
                elif (k == 3):
                    f = open("%s.sdf" % output, 'a')
                    f.write(" {}  ".format(str(i[k])))
                elif (k==len(i)-1):
                    i[k] = 0 # Replace last item with zero
                    f = open("%s.sdf" % output, 'a')
                    f.write("{}\n".format(str(i[k])))
                else:
                    f = open("%s.sdf" % output, 'a')
                    f.write("{}  ".format(str(i[k])))
    
        else:
             for k in range (len(i)):
                if (k==len(i)-1):
                    f = open("%s.sdf" % output, 'a')
                    f.write("  {}\n".format(str(i[k])))
                else:
                    try:
                        if (int(i[k]) < 10):  
                            f = open("%s.sdf" % output, 'a')
                            f.write("  {}".format(str(i[k])))
                        else:
                            f = open("%s.sdf" % output, 'a')
                            f.write(" {}".format(str(i[k])))
                    except:
                        f = open("%s.sdf" % output, 'a')
                        f.write("{}".format(str(i[k])))
                     
        row = row + 1
         
    return ("%s.sdf" % output)      