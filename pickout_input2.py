import pandas as pd
import re

path="/export/dahlefs/work/Metagenomes_chimneys_2020_workfolder/02_HUMAN_Decontam/"
#path=input("Give the full path to the sample files:")

# Define a function for adding lines to a file
def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a") as file_object:
        file_object.write(text_to_append)


# Read CSV
data = pd.read_csv("SampleInfoHakon.csv", sep='\t',header=(0))
data['Site'] = data['Site'].str.replace(' ','_')
data['SampleNames'] = data['SampleNames'].str.replace('_','-')

print("Available meta columns to choose from:")
print(data.columns)

# User inputs to gather data from + filename endings
input1=input("Which metadata column are we picking from?: ")

print((data[input1].unique()))

input2=input("Which value in that category to use?: ")
input2a=input("Any other metadata columns to pick within?(y/n): ")
if input2a=='y':
    print(data.columns)
    input2b=input("Which metadata column are we picking from?: ")
    print((data[input2b].unique()))
    input2c=input("Which value in that category to use?: ")

input3="-cleanR1.fq"
input4="-cleanR2.fq"
#input3=input("If there is a file type ending you want to make a list with, what is it? Otherwise, leave this blank: ")
#input4=input("If these are for paired files input the second file ending type here: ")

rprim=[]
cnt= 0	# set counter for moving through samples
for i in data.index:	# Indicate which row we're in using the dataframe index column
    sample = data.at[cnt,'SampleNames']	# designate which sample we are looking at
    category = data.at[cnt,input1]
    if input2a=='y':
        category2= data.at[cnt,input2b]
        if category==input2 and category2==input2c:	# Different data categories
            rprim.append(path+sample)
        elif input1==input2b:		# Same data category, different entries
            rprim.append(path+sample)
    elif category==input2:
        rprim.append(path+sample)	# append the sample name to the rprim list
    cnt=cnt+1
if input3=="" and input4=="":
    print("\n")
    to_string= (',').join(map(str, rprim))
    print(to_string)
elif input3 !="" and input4 =="":
    print("\n") # Make some vertical space  before the output
    to_string= (input3+',').join(map(str, rprim))	# Create a string of filenames from the list
    print(to_string+input3)
elif input3!="" and input4!="":
    print("\n") # Make some vertical space  before the output
    to_string= (input3+',').join(map(str, rprim)) 
    to_string2= (input4+',').join(map(str, rprim))      # Create a string of filenames from the list
    print("Writing line to file:")
    try:			# A second categorical input isn't always defined. Here we check if it exists.
        input2c
    except NameError:
         a=(input2)
    else:
         a=(input2+"/"+input2c)
    b=(to_string+input3)
    c=(to_string2+input4)
    print(a, b, c)
    line= "{}\t{}\t{}\n".format(a, b, c)
    append_new_line('megahitSamples.txt', line)
    
#    print("\n")
#    to_string= (input4+',').join(map(str, rprim))       # Create a string of filenames from the list
#    append_new_line('sampleR2s.txt', (to_string+input4))
#    print("Writing to R2 file:")
#    print(to_string+input4)
