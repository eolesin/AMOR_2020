import pandas as pd
import re

# We have some differences in the Excel file names and the actual seq file names
seqinfo_df = pd.read_csv("SeqInfo_EDedit.csv", sep='\t',header=(0), encoding='latin-1')
seqinfo_df['SampleNames'] = seqinfo_df['SampleNames'].str.replace('_','-')

# This was also true for the descriptions file
descr_df = pd.read_csv("SampleInfoHakon.csv", sep='\t',header=(0), encoding='latin-1')
descr_df['SampleNames'] = descr_df['SampleNames'].str.replace('_','-')

# Outer merge these on the SampleNames column
r= pd.merge(seqinfo_df,descr_df, on='SampleNames', how='outer')

# Bringing in the human contamination seq file size ls -l outputs
human_df = pd.read_csv("humanfilesizes.csv", sep='\t',header=(0))

rprim=[] # empty dataframe to store complete data in 
cnt= 0	# set counter for moving through samples
for i in r.index:	# Indicate which row we're in using the dataframe index column
        sample = r.at[cnt,'SampleNames']	# designate which sample we are looking at
        hum_cnt=0				# start the counter for the human filesize file
        for x in human_df.index:
                line = human_df.at[hum_cnt,'filename']
                size= human_df.at[hum_cnt,'filesize']
                if sample in line:
                        rprim.append([sample, size])	# append the information into the rprim dataframe
                hum_cnt=hum_cnt+1
        cnt=cnt+1
rprim_df=pd.DataFrame(rprim, columns= ["SampleNames", "HumanContamFileSize"]) # Create a dataframe from the loop output, give it headers.

rfinal= pd.merge(r, rprim_df, on='SampleNames', how='outer') # Merge the file with the previously merged info
rfinal.to_csv('SampleMaster.csv', sep="\t", index=False)
