import numpy as np
import pandas as pd
file_path = '/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/'

# Create numpy arrays for raw data- using cell_values_averaged.npy
# Data created: unfiltered_AS_data.npy, unfiltered_AS_weeks.npy, unfiltered_AD_data.npy, unfiltered_AD_weeks.npy
# region
def create_numpy_arrays_raw_data():
    values = np.load('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/cell_values_averaged.npy', allow_pickle=True, fix_imports=True)
    labels = values[0]
    ad_values = values[values[:,0].astype('U2') == 'AD']

    ad_values = ad_values[np.argsort([s[-4:]+s[-10:-8]+s[-7:-5] for s in ad_values[:,0].tolist()])]
    ad_values = np.vstack([labels, ad_values])
    #seperate the first column by '_'
    ad_weeks = np.array([s.split('_') for s in ad_values[1:,0]])
    ad_weeks = ad_weeks[:,1]
    #remove sample name and n column
    ad_values = ad_values[:,2:]

    

    #np array for all AS data
    as_values = values[values[:,0].astype('U2') == 'AS']
    as_values = as_values[np.argsort([s[-4:]+s[-10:-8]+s[-7:-5] for s in as_values[:,0].tolist()])]
    as_values = np.vstack([labels, as_values])
    as_weeks = np.array([s.split('_') for s in as_values[1:,0]])
    as_weeks = as_weeks[:,1]
    #remove sample name and n column
    as_values = as_values[:,2:]

    return ad_values, ad_weeks, as_values, as_weeks

ad_values, ad_weeks, as_values, as_weeks = create_numpy_arrays_raw_data()
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/unfiltered_AS_data.npy', as_values)
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/unfiltered_AS_weeks.npy', as_weeks)
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/unfiltered_AD_data.npy', ad_values)
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/unfiltered_AD_weeks.npy', ad_weeks)

#endregion


# Create a list of all OTU's growing in AD system
# Data created: filtered_growing_AD_otus.npy, filtered_dying_AD_otus.npy, filtered_growing_AS_otus.npy, filtered_dying_AS_otus.npy
# region
def create_list_of_growing_OTUs():
    filtered_AD_data = np.load(f'{file_path}filtered_AD_data_otu.npy', allow_pickle=True, fix_imports=True)

    #in each array, if the second value is greater than 0, then put the third value into a list
    AD_OTUs_growing = []
    AD_OTUs_dying = []
    #iterate through a numpy array with a list in each cell
    for i in range(filtered_AD_data.shape[0]):
        for j in range(len(filtered_AD_data[i])):
            if filtered_AD_data[i][j][1] > 0:
                AD_OTUs_growing.append(filtered_AD_data[i][j][2])
            if filtered_AD_data[i][j][1] < 0:
                AD_OTUs_dying.append(filtered_AD_data[i][j][2])
    AD_OTUs_growing = list(set(AD_OTUs_growing))
    AD_OTUs_dying = list(set(AD_OTUs_dying))

    filtered_AS_data = np.load(f'{file_path}filtered_AS_data_otu.npy', allow_pickle=True, fix_imports=True)

    #in each array, if the second value is greater than 0, then put the third value into a list
    AS_OTUs_growing = []
    AS_OTUs_dying = []
    #iterate through a numpy array with a list in each cell
    for i in range(filtered_AS_data.shape[0]):
        for j in range(len(filtered_AS_data[i])):
            if filtered_AS_data[i][j][1] > 0:
                AS_OTUs_growing.append(filtered_AS_data[i][j][2])
            if filtered_AS_data[i][j][1] < 0:
                AS_OTUs_dying.append(filtered_AS_data[i][j][2])
    AS_OTUs_growing = list(set(AS_OTUs_growing))
    AS_OTUs_dying = list(set(AS_OTUs_dying))
    
    return AD_OTUs_growing, AD_OTUs_dying, AS_OTUs_growing, AS_OTUs_dying

AD_OTUs_growing, AD_OTUs_dying, AS_OTUs_growing, AS_OTUs_dying = create_list_of_growing_OTUs()

np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_growing_AD_otus.npy', AD_OTUs_growing)
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_dying_AD_otus.npy', AD_OTUs_dying)

np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_growing_AS_otus.npy', AS_OTUs_growing)
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_dying_AS_otus.npy', AS_OTUs_dying)
#endregion

# For R data- all of these are exported to create phyloseq object (metadata from wwtp_parameters.py)
#region
#Data created: abundance_table.csv, taxonomy_table.xlsx
#region
abundance_table = np.load(f'{file_path}abundance_table.npy', allow_pickle=True, fix_imports=True)
#transposed abundance table
abundance_table = abundance_table.transpose()
#export to excel
np.savetxt('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/abundance_table.csv', abundance_table, delimiter=',', fmt='%s')

import pandas as pd
#taxonomy table
genus_table = np.load(f'{file_path}genus_table.npy', allow_pickle=True, fix_imports=True)
#turn to pd dataframe
genus_table = pd.DataFrame(genus_table)
genus_table = genus_table[['OTU','domain','phylum','class','order','family','genus','species']]
genus_table.to_excel('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/taxonomy_table.xlsx')
#endregion

# Create p-value otu table for R data
# AS Data created: filtered_as_data.xlsx
#region
#def create_numpy_arrays_growing_OTUs():
filtered_AS_data_otu = np.load(f'{file_path}filtered_AS_data_otu.npy', allow_pickle=True, fix_imports=True)

#mask for SRT_AS > 0
filtered_AS_data_growing = filtered_AS_data_otu[:,:,1] > 0
#mask for p-value > 0.0005
filtered_AS_data_high_p = filtered_AS_data_otu[:,:,3] > 0.0005

combined_mask = filtered_AS_data_growing & filtered_AS_data_high_p 

# Apply the combined mask to the array
filtered_array =  filtered_AS_data_otu[combined_mask]
filtered_array = pd.DataFrame(filtered_array)
filtered_array.columns = ['Week', 'SRT', 'OTU', 'p-value']
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_as_array.npy', filtered_array)

dates = np.load(f'{file_path}week-date-combo.npy', allow_pickle=True, fix_imports=True)
#filter for column 5 to only be 'AS'
as_dates = dates[:,[0,1,2]]
as_dates = pd.DataFrame(as_dates)
as_dates.columns = ['OTU','Date', 'Week']
#drop duplicates
as_dates = as_dates.drop_duplicates()

filtered_AS_data_growing = pd.merge(filtered_array, as_dates, on=['Week','OTU'], how='inner')

original = np.load(f'{file_path}abundance_table.npy', allow_pickle=True, fix_imports=True)
original = pd.DataFrame(original)
#make first row column names
original.columns = original.iloc[0]
#drop first row
original = original.drop(original.index[0])
#filter first coolumn to only be 'AS' in first two letters
original = original[original['Sample'].str[:2] == 'AS']
#keep only first column
original = original.iloc[:,0]
#split by '_'
original= original.str.split('_', expand=True)
original = original.iloc[:,1]
#drop duplicates
original = original.drop_duplicates()

#filter filtered_AS_data_growing to only include dates that are in original
filtered_AS_data_growing = filtered_AS_data_growing[filtered_AS_data_growing['Date'].isin(original)]

otu_table = filtered_AS_data_growing[['Date','OTU','p-value']]
#go from long to wide
otu_table = otu_table.pivot(index='OTU', columns='Date', values='p-value')
#fill missing values with 0
otu_table = otu_table.fillna(0)
#add 'AS_' to the beginning of each column
otu_table.columns = 'AS_' + otu_table.columns
otu_table = otu_table.to_excel('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/otu_table.xlsx')

as_data = np.load('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/as_metadata.npy', allow_pickle=True, fix_imports=True)
#put into pandas dataframe
as_data = pd.DataFrame(as_data)
#make first row column names
as_data.columns = as_data.iloc[0]
#remove first row
as_data = as_data[1:]
#add 'AS_' to the beginning of Date column
as_data['Date'] = 'AS_' + as_data['Date']
#Change 'Date' to 'Sample'
as_data = as_data.rename(columns={'Date':'Sample'})
#export to excel
as_data.to_excel('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/as_metadata.xlsx')

#endregion

# AD Data
# Data created: filtered_ad_data.xlsx
#region
filtered_AD_data_otu = np.load(f'{file_path}filtered_AD_data_otu.npy', allow_pickle=True, fix_imports=True)
#mask for SRT_AS > 0
filtered_AD_data_growing = filtered_AD_data_otu[:,:,0] > 0
#mask for p-value > 0.0005
filtered_AD_data_high_p = filtered_AD_data_otu[:,:,1] > 0.0005

combined_mask = filtered_AD_data_growing & filtered_AD_data_high_p 

# Apply the combined mask to the array
filtered_array =  filtered_AD_data_otu[combined_mask]
filtered_array = pd.DataFrame(filtered_array)
filtered_array.columns = ['Week', 'SRT', 'OTU', 'p-value']
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_ad_array.npy', filtered_array)

dates = np.load(f'{file_path}week-date-combo.npy', allow_pickle=True, fix_imports=True)
#filter for column 5 to only be 'AS'
ad_dates = dates[:,[0,1,2]]
ad_dates = pd.DataFrame(as_dates)
ad_dates.columns = ['OTU','Date', 'Week']
#drop duplicates
ad_dates = ad_dates.drop_duplicates()

filtered_AD_data_growing = pd.merge(filtered_array, ad_dates, on=['Week','OTU'], how='inner')

original = np.load(f'{file_path}abundance_table.npy', allow_pickle=True, fix_imports=True)
original = pd.DataFrame(original)
#make first row column names
original.columns = original.iloc[0]
#drop first row
original = original.drop(original.index[0])
#filter first coolumn to only be 'AS' in first two letters
original = original[original['Sample'].str[:2] == 'AD']
#keep only first column
original = original.iloc[:,0]
#split by '_'
original= original.str.split('_', expand=True)
original = original.iloc[:,1]
#drop duplicates
original = original.drop_duplicates()

#filter filtered_AS_data_growing to only include dates that are in original
filtered_AD_data_growing = filtered_AD_data_growing[filtered_AD_data_growing['Date'].isin(original)]

otu_table = filtered_AD_data_growing[['Date','OTU','p-value']]
#go from long to wide
otu_table = otu_table.pivot(index='OTU', columns='Date', values='p-value')
#fill missing values with 0
otu_table = otu_table.fillna(0)
#add 'AD_' to the beginning of each column
otu_table.columns = 'AD_' + otu_table.columns
otu_table = otu_table.to_excel('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/filtered_ad_otu_table.xlsx')

ad_data = np.load('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/ad_metadata.npy', allow_pickle=True, fix_imports=True)
#put into pandas dataframe
ad_data = pd.DataFrame(ad_data)
#make first row column names
ad_data.columns = ad_data.iloc[0]
#remove first row
ad_data = ad_data[1:]
#add 'AD_' to the beginning of Date column
ad_data['Date'] = 'AD_' + ad_data['Date']
#Change 'Date' to 'Sample'
ad_data = ad_data.rename(columns={'Date':'Sample'})
#export to excel
ad_data.to_excel('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Codes_R/ad_metadata.xlsx')
#endregion
#endregion

# create a numpy array with filtered AS Data with OTU's that appear in more than 35 weeks
# Data created: filtered_AS_data_over35_weeks_otu.npy
# region
filtered_as_data = np.load(f'{file_path}filtered_as_array.npy', allow_pickle=True, fix_imports=True)
filtered_as_data = pd.DataFrame(filtered_as_data)
filtered_as_data.columns = ['Week','SRT','OTU','p-value']
#filter out zotus that appear less than 35 times
super_filtered_as_data = filtered_as_data.groupby('OTU').filter(lambda x: len(x) > 35)
#turn into numpy array
super_filtered_as_data = super_filtered_as_data.to_numpy()
#vstack column names
super_filtered_as_data = np.vstack([filtered_as_data.columns, super_filtered_as_data])
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_AS_data_over35_weeks_otu.npy', super_filtered_as_data)

# endregion

# create a numpy array with filtered AD Data with OTU's that appear in more than 35 weeks
# Data created: filtered_AD_data_over35_weeks_otu.npy
# region
filtered_ad_data = np.load(f'{file_path}filtered_ad_array.npy', allow_pickle=True, fix_imports=True)
filtered_ad_data = pd.DataFrame(filtered_ad_data)
filtered_ad_data.columns = ['Week','SRT','OTU','p-value']
#filter out zotus that appear less than 35 times
super_filtered_ad_data = filtered_ad_data.groupby('OTU').filter(lambda x: len(x) > 35)
#turn into numpy array
super_filtered_ad_data = super_filtered_ad_data.to_numpy()
#vstack column names
super_filtered_ad_data = np.vstack([filtered_ad_data.columns, super_filtered_ad_data])
np.save('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/filtered_AD_data_over35_weeks_otu.npy', super_filtered_ad_data)

# endregion