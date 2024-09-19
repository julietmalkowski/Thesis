import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#65 % of my data contains functional group classifications - can see this by printing number of empty rows

file_path = '/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/'

AD_metadata = np.load(f'{file_path}AD_metadata.npy', allow_pickle=True, fix_imports=True)
AS_metadata = np.load(f'{file_path}AD_metadata.npy', allow_pickle=True, fix_imports=True)

as_data_midas = np.load(f'{file_path}filtered_AS_data_midas.npy', allow_pickle=True, fix_imports=True)
ad_data_midas = np.load(f'{file_path}filtered_AD_data_midas.npy', allow_pickle=True, fix_imports=True)

midas_data = {'AS': as_data_midas, 'AD': ad_data_midas}
for key, filtered_AS_data_midas in midas_data.items():
    filtered_AS_data_midas = pd.DataFrame(filtered_AS_data_midas)    
    #make column names first row
    filtered_AS_data_midas.columns = filtered_AS_data_midas.iloc[0]
    filtered_AS_data_midas = filtered_AS_data_midas[1:]
    #groupby week and genus find mean of SRT_AS and unique value of all other columns
    filtered_AS_data_midas = filtered_AS_data_midas.groupby(['week', 'genus']).agg({(f'SRT_{key}'): 'mean', 'Chemoautotroph/mixotroph': 'unique', 'AOB': 'unique', 
                                                                                    'NOB': 'unique', 'Anammox': 'unique', 'Aerobic heterotroph': 'unique', 'PAO': 'unique',
                                                                                    'GAO': 'unique', 'Nitrite_reduction': 'unique', 'Sulfate_reduction': 'unique', 'Fermentation': 'unique',
                                                                                    'Acetogen': 'unique',  'Methanogen': 'unique'}).reset_index()

    #Plot the groups of genera for all weeks based on functional group
    filtered_AS_data_midas = filtered_AS_data_midas.astype({ 'Chemoautotroph/mixotroph': 'int', 'AOB': 'int', 
                                                    'NOB': 'int', 'Anammox': 'int', 'Aerobic heterotroph': 'int', 'PAO': 'int',
                                                    'GAO': 'int', 'Nitrite_reduction': 'int', 'Sulfate_reduction': 'int', 'Fermentation': 'int',
                                                    'Acetogen': 'int',  'Methanogen': 'int'})
    #turn into numpy array
    cols = filtered_AS_data_midas.columns
    as_midas_genus_per_week = filtered_AS_data_midas.to_numpy()
    cols = cols.to_numpy()
    vstacked = np.vstack((cols, as_midas_genus_per_week))
    np.save(f'{file_path}{key}_midas_genus_per_week', vstacked)
    
    growing = filtered_AS_data_midas[filtered_AS_data_midas['SRT_AS'] > 0]
    dying = filtered_AS_data_midas[filtered_AS_data_midas['SRT_AS'] < 0]
    # Sum each of the last 12 columns in growing and dying
    growing_sums = growing.iloc[:, -12:].sum()
    dying_sums = dying.iloc[:, -12:].sum()

    # Plot the sums as barplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    ax1.bar(growing_sums.index, growing_sums.values)
    ax1.set_xlabel('Functional Group')
    ax1.set_xticklabels(growing_sums.index, rotation=90)
    ax1.set_ylabel('Count')
    ax1.set_title(f'{key} Functional Group Counts in Active Population')

    ax2.bar(dying_sums.index, dying_sums.values)
    ax2.set_ylabel('Count')
    ax2.set_xlabel('Functional Group')
    ax2.set_xticklabels(dying_sums.index, rotation=90)
    ax2.set_title(f'{key} Functional Group Counts in Inactive Population')
    fig.tight_layout()
    fig.savefig(f'/Users/julietmalkowski/Desktop/Thesis/Final_Figures/{key}_Functional_Group_Counts.png', dpi=300)
    plt.show()


#count the number of each functional group when SRT < 0 and SRT > 0
filtered_AS_data_midas = np.load(f'{file_path}as_midas_genus_per_week.npy', allow_pickle=True, fix_imports=True)
filtered_AS_data_midas = pd.DataFrame(filtered_AS_data_midas)
filtered_AS_data_midas.columns = filtered_AS_data_midas.iloc[0]
filtered_AS_data_midas = filtered_AS_data_midas[1:]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# # Aerobic Heterotrophs
# bod_cod_load_removed = AS_metadata[:,[-1,9,11]]
# bod_cod_load_removed = bod_cod_load_removed[1:,:]
# ax1.plot(bod_cod_load_removed[:,0], bod_cod_load_removed[:,1], label='BOD/CBOD_Load_PE_Removed')
# ax1.set_xlim(25, 74)
# ax1.set_ylabel('Concentration')
# ax1.set_title('BOD/CBOD Removal over time in AS system')
# ax1.legend()

# filtered_AS_data_midas = filtered_AS_data_midas[filtered_AS_data_midas['Aerobic heterotroph'] == 1]
# filtered_AS_data_midas = filtered_AS_data_midas[['week', 'genus', 'SRT_AS']]
# for i in range(len(filtered_AS_data_midas['genus'].unique())):
#     ax2.plot(filtered_AS_data_midas[filtered_AS_data_midas['genus'] == filtered_AS_data_midas['genus'].unique()[i]]['week'],
#              filtered_AS_data_midas[filtered_AS_data_midas['genus'] == filtered_AS_data_midas['genus'].unique()[i]]['SRT_AS'])
# ax2.set_xlim(25, 74)
# ax2.set_xlabel('Time (weeks)')
# ax2.set_ylabel('SRT Time (d)')
# ax2.set_title('SRT over time in AS system for Aerobic heterotrophs')
# fig.savefig('/Users/julietmalkowski/Desktop/Thesis/Final_Figures/Aerobic_Heterotrophs_Over_Time.png', dpi=300)



## Nitrogen Reducers
# bod_cod_load_removed = AS_metadata[:,[-1,19,20]]
# bod_cod_load_removed = bod_cod_load_removed[1:,:]
# ax1.plot(bod_cod_load_removed[:,0], bod_cod_load_removed[:,1], label='Ammonia Removed')
# ax1.set_ylabel('Concentration')
# ax1.set_title('Ammonia Removed over time in AS system')
# ax1.legend()


# filtered_AS_data_midas = filtered_AS_data_midas[(filtered_AS_data_midas['Nitrite_reduction'] == 1)]
# filtered_AS_data_midas = filtered_AS_data_midas[['week', 'genus', 'SRT_AS']]
# for i in range(len(filtered_AS_data_midas['genus'].unique())):
#     ax2.plot(filtered_AS_data_midas[filtered_AS_data_midas['genus'] == filtered_AS_data_midas['genus'].unique()[i]]['week'],
#              filtered_AS_data_midas[filtered_AS_data_midas['genus'] == filtered_AS_data_midas['genus'].unique()[i]]['SRT_AS'])
#              #label=filtered_AS_data_midas['genus'].unique()[i])
# ax2.set_xlabel('Time (weeks)')
# ax2.set_ylabel('SRT')
# ax2.set_title('SRT over time in AS system for Nitrogen Reducers')
# fig.savefig('/Users/julietmalkowski/Desktop/Thesis/Final_Figures/Nitrogen_Reducers_Over_Time.png', dpi=300)


# Fermenters
filtered_AD_data_midas = np.load(f'{file_path}ad_midas_genus_per_week.npy', allow_pickle=True, fix_imports=True)
filtered_AD_data_midas = pd.DataFrame(filtered_AD_data_midas)
filtered_AD_data_midas.columns = filtered_AD_data_midas.iloc[0]
filtered_AD_data_midas = filtered_AD_data_midas[1:]
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))


bod_cod_load_removed = AD_metadata[:,[-1,7]]
bod_cod_load_removed = bod_cod_load_removed[1:,:]
ax1.plot(bod_cod_load_removed[:,0], bod_cod_load_removed[:,1], label='Methane Gas')
ax1.set_ylabel('Concentration')
ax1.set_title('Methane gas over time in AS system')
ax1.legend()


filtered_AD_data_midas = filtered_AD_data_midas[(filtered_AD_data_midas['Fermentation'] == 1)]
filtered_AD_data_midas = filtered_AD_data_midas[['week', 'genus', 'SRT_AS']]
#filtered_AD_data_midas['SRT_AS'] = filtered_AD_data_midas['SRT_AS'].apply(lambda x: 1/x if x != 0 else 0)

for i in range(len(filtered_AD_data_midas['genus'].unique())):
    ax2.plot(filtered_AD_data_midas[filtered_AD_data_midas['genus'] == filtered_AD_data_midas['genus'].unique()[i]]['week'],
             filtered_AD_data_midas[filtered_AD_data_midas['genus'] == filtered_AD_data_midas['genus'].unique()[i]]['SRT_AS'])
             #label=filtered_AS_data_midas['genus'].unique()[i])
ax2.set_xlabel('Time (weeks)')
ax2.set_ylabel('SRT')
ax2.set_title('SRT over time in AS system for Fermenters')
fig.savefig('/Users/julietmalkowski/Desktop/Thesis/Final_Figures/Fermenters_Over_Time.png', dpi=300)

