import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# filename = path.join(mkdtemp(), 'newfile.dat')
# fp = np.memmap(filename, dtype='float32', mode='w+', shape=(3,4))

#Load data
srt_data = np.load('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/genus_table.npy', allow_pickle=True, fix_imports=True)
######################################################### Removing outliers from each week ##############################################################
class Growth_Groups:
    def __init__(self, df):
        self.df = df
        self.SRT_AS, self.AS_outliers, self.SRT_AD, self.AD_outliers = self.create_SRT_dataframe()

    def create_SRT_dataframe(self):
        #print(self.df.shape)
        df_as = self.df[:,[0,21,27]]
        df_ad = self.df[:,[0,20,27]]
        print(df_as)
        dfs = {'AS system' : df_as, 'AD system' : df_ad}

        for key, df_new in dfs.items():
            number_of_weeks = np.max(df_new[:,0].astype(int)) - np.min(df_new[:,0].astype(int))
            df_final = np.empty((number_of_weeks+1, 18155, 3), dtype=object)
            outliers_final = {}

            for i, week in enumerate(range(np.min(df_new[:,0].astype(int)), np.max(df_new[:,0].astype(int))+1)):
                week_data = df_new[df_new[:,0].astype(int) == week]
                week_data_to_filter = week_data[week_data[:,1].astype(float) != 0]

                try:
                    q1 = np.percentile(week_data_to_filter[:,1].astype(float), 25)
                    q3 = np.percentile(week_data_to_filter[:,1].astype(float), 75)
                    iqr = q3 - q1
                    lower_bound = q1 - 1.5 * iqr
                    upper_bound = q3 + 1.5 * iqr

                    outliers = week_data[(week_data[:,1].astype(float) < lower_bound) | (week_data[:,1].astype(float) > upper_bound)]
                    
                    week_data[(week_data[:,1].astype(float) < lower_bound) | (week_data[:,1].astype(float) > upper_bound)] = np.nan
                    outliers_final[week] = outliers[:,2]

                except:
                    print(f'week missing data: {week}')

                df_final[i] = week_data
                outliers_final[week] = outliers[:,2]
                # Plot boxplot for the week
                week_srt_times = week_data[:,1].astype(float)
                week_srt_times = week_srt_times[~np.isnan(week_srt_times)]
        
            outliers = list(outliers_final.values())
            # Concatenate all arrays into one array
            outliers = np.concatenate(outliers)

            #remove nan values
            outliers = outliers[~pd.isnull(outliers)]
            unique_strings, counts = np.unique(outliers, return_counts=True)

            # Sort the counts in descending order and get the indices
            sorted_indices = np.argsort(counts)[::-1]

            # Select the top 20 unique string values and their counts
            top_20_indices = sorted_indices[:20]
            top_20_strings = unique_strings[top_20_indices]
            top_20_counts = counts[top_20_indices]

            # Plot the frequencies of the top 20 unique string values
            plt.figure(figsize=(10, 6))
            plt.bar(top_20_strings, top_20_counts)
            plt.xlabel('Genus')
            plt.ylabel('Frequency')
            plt.title(f'Top 20 Genera with the Most Frequent Outliers in {key}')
            plt.xticks(rotation=90)
            plt.tight_layout()
            #plt.show()
            plt.savefig(f'/Users/julietmalkowski/Desktop/Thesis/top_20_genera_outliers_{key[0:2]}.png')
            plt.close()

            if key == 'AS system':
                df_final_AS = df_final
                outliers_final_AS = outliers_final
            else:
                df_final_AD = df_final
                outliers_final_AD = outliers_final
        
        return df_final_AS, outliers_final_AS, df_final_AD, outliers_final_AD
    
    def create_srt_graphs(self):
        #Aesthetics:
        sns.set_style("ticks")
        sns.color_palette("Set2")

        df_as = self.df[:,[0,21,27]]
        df_ad = self.df[:,[0,20,27]]
        df_as = pd.DataFrame(df_as, columns = ['week', 'SRT', 'Genus'])
        df_as['week'] = df_as['week'].astype(int)
        df_as['SRT'] = df_as['SRT'].astype(float)
        df_ad = pd.DataFrame(df_ad, columns = ['week', 'SRT', 'Genus'])
        df_ad['week'] = df_ad['week'].astype(int)
        df_ad['SRT'] = df_ad['SRT'].astype(float)
        dfs = {'AS System' : df_as, 'AD System' : df_ad}
        
        plt.figure(figsize=(12, 10)) 
        for key, df_new in dfs.items():
            if key == 'AS System':
                ylims_as = (-10, 10)
                #as_plot = sns.violinplot(data = df_new , x = 'week', y = 'SRT', inner=None, zorder=2, linewidth=1, scale='width', palette= ['lightskyblue'])
                as_plot = sns.stripplot(data = df_new , x = 'week', y = 'SRT', jitter=True, zorder=1, size = 1, color= 'orange')
                as_plot.set(ylim = ylims_as)
            
            else:
                plt.figure(figsize=(12, 10)) 
                ylims_ad = (-100, 100)
                #ad_plot = sns.violinplot(data = df_new , x = 'week', y = 'SRT', inner=None, zorder=2, linewidth=1, scale='width', palette= ['lightskyblue'])
                ad_plot = sns.stripplot(data = df_new , x = 'week', y = 'SRT', jitter=True, zorder=1, size = 1, color= 'orange')
                ad_plot.set(ylim = ylims_ad)
              
            plt.title(f'Weekly SRT of {key}')
            plt.xlabel("week", fontsize=15)
            plt.ylabel("SRT Time [d]", fontsize=15)
            plt.tight_layout()
            plt.show()
            plt.savefig(f'/Users/julietmalkowski/Desktop/Thesis/srt_graph_{key[0:2]}.png', dpi = 300)
            plt.close()

    def create_growth_rate_graphs(self):
        #Aesthetics:
        sns.set_style("ticks")
        sns.color_palette("Set2")

        df_as = self.df[:,[0,21,27]]
        df_ad = self.df[:,[0,20,27]]
        df_as = pd.DataFrame(df_as, columns = ['week', 'SRT', 'Genus'])
        df_as['week'] = df_as['week'].astype(int)
        df_as['SRT'] = df_as['SRT'].astype(float)
        df_as['Growth Rate'] = df_as['SRT'].apply(lambda x: 1/x if x != 0 else 0)
        #remove SRT column
        df_as = df_as.drop(columns = ['SRT'])
        df_ad = pd.DataFrame(df_ad, columns = ['week', 'SRT', 'Genus'])
        df_ad['week'] = df_ad['week'].astype(int)
        df_ad['SRT'] = df_ad['SRT'].astype(float)
        df_ad['Growth Rate'] = df_ad['SRT'].apply(lambda x: 1/x if x != 0 else 0)
        #remove SRT column
        df_ad = df_ad.drop(columns = ['SRT'])
        dfs = {'AS System' : df_as, 'AD System' : df_ad}
        
        plt.figure(figsize=(12, 10)) 
        for key, df_new in dfs.items():
            if key == 'AS System':
                ylims_as = (-100, 100)
                as_plot = sns.violinplot(data = df_new , x = 'week', y = 'Growth Rate', inner=None, zorder=2, linewidth=1, scale='width', palette= ['lightskyblue'])
                as_plot = sns.stripplot(data = df_new , x = 'week', y = 'Growth Rate', jitter=True, zorder=1, size = 1, color= 'orange')
                as_plot.set(ylim = ylims_as)
            
            else:
                plt.figure(figsize=(12, 10)) 
                ylims_ad = (-15, 5)
                ad_plot = sns.violinplot(data = df_new , x = 'week', y = 'Growth Rate', inner=None, zorder=2, linewidth=1, scale='width', palette= ['lightskyblue'])
                ad_plot = sns.stripplot(data = df_new , x = 'week', y = 'Growth Rate', jitter=True, zorder=1, size = 1, color= 'orange')
                ad_plot.set(ylim = ylims_ad)
              
            plt.title(f'Weekly Growth Rate of {key}')
            plt.xlabel("week", fontsize=15)
            plt.ylabel("Growth Rate [1/d]", fontsize=15)
            plt.tight_layout()
            #plt.show()
            plt.savefig(f'/Users/julietmalkowski/Desktop/Thesis/Final_Figures/growth_rate_graph_{key[0:2]}.png', dpi = 300)
            plt.close()
filtered_data = Growth_Groups(srt_data)
filtered_AS_data = filtered_data.SRT_AS
np.save('/Users/julietmalkowski/Desktop/Thesis/filtered_AS_data.npy', filtered_AS_data)
filtered_AD_data = filtered_data.SRT_AD
np.save('/Users/julietmalkowski/Desktop/Thesis/filtered_AD_data.npy', filtered_AD_data)
graphs = filtered_data.create_srt_graphs()
graphs = filtered_data.create_growth_rate_graphs()




## OTU Version
asv_srt_table = np.load('/Users/julietmalkowski/Desktop/Thesis/Final_Codes/Final_Objects/asv_srt_table.npy', allow_pickle=True, fix_imports=True)
def create_SRT_dataframe(self):
    #print(self.df.shape)
    df_as = self[1:,[0,21,1,15]]
    df_ad = self[1:,[0,20,1,13]]
    print(df_as)
    dfs = {'AS system' : df_as, 'AD system' : df_ad}

    for key, df_new in dfs.items():
        number_of_weeks = np.max(df_new[:,0].astype(int)) - np.min(df_new[:,0].astype(int))
        df_final = np.empty((number_of_weeks+1, 18155, 4), dtype=object)
        outliers_final = {}

        for i, week in enumerate(range(np.min(df_new[:,0].astype(int)), np.max(df_new[:,0].astype(int))+1)):
            week_data = df_new[df_new[:,0].astype(int) == week]
            week_data_to_filter = week_data[week_data[:,1].astype(float) != 0]

            try:
                q1 = np.percentile(week_data_to_filter[:,1].astype(float), 25)
                q3 = np.percentile(week_data_to_filter[:,1].astype(float), 75)
                iqr = q3 - q1
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr

                outliers = week_data[(week_data[:,1].astype(float) < lower_bound) | (week_data[:,1].astype(float) > upper_bound)]
                
                week_data[(week_data[:,1].astype(float) < lower_bound) | (week_data[:,1].astype(float) > upper_bound)] = np.nan
                outliers_final[week] = outliers[:,2]

            except:
                print(f'week missing data: {week}')

            df_final[i] = week_data
            outliers_final[week] = outliers[:,2]
            # Plot boxplot for the week
            week_srt_times = week_data[:,1].astype(float)
            week_srt_times = week_srt_times[~np.isnan(week_srt_times)]
    
        outliers = list(outliers_final.values())
        # Concatenate all arrays into one array
        outliers = np.concatenate(outliers)

        #remove nan values
        outliers = outliers[~pd.isnull(outliers)]
        unique_strings, counts = np.unique(outliers, return_counts=True)

        # Sort the counts in descending order and get the indices
        sorted_indices = np.argsort(counts)[::-1]

        # Select the top 20 unique string values and their counts
        top_20_indices = sorted_indices[:20]

        if key == 'AS system':
            df_final_AS = df_final
            outliers_final_AS = outliers_final
        else:
            df_final_AD = df_final
            outliers_final_AD = outliers_final
    
    return df_final_AS, outliers_final_AS, df_final_AD, outliers_final_AD

final_AS, outliers_final_AS, final_AD, outliers_final_AD = create_SRT_dataframe(asv_srt_table)

np.save('/Users/julietmalkowski/Desktop/Thesis/filtered_AS_data_otu.npy', final_AS)
np.save('/Users/julietmalkowski/Desktop/Thesis/filtered_AD_data_otu.npy', final_AD)
