# imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from chembl_webresource_client.new_client import new_client
from rdkit.Chem import Lipinski, Descriptors
from rdkit import Chem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import r2_score, mean_squared_error
# import lazypredict
from collections import Counter
# from lazypredict.Supervised import LazyRegressor
import  os
from sklearn.tree import DecisionTreeRegressor

def model_scoring(y_test, y_pred):
    """Calculates the statistics necessary to score a ML models performance

    Parameters
    ----------
    y_test : PD_array
        Test values split from data using sklearn.train_test_split
    y_pred : PD_array
        Predicted values generated from model.predict(X_test)

    Returns
    -------
    Statistics : PD_series
        Array of values coorisponding to MSE, R2, ASE

    """

    from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

    MSE = float(round(mean_squared_error(y_test, y_pred),2))
    R2 = float(round(r2_score(y_test, y_pred),2))
    MAE = float(round(mean_absolute_error(y_test, y_pred),2))

    Statistics = [{'MSE' : MSE}, {'R2' : R2}, {'MAE' : MAE}]

    return Statistics

def lipinski(smiles, verbose=False):
    """Determines the lipinski descriptors for each molecule and creates a list
    containing each of the descriptors

    Parameters
    ----------
    smiles : 1D_Array
    verbose :
         (Default value = False)

    Returns
    -------
    Descriptors: PD_array
            Array of values MolWt, MolLogP, NumHDonors, NumHAcceptors coorisponding to each
        cannonical smile in the smiles array
        
        MolWt: int
        MolLogP: int
        NumHDonors: int
        NumHAcceptors: int

    """
    moldata=[]
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem)
        moldata.append(mol)
    
    baseData = []
    i=0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
    
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHAcceptors,
                        desc_NumHDonors])
        
        if (i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1
    columnNames=["MW", "LogP", "NumDonors", "NumAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors

def pIC50(input):
    """Takes normalized standard values and converts them to the pIC50 for chemical analysis
    further on in the program

    Parameters
    ----------
    input : PD_dataframe
        PD_dataframe with column 'standard_value_norm' generated from the function norm_value
        

    Returns
    -------
    x : Array pIC50 values of class float

    """

    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9)
        pIC50.append(-np.log10(molar))
    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', axis=1)

    return x

def norm_value(input):
    """Function to normalize all standard values greater than 1e6 to 1e6 as to not throw
    negative log values when calculating the IC50 value.

    Parameters
    ----------
    input : PD_dataframe
        Array with column 'standard_value' 

    Returns
    -------
    x : input array with standard_value_norm column of class float

    """

    norm = []

    for i in input['standard_value']:
        if i > 100000000:
            i = 100000000
        norm.append(i)
    input['standard_value_norm'] = norm
    x = input.drop('standard_value', axis=1)

    return x

# Target search
target = new_client.target

# query = input('Target query:')
query = 'Acetylcholinesterase'
target_query = target.search(query)
targets = pd.DataFrame.from_dict(target_query)
print(targets)

# selected_target = targets.target_chembl_id[int(input('Select Target Index:)'))]
selected_target = targets.target_chembl_id[0]
# print(selected_target)

# Select target protein
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50')
df = pd.DataFrame.from_dict(res)
# print(df)

current_directory = os.getcwd()
new_folder = os.path.join(current_directory, f'{query}_{selected_target}_drug_discovery')

# Create the folder
for i in range(11):
    try:
        if i == 0:
            os.makedirs(new_folder)
        else:
            new_folder_with_suffix = os.path.join(current_directory, f"{query}_{selected_target}_drug_discovery_{i:02}")
            os.makedirs(new_folder_with_suffix)
            new_folder = new_folder_with_suffix
        print(f"Folder '{new_folder}' created successfully.")
        break 
    except FileExistsError:
        if i == 10:
            print(f"All folder names '{query}_{selected_target}_drug_discovery' to '{query}_{selected_target}_drug_discovery_10' already exist.")
    except Exception as e:
        print(f"An error occurred: {e}")
        break

# Export original data to csv
df.to_csv(os.path.join(new_folder, f'{query}_01_Bioactivity_data.csv'), index=False)

# process data and get rid of missing values in dataset
df2 = pd.read_csv(os.path.join(new_folder, f'{query}_01_Bioactivity_data.csv'))
df2 = df2[df2.standard_value.notna()]
# print(df2.isna().sum())

# create a bioactivty class column
bioactivity_class = []
for i in df2.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append("inactive")
    elif float(i) <= 1000:
        bioactivity_class.append("active")
    else:
        bioactivity_class.append("intermediate")

# Make a compressed dataset with just the values we would like for our analysis
select = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df3 = df2[select].reset_index(drop=True)


df4 = pd.concat([df3, pd.DataFrame(bioactivity_class)], axis=1).dropna(subset=['canonical_smiles']).rename(columns = {0:'bioactivity_class'})
df4 = df4[df4['standard_value'] != 0].reset_index(drop=True)
# print(df4.isna().sum())

#Export preprocessed data to csv
df4.to_csv(os.path.join(new_folder, f'{query}_02_bioactivity_data_preprocessed.csv'), index=False)

df_lipinski = lipinski(df4.canonical_smiles).reset_index(drop=True)

# combine both data sets to make a complete dataset for analysis
df_combined = pd.concat([df4, df_lipinski], axis=1)
# print(df_combined)

# calculate the pIC50 value utilizing the functions above
df_norm = norm_value(df_combined)
df_final = pIC50(df_norm)

# Export data with calculated pIC50 
df_final.to_csv(os.path.join(new_folder,f'{query}_03_bioactivity_data_final_pIC50.csv'), index=False)

# Isolate active and inactive compounds for exploratory analysis
df_2_class = df_final[df_final['bioactivity_class'] != 'intermediate']

# Export data with active and inactive compounds excluding intermediates
df_2_class.to_csv(os.path.join(new_folder, f'{query}_04_bioactivity_data_2class.csv'), index=False)


def mannwhitney(descriptor, verbose=False):
    """Calculates mannwhitneyu statistics for each descriptor against bioactivity class
    for active and inactive compounds.

    Parameters
    ----------
    descriptor :
        
    verbose :
         (Default value = False)

    Returns
    -------

    """
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu

    seed(1)

    selection = [descriptor, 'bioactivity_class']
    df = df_2_class[selection]
    active = df[df['bioactivity_class'] == 'active']
    active = active[descriptor]

    selection = [descriptor, 'bioactivity_class']
    df = df_2_class[selection]
    inactive = df[df['bioactivity_class'] == 'inactive']
    inactive = inactive[descriptor]

    stat, p = mannwhitneyu(active, inactive)
    print('Statistics=%.3f, p=%.3f' % (stat, p))

    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distrubution (fail to reject H0)'
    else:
        interpretation = 'Different distrubution (reject H0)'
    
    results = pd.DataFrame({'Descriptor':descriptor, 
                            'Statistics':stat,
                            'p':p,
                            'alpha':alpha,
                            'Interpretation':interpretation}, index=[0])
    filename = f'mannwhitneyu_{descriptor}.csv'
    results.to_csv(os.path.join(new_folder, filename))
    return results

# Set graphical style for charts and graphs
sns.set_palette('deep')
sns.set_style("ticks")

# Frequency of bioactivity class (active vs inactive)
plt.figure(figsize=(5.5,5.5))
sns.countplot(x='bioactivity_class', data=df_2_class, edgecolor='black', hue='bioactivity_class')
plt.xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
plt.ylabel('Frequency', fontsize=14, fontweight='bold')
plt.show()
plt.savefig(os.path.join(new_folder, f'{query}_plot_bioactivity_class.pdf'))

# Molecular weight vs LogP of inactive and active compounds
plt.figure(figsize=(5.5,5.5))
sns.scatterplot(x='MW', y='LogP', size='pIC50', data=df_2_class, edgecolor='black', hue='bioactivity_class', alpha=0.7)
plt.xlabel('Molecular Weight', fontsize=14, fontweight='bold')
plt.ylabel('LogP', fontsize=14, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
plt.show()
plt.savefig(os.path.join(new_folder, f'{query}_plot_MW_vs_LogP.pdf'))

p_value = mannwhitney('MW')['p']
print(p_value)

# Muliplot analysis of lipinski descriptors 
fig, axs = plt.subplots(2,3, figsize=(18,10))
sns.boxplot(x='bioactivity_class', y='pIC50', data=df_2_class, hue='bioactivity_class', ax= axs[0,0])
axs[0,0].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
axs[0,0].set_ylabel('pIC50', fontsize=14, fontweight='bold')
# axs[0,0].text(0.05, 0.95, f'MSE = {p_value:.2f}', 
#             transform=axs[0].transAxes, fontsize=12, verticalalignment='top',
#             bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

sns.boxplot(x='bioactivity_class', y='MW', data=df_2_class, hue='bioactivity_class', ax=axs[0,1])
axs[0,1].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
axs[0,1].set_ylabel('MW', fontsize=14, fontweight='bold')

sns.boxplot(x='bioactivity_class', y='LogP', data=df_2_class, hue='bioactivity_class', ax=axs[1,0])
axs[1,0].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
axs[1,0].set_ylabel('LogP', fontsize=14, fontweight='bold')

sns.boxplot(x='bioactivity_class', y='NumAcceptors', data=df_2_class, hue='bioactivity_class', ax=axs[1,1])
axs[1,1].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
axs[1,1].set_ylabel('Number of H+ Acceptors', fontsize=14, fontweight='bold')

sns.boxplot(x='bioactivity_class', y='NumDonors', data=df_2_class, hue='bioactivity_class', ax=axs[0,2])
axs[0,2].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
axs[0,2].set_ylabel('Number of H+ Donors', fontsize=14, fontweight='bold')

sns.scatterplot(x='MW', y='LogP', size='pIC50', data=df_2_class, edgecolor='black', hue='bioactivity_class', alpha=0.5, ax=axs[1,2])
axs[1,2].set_xlabel('Molecular Weight', fontsize=14, fontweight='bold')
axs[1,2].set_ylabel('LogP', fontsize=14, fontweight='bold')
axs[1,2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)



###This Section is the ML Model. Need to use the bash file here to run the PaDel descriptors
selection = ['canonical_smiles', 'molecule_chembl_id']
df_selection = df_final[selection]
df_selection.to_csv(os.path.join(new_folder, f'{query}_{selected_target}_molecule.smi'), sep='\t', \
                    index=False, header=False)

# Export path directory to .txt document to be read in by our BASH file which will run the PaDEL program
with open(os.path.join('./PaDEL-Descriptor', 'output_dir.txt'), 'w') as f:
    f.write(new_folder)

#Run BASH file here############

# Import the generated fingerprints for all compounds in our data
df_descriptors = pd.read_csv(os.path.join(new_folder, 'descriptors_output.csv'))
# print(df_descriptors.columns)

# Train the Model

df_X = df_descriptors.drop(columns='Name')
print(df_X.shape)
df_y = df_final['pIC50']
print(df_y.shape)

random = 42

# Review this portion for understanding removal of low variance features
selection = VarianceThreshold(threshold=(.8 * (1-.8)))
X = selection.fit_transform(df_X)
# print(df_y.shape, X.shape)

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, df_y, test_size=0.2, random_state=random)
# print(X_train.shape, y_train.shape, X_test.shape, y_test.shape)

#Establish model and train it
model = RandomForestRegressor(n_estimators=100, random_state=random)
model_train = model.fit(X_train, y_train)
Y_pred = model.predict(X_test)
# print(Y_pred.shape)

# Lets establish a different model that performed better in our model analysis:
model2 = GradientBoostingRegressor(random_state=random, n_estimators=100)
model2.fit(X_train, y_train)
y2_pred = model2.predict(X_test)

MSE2 = mean_squared_error(y_test, y2_pred)
r22 = r2_score(y_test, y2_pred)
print(f"MSE: {MSE2:.2f}, R2: {r22:.2f}")
# How well did the model perform:
MSE = float(mean_squared_error(y_test, Y_pred))
r2 = float(r2_score(y_test, Y_pred))
print(f"MSE: {MSE:.2f}, R2: {r2:.2f}")



model3 = DecisionTreeRegressor(random_state=random)
model3.fit(X_train, y_train)
y3_pred = model3.predict(X_test)
MSE3 = mean_squared_error(y_test, y3_pred)
r23 = r2_score(y_test, y3_pred)
print(f"MSE: {MSE3:.2f}, R2: {r23:.2f}")

# Create visual of true values vs predicted values to visualize linear coorilation. 
sns.set(color_codes=True)
sns.set_style('white')

fig, axs = plt.subplots(1,3,figsize=(15,5))
sns.regplot(x=y_test, y=Y_pred,ax=axs[0], scatter_kws={'alpha':0.4, 'color':'orange'})
axs[0].set_title('Random Forest', fontsize=14, fontweight='bold')
axs[0].set_xlabel('Experiemental pIC50', fontsize=14, fontweight='bold')
axs[0].set_ylabel('True pIC50', fontsize=14, fontweight='bold')
axs[0].set_xlim(0, 12)
axs[0].set_ylim(0, 12)
axs[0].text(0.05, 0.95, f'MSE = {MSE:.2f}\nR2 = {r2:.2f}', 
            transform=axs[0].transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

sns.regplot(x=y_test, y=y2_pred, ax=axs[1], scatter_kws={'alpha':0.4, 'color':'green'})
axs[1].set_title('Gradient Boosting Regressor', fontsize=14, fontweight='bold')
axs[1].set_xlabel('Experiemental pIC50', fontsize=14, fontweight='bold')
axs[1].set_ylabel('True pIC50', fontsize=14, fontweight='bold')
axs[1].set_xlim(0, 12)
axs[1].set_ylim(0, 12)
axs[1].text(0.05, 0.95, f'MSE = {MSE2:.2f}\nR2 = {r22:.2f}', 
            transform=axs[1].transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

sns.regplot(x=y_test, y=y3_pred,ax=axs[2], scatter_kws={'alpha':0.4, 'color':'orange'})
axs[2].set_title('Decision Tree Regressor', fontsize=14, fontweight='bold')
axs[2].set_xlabel('Experiemental pIC50', fontsize=14, fontweight='bold')
axs[2].set_ylabel('True pIC50', fontsize=14, fontweight='bold')
axs[2].set_xlim(0, 12)
axs[2].set_ylim(0, 12)
axs[2].text(0.05, 0.95, f'MSE = {MSE3:.2f}\nR2 = {r23:.2f}', 
            transform=axs[2].transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))


clf = LazyRegressor(verbose=0, ignore_warnings=True, custom_metric=None)
train, test = clf.fit(X_train, X_test, y_train, y_test)
print(train)
print(test)

# Visualize the data

# train['R-Squared'] = [0 if i < 0 else i for i in train.iloc[:,0] ]

# print(train['R-Squared'])

sns.set(style="whitegrid")

fig, axs = plt.subplots(1, 3, figsize=(20,10))
fig.suptitle('Model Performance Metrics', fontsize=20, fontweight='bold')


axs[0] = sns.barplot(y=train.index, x="R-Squared", data=train, ax=axs[0])
axs[0].set(xlim=(0,1))
axs[0].set_xlabel("R-Squared", fontsize=15, fontweight='bold')
axs[0].set_ylabel("Model", fontsize=15, fontweight='bold')

axs[1] = sns.barplot(y=train.index, x="Time Taken", data=train, ax=axs[1])
axs[1].set(xlim=(0,10))
axs[1].set_xlabel("Time Taken", fontsize=15, fontweight='bold')
axs[1].set_ylabel("Model", fontsize=15, fontweight='bold')

axs[2] = sns.barplot(y=train.index, x="RMSE", data=train, ax=axs[2])
axs[2].set(xlim=(0,10))
axs[2].set_xlabel("RMSE", fontsize=15, fontweight='bold')
axs[2].set_ylabel("Model", fontsize=15, fontweight='bold')
plt.tight_layout()