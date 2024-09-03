import tkinter as tk
from tkinter import ttk, messagebox
import pandas as pd
from chembl_webresource_client.new_client import new_client
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image, ImageTk
import io
import os
import numpy as np
from rdkit.Chem import Lipinski, Descriptors
from rdkit import Chem

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


class ChEMBLApp:
    def __init__(self, root):
        self.root = root
        self.root.title("ChEMBL Data Explorer")
        self.root.geometry("1850x700") 
        self.create_widgets()

    def create_widgets(self):
        # Frame for Search and Results
        self.left_frame = tk.Frame(self.root, width=700)
        self.left_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, pady=10)

        # Frame for Plot
        self.right_frame = tk.Frame(self.root, width=2000)
        self.right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, padx=10, pady=10)

        # Search Query Entry
        self.query_label = tk.Label(self.left_frame, text="Enter Target Query:")
        self.query_label.pack(pady=10)

        self.query_entry = tk.Entry(self.left_frame, width=50)
        self.query_entry.pack(pady=5)

        self.search_button = tk.Button(self.left_frame, text="Search", command=self.search_targets)
        self.search_button.pack(pady=10)

        # Results Treeview
        self.results_label = tk.Label(self.left_frame, text="Search Results:")
        self.results_label.pack(pady=10)

        self.results_tree = ttk.Treeview(self.left_frame, columns=('pref_name', 'organism', 'target_chembl_id', 'target_type'), show='headings')
        self.results_tree.heading('pref_name', text='Pref Name')
        self.results_tree.heading('organism', text='Organism')
        self.results_tree.heading('target_chembl_id', text='Chembl ID')
        self.results_tree.heading('target_type', text='Target Type')
        self.results_tree.pack(pady=5, fill=tk.BOTH, expand=True)

        self.select_button = tk.Button(self.left_frame, text="Select Target", command=self.select_target)
        self.select_button.pack(pady=10)

        # Progress Bar
        self.progress = ttk.Progressbar(self.left_frame, orient="horizontal", length=300, mode="determinate")
        self.progress.pack(pady=20)

        # Plot Area
        self.plot_frame = tk.Frame(self.right_frame)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        self.figures = {}

    def search_targets(self):
        query = self.query_entry.get()
        if not query:
            messagebox.showerror("Input Error", "Please enter a target query.")
            return
        
        self.results_tree.delete(*self.results_tree.get_children())
        self.progress.start()
        
        # Fetch results from ChEMBL
        try:
            target = new_client.target
            target_query = target.search(query)
            self.targets_df = pd.DataFrame.from_dict(target_query)
            
            if all(col in self.targets_df.columns for col in ['pref_name', 'organism', 'target_chembl_id', 'target_type']):
                for idx, row in self.targets_df.iterrows():
                    self.results_tree.insert('', 'end', values=(row['pref_name'], row['organism'], row['target_chembl_id'], row['target_type']))
            else:
                messagebox.showerror("Data Error", "One or more expected columns are not present in the search results.")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to fetch data: {e}")
        
        self.progress.stop()

    def select_target(self):
        selected_item = self.results_tree.selection()
        if not selected_item:
            messagebox.showerror("Selection Error", "Please select a target from the list.")
            return
        
        selected_item_id = self.results_tree.item(selected_item)['values'][2]  
        
        self.progress.start()
        
        # Fetch activity data for the selected target
        try:
            activity = new_client.activity
            res = activity.filter(target_chembl_id=selected_item_id).filter(standard_type='IC50')
            self.df = pd.DataFrame.from_dict(res)
            
            self.process_data()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to fetch activity data: {e}")
        
        self.progress.stop()

    def process_data(self):
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        # Initialize the progress bar
        self.progress['value'] = 0
        self.progress['maximum'] = 100
        self.progress.start()

        try:
            new_folder = 'output'
            os.makedirs(new_folder, exist_ok=True)

            # Export original data to CSV
            df_path = os.path.join(new_folder, f'{self.query_entry.get()}_01_Bioactivity_data.csv')
            self.df.to_csv(df_path, index=False)

            # Update progress
            self.progress['value'] = 10
            self.root.update_idletasks()

            # Process data and get rid of missing values
            df2 = pd.read_csv(df_path)
            df2 = df2[df2.standard_value.notna()]

            # Create bioactivity class column
            bioactivity_class = []
            for i in df2.standard_value:
                if float(i) >= 10000:
                    bioactivity_class.append("inactive")
                elif float(i) <= 1000:
                    bioactivity_class.append("active")
                else:
                    bioactivity_class.append("intermediate")

            # Update progress
            self.progress['value'] = 20
            self.root.update_idletasks()

            # Make a compressed dataset with just the values for analysis
            select = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
            df3 = df2[select].reset_index(drop=True)
            df4 = pd.concat([df3, pd.DataFrame(bioactivity_class)], axis=1).dropna(subset=['canonical_smiles']).rename(columns={0: 'bioactivity_class'})
            df4 = df4[df4['standard_value'] != 0].reset_index(drop=True)

            # Update progress
            self.progress['value'] = 30
            self.root.update_idletasks()

            # Calculate Lipinski descriptors
            df_lipinski = lipinski(df4.canonical_smiles).reset_index(drop=True)
            
            # Combine datasets for analysis
            df_combined = pd.concat([df4, df_lipinski], axis=1)
            df_norm = norm_value(df_combined)
            df_final = pIC50(df_norm)

            # Update progress
            self.progress['value'] = 50
            self.root.update_idletasks()

            # Isolate active and inactive compounds for exploratory analysis
            df_2_class = df_final[df_final['bioactivity_class'] != 'intermediate']

            # Update progress
            self.progress['value'] = 60
            self.root.update_idletasks()

            # Plotting
            sns.set_palette('deep')
            sns.set_style("ticks")

            # Multi-plot analysis of Lipinski descriptors
            fig, axs = plt.subplots(2, 3, figsize=(10, 6))
            sns.boxplot(x='bioactivity_class', y='pIC50', data=df_2_class, hue='bioactivity_class', ax=axs[0, 0])
            axs[0, 0].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
            axs[0, 0].set_ylabel('pIC50', fontsize=14, fontweight='bold')

            sns.boxplot(x='bioactivity_class', y='MW', data=df_2_class, hue='bioactivity_class', ax=axs[0, 1])
            axs[0, 1].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
            axs[0, 1].set_ylabel('MW', fontsize=14, fontweight='bold')

            sns.boxplot(x='bioactivity_class', y='LogP', data=df_2_class, hue='bioactivity_class', ax=axs[1, 0])
            axs[1, 0].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
            axs[1, 0].set_ylabel('LogP', fontsize=14, fontweight='bold')

            sns.boxplot(x='bioactivity_class', y='NumAcceptors', data=df_2_class, hue='bioactivity_class', ax=axs[1, 1])
            axs[1, 1].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
            axs[1, 1].set_ylabel('Number of H+ Acceptors', fontsize=14, fontweight='bold')

            sns.boxplot(x='bioactivity_class', y='NumDonors', data=df_2_class, hue='bioactivity_class', ax=axs[0, 2])
            axs[0, 2].set_xlabel('Bioactivity Class', fontsize=14, fontweight='bold')
            axs[0, 2].set_ylabel('Number of H+ Donors', fontsize=14, fontweight='bold')

            sns.scatterplot(x='MW', y='LogP', size='pIC50', data=df_2_class, edgecolor='black', hue='bioactivity_class', alpha=0.5, ax=axs[1, 2])
            axs[1, 2].set_xlabel('Molecular Weight', fontsize=14, fontweight='bold')
            axs[1, 2].set_ylabel('LogP', fontsize=14, fontweight='bold')
            axs[1, 2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
            plt.tight_layout()

            buf = io.BytesIO()
            plt.savefig(buf, format='png')
            buf.seek(0)
            
            # Display the plot in Tkinter
            img = Image.open(buf)
            photo = ImageTk.PhotoImage(img)
            
            img_label = tk.Label(self.plot_frame, image=photo)
            img_label.image = photo 
            img_label.pack(fill=tk.BOTH, expand=True)
            
            plt.close()

            # Update progress
            self.progress['value'] = 100
            self.root.update_idletasks()

        except Exception as e:
            messagebox.showerror("Processing Error", f"An error occurred while processing data: {e}")

        finally:
            # Stop the progress bar
            self.progress.stop()


if __name__ == "__main__":
    root = tk.Tk()
    app = ChEMBLApp(root)
    root.mainloop()
