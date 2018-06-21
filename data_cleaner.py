+import rdkit
+from rdkit import Chem
+import json
+from collections import Counter
+from rdkit.Chem import AllChem
+from rdkit.Chem import Draw
+import pandas as pd
+
+c = open('t_tmp')
+js = c.read()
+reaction_dicts = json.loads(js)
+reaction_df = pd.DataFrame(reaction_dicts)
+new_df = pd.DataFrame()
+
+new_df['ExperimentalYield'] = reaction_df['ExperimentalYield']
+new_df['Reagents'] = reaction_df.loc[:, "Reagents"].apply(lambda l: [e["SmilesFormula"] for e in l])
+new_df['Products'] = reaction_df.loc[:, "Products"].apply(lambda l: [e["SmilesFormula"] for e in l])
+new_rct_series = new_df["Reagents"].apply(lambda l: ".".join(l))
+new_prd_series = new_df["Products"].apply(lambda l: ".".join(l))
+
+def canon(df):
+    df_1 = []    
+    for smiles in df:
+        try:
+            canon = Chem.CanonSmiles(smiles)
+            df_1.append(canon)
+        except:    
+            mol = Chem.MolFromSmiles(smiles, sanitize=False)
+            rea = AllChem.ReactionFromSmarts("[n:1]([O-:2]):[o:3]>>[N+:1]([O-:2])=[O:3]")
+            rea.Initialize()
+            mol = Chem.MolToSmiles(rea.RunReactants([mol])[0][0])
+            try:
+                canon = Chem.CanonSmiles(mol)
+                df_1.append(canon)
+            except:
+                if '[n+]' in mol:
+                     mol = mol.replace('[n+]', '[N+]') 
+                if '(:o)' in mol:
+                    mol = mol.replace('(:o)', '(=O)')
+                if 'o' in mol:
+                     mol = mol.replace('o', 'O')
+                canon = Chem.CanonSmiles(mol)
+                df_1.append(canon)
+    return df_1
+    
+new_df['ReagentsCanon'] = new_df['Reagents'].apply(canon)
+new_df['ProductsCanon'] = new_df['Products'].apply(canon)
+
+dict_prod = dict()
+reaction_dicts = []
+for row in new_df.iterrows():
+    for elem in row[1]['ProductsCanon']:
+        reaction_dict = {'ProductCanon':elem, 'ReagentsCanon': row[1]['ReagentsCanon'], 'ExperimentalYield': row[1]['ExperimentalYield']}
+        reaction_dicts.append(reaction_dict)
+
+reaction_df = pd.DataFrame(reaction_dicts)
+new_reag = reaction_df["ReagentsCanon"].apply(lambda f: ".".join(f))
+reaction_df.loc[:, "ReactionSmiles"] = new_reag + ">>" + reaction_df["ProductCanon"]
+del reaction_df['ProductCanon']
+del reaction_df['ReagentsCanon']
+
+def group_reactions(yields):
+
+        if any([elem > 0 for elem in yields]):
+            return 1
+        elif any([elem == 0 for elem in yields]):
+            return 0 
+        else:
+            return -1
+            
+
+reaction_df = reaction_df.groupby('ReactionSmiles').agg({"ExperimentalYield": group_reactions})
+reaction_df = reaction_df.rename(columns={'ExperimentalYield': 'Flag'})
+reaction_df
