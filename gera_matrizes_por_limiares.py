import numpy as np
import pandas as pd
from datetime import datetime


dateTimeObj = datetime.now()
print(dateTimeObj)

limite = 101

for j in range(1, limite, 10):
    limiar=161.01
    
    # read the large csv file with specified chunksize 
    df_chunk = pd.read_csv('flow_matrix.csv', header=None,sep=';', skiprows=0, encoding='utf-8', chunksize=1000)
    
    chunk_list=[]
    
    # Each chunk is in df format
    for chunk in df_chunk: 
#        print('um chunk')
    
        for i in range(0, len(chunk.columns)):
            chunk[i] = np.where((chunk[i] <= limiar ), 0, 1)
        chunk_list.append(chunk)
    
    data = pd.concat(chunk_list)
    
    #Simetrizando a matriz
    data = data +data.T
    
    data.to_csv('adj_matrix_roads'+str(limiar)+'.csv') 
    dateTimeObj = datetime.now()
    print(dateTimeObj) 
