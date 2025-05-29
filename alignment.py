import pandas as pd
import numpy as np


def Create_Table(X,Y):
    df = pd.read_csv("gfg.txt", sep=" ")

#Aligns the two fasta files
def aligner(df,blo62):
    num_rows, num_cols = df.shape
    df.iloc[0,0]=1
    
    for y in range(num_cols):
        for x in range(num_rows):
            if y == 0 and x== 0:
                continue
            elif y == 0:
                df.iloc[x,y]= df.iloc[x-1,y]+blo62.loc[df.index[x],"*"]
                
            elif x == 0:

                df.iloc[x,y]= df.iloc[x,y-1]+blo62.loc["*",df.columns[y]]
            else:
                
                mach = df.iloc[x-1,y-1]+blo62.loc[df.index[x],df.columns[y]]
                skip_l = df.iloc[x-1,y]+blo62.loc[df.index[x],"*"]
                skip_u = df.iloc[x,y-1]+blo62.loc["*",df.columns[y]]
                df.iloc[x,y]=max(mach,skip_l,skip_u)
    
    print(df.iloc[num_rows-1,num_cols-1])
    return_alignment(df,blo62)
    print(df)
 
def return_alignment(df,blo62):
    string1=""
    string2=""
    row, col = df.shape
    row=row-1
    col=col-1
    while row != 0 and col != 0:
        if row == 0 and col != 0:
            
            col=col-1
            string1=string1+"_"
        elif col == 0 and row != 0:
            string2=string2+"_"
            row=row-1
        elif df.iloc[row,col] - blo62.loc[df.index[row],df.columns[col]] == df.iloc[row-1,col-1]:
            string1=string1+df.index[row]
            string2=string2+df.columns[col]
            row=row-1
            col=col-1

        elif df.iloc[row,col] - blo62.loc[df.index[row],"*"] == df.iloc[row-1,col]:

            string1 = string1 + df.index[row]
            string2=string2+"_"
            row = row-1
        else:
            string2 = string2 + df.columns[col]
            string1=string1+"_"
            col = col-1
    string1 = string1[::-1]
    string2 = string2[::-1]
    print(string1)
    print(string2)
if __name__ == "__main__":
    blo62_path = input("Enter the path to the BLOSUM62 matrix file: ").strip('"')
    seq1_path = input("Enter the path to the first FASTA file (seq1): ").strip('"')
    seq2_path = input("Enter the path to the second FASTA file (seq2): ").strip('"')
    blo62 = pd.read_table(blo62_path, delimiter="\t", index_col=0).astype('float64')
    
    
    with open(seq1_path, 'r') as file:
        file.readline()
        text = "".join(line.strip() for line in file)
    X = "*" + text

    # Read second sequence
    with open(seq2_path, 'r') as file:
        file.readline()
        text = "".join(line.strip() for line in file)
    Y = "*" + text
    
    df = pd.DataFrame(index=list(Y), columns=list(X))
    df = df.astype('float64')
    blo62 = blo62.astype('float64')
    aligner(df,blo62)


 
