import argparse
import sys
import os
import pandas as pd
import gc

def process_batch(v_batch, snp):
   
    v_batch['DPWhole'] = v_batch['info'].str.split(';').str[0]
    v_batch['QSWhole'] = v_batch['info'].str.split(';').str[2]
    v_batch['DP'] = v_batch['DPWhole'].str.split('=').str[1]
    v_batch['QS_All'] = v_batch['QSWhole'].str.split('=').str[1]
    v_batch['QS'] = v_batch['QS_All'].str.split(',').str[0]
    v_batch['QS'] = v_batch['QS'].astype('float')
    v_batch['DP'] = v_batch['DP'].astype('float')
   
    #formula to find the variance
    v_batch['%_variance'] = round((v_batch['QS']/ v_batch['DP']) * 100, 2)
   
    v_batch['pos'] = v_batch['pos'].astype('int')
    v_batch = v_batch.drop(columns=['qual','filter','info','tags','taginfo','DPWhole','QSWhole','QS_All'])
   
    # Merge the two DataFrames on the primary key column
    merged_df = pd.merge(v_batch, snp, on=['chr', 'pos'], how='outer', indicator=True)

    # Filter out the rows present in v_batch but not in snp
    rows_only_in_v_batch = merged_df[merged_df['_merge'] == 'left_only']
   
    rows_only_in_v_batch = rows_only_in_v_batch.drop(columns=['strand','_merge'])
    return rows_only_in_v_batch

def process_files(file_name):    
   
    #to collect data processed in batches
    output_dfs = []
   
    # process data in batches
    for chunk_v in pd.read_csv(file_name, chunksize=10000, header=None, sep='\t',
                                names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info',
                                       'tags', 'taginfo']):
        filtered_df = chunk_v.dropna(how='any')
        filtered_df = filtered_df.iloc[1:]

        processed_batch = process_batch(filtered_df, snp)
        output_dfs.append(processed_batch)
        gc.collect()

    final_output_df = pd.concat(output_dfs, ignore_index=True)
    print ("Done with file ",file_name)
    return final_output_df

#finding common data between files
def findCommon(f1,f2,f3):
    common_f1_f2 = pd.merge(f1, f2, how='inner', on=list(f1.columns))
    f1 = None
    f2 = None
    common_values = pd.merge(common_f1_f2, f3, how='inner', on=list(common_f1_f2.columns))
    f3 = None
    return common_values

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument('--m1', required=True, help='vcf file converted to csv containing the output of mpileup1')
    ap.add_argument('--m2', required=True, help='vcf file converted to csv containing the output of mpileup2')
    ap.add_argument('--m3', required=True, help='vcf file converted to csv containing the output of mpileup3')
   
    ap.add_argument('--nm1', required=True, help='vcf file converted to csv containing the output of mpileup nm1')
    ap.add_argument('--nm2', required=True, help='vcf file converted to csv containing the output of mpileup nm2')
    ap.add_argument('--nm3', required=True, help='vcf file converted to csv containing the output of mpileup nm3')

    ap.add_argument('--snp', required=True, help='bed file containing the known SNPs')
    ap.add_argument('--o', required=True, help='path for output file')
    args = ap.parse_args()
   
    #reading the snp file
    snp = pd.read_csv(args.snp, header=None, sep='\t', names=['chr', 'pos-1', 'pos', 'strand'])
    snp = snp.drop(['pos-1'], axis=1)
   
    print("processing Ms")
    #processing Ms
    m1 = process_files(args.m1)
    m2 = process_files(args.m2)
    m3 = process_files(args.m3)
   
    print("finding commons Ms")
    common_ms = findCommon(m1,m2,m3)
   
    #running garbage collection after processing m
    gc.collect()
   
    print("processing NMs")
    #processing NMs
    nm1 = process_files(args.nm1)
    nm2 = process_files(args.nm2)
    nm3 = process_files(args.nm3)
   
    print("finding commons NMs")
    common_nms = findCommon(nm1,nm2,nm3)
   
    #running garbage collection after processing nm
    gc.collect()
   
    # Merge the two DataFrames on the primary key column
    merged_both = pd.merge(common_ms, common_nms, on=['chr', 'pos'], how='outer', indicator=True)
   
    common_ms = None
    common_nms = None

    # Filter out the rows present in v_batch but not in snp
    rows_only_in_ms = merged_both[merged_both['_merge'] == 'left_only']
    merged_both = None
    gc.collect()
   
    rows_only_in_ms.to_csv(args.o)
   
