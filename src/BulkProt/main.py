# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:57:56 2024

@author: u03132tk
"""
from BulkProt.cli import read_args, check_args
from BulkProt.functions import BulkProt
import time

def main(arg_list: list[str] | None = None):
    print ('Reading args...')
    args, parser = read_args(arg_list)
    print ('Checking args...')
    check_args(args, parser) 
    print ('Running BulkProt...')
    start = time.time()
    BulkProt(args.csv_path, 
             args.fields,#.split(','), 
             args.organism_id,
             args.seed_only,
             args.excel_compatible)
    parser.exit(f'Successfully completed run in {int(time.time() - start)} '\
                'seconds - thank you for using BulkProt')
    
#if __name__ == '__main__':
#    main()