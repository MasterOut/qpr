# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 12:11:08 2015

@author: johannes
"""
#
import numpy as np
import qpr_fundamentials as fun
import matplotlib as plt
#import timeit
#import time
#import csv
#==============================================================================
# From qpr_main_new
#==============================================================================
def print_its_res(its, **kwargs):
    """
    Print function for all informations
    """
    line_string = get_string_format()
    value_dict = {}
    
    for i in range(0, its):
        for key, value in kwargs.iteritems():
            if np.size(value) > 1:
                value_dict[key] = value[i]
            else:
                value_dict[key] = value
        print line_string.format(**value_dict)
        
def print_processed_res(**kwargs):
    """
    Print function for processed informations
    """
    line_string = get_string_format(a_occ=True)
    value_dict = {}
    FILE.write(line_string+"\n")
    for key in kwargs.keys():
        if key == "h":
            value_dict["h"] = "std"
        elif key == "a":
            value_dict["a"] = kwargs["a"][0][0]
            value_dict["a_occ"] = kwargs["a"][0][1]
        else:
            value_dict[key] = kwargs[key]
#            value_dict[key] = float(kwargs[key])
    line_string = line_string.format(**value_dict)
#    print(line_string)
    FILE.write(line_string)
    """ Additional lines for a and occurence of a """
    value_dict = {} # clear value_dict to print only a and a_occ
    line_string = "\n"
    for i in range(1, len(kwargs["a"])):
        line_string += 4*"\t"+"{a}\t{a_occ}" +"\n"
        value_dict["a"] = kwargs["a"][i][0]
        value_dict["a_occ"] = kwargs["a"][i][1]
        line_string = line_string.format(**value_dict)
#    print(line_string)
    FILE.write(line_string+"\n")
    
def get_string_format(a_occ=False):   
    if a_occ:
        column_list = ["{nr}", "{p:.4}", "{pdb:.4}", "{h}", "{a}", "{a_occ}", "{cr:.2}", "{time:.4}", "{its}", "{cr_ref:.2}"]
             
    else:
        column_list = ["{nr}", "{p:.2}", "{pdb:.1}", "{h}", "{a}", "{cr:.2}", "{cr_ref:.2}"]
    
    string = len(column_list)*"{}\t"
    line_string = string.format(*column_list)   
    return line_string
    
def print_par_dictict(par_dict):
    string = ""
    for key in par_dict:
        string += key + ":\t{" + key+  "}\n"
    print(string.format(**par_dict))
#==============================================================================
# 
#==============================================================================
    
    
    
class file_dumper:
    def __init__(self, filename, cmdl_print=True, debug_print=False):
        self.filename = filename
        self.__cmdl_print = cmdl_print
        self.__deb_print = debug_print
        self.f = open(filename, "w")
        
    def def_columns(self, column_list):
        """
        column_tuple_list is a list of tuples in the form: (column name, position)
        """
        self.__disassamble_column_list(column_list)
        
        self.__sort_columns(column_list, self.__col_pos, self.__col_names)
        
        if self.__deb_print:
            self.__print_columns()
        
    def print_row(self, **kwargs):
        self.__check_column(**kwargs)
        col_names = self.__col_names # local copy
        
        for index, col_name in enumerate(col_names):
            if col_name in kwargs:  # wenn diese spalte beschrieben werden soll
                
                col_names[index] = "{" +col_name+ "}\t"
#                print index, col_name, col_names
            else:
                col_names[index] = "\t"
        col_str = "".join(col_names).format(**kwargs)
        print col_str
#        self.__print_columns()
    
    
    def __check_multirow(self, **kwargs):
        multirow = (1, )
        for key in kwargs.keys():
            if len(kwargs[key]) > multirow[0]:
                multirow = (len(kwargs[key]), key)
        return multirow
            
            
    def __check_column(self, **kwargs):
        """
        Checks wheather a given keyword=argument tuple is already defined in columns
        """
        for key in kwargs.keys():
            if not key in self.__col_names:
                if self.__deb_print:
                    print "A new column called '{}' was added at position '{}'".format(key, self.__add_column(key))
                
                
    def __add_column(self, name):
        self.__col_names.append(name)
        col_max = max(self.__col_pos)
        self.__col_pos.append(col_max+1)
        self.col_tuples.append((name, col_max+1))
        return col_max+1
    
    def __print_columns(self):
        print "Names:\t", self.__col_names
        print "pos:\t", self.__col_pos
        print "sorted tuple:\t", self.col_tuples
    
    def __disassamble_column_list(self, column_list):
        col_names = []
        col_pos = []
               
        for col in column_list:
            name = col[0]
            pos = col[1]
            
            # Check for double entries
            if self.__deb_print:
                if name in col_names:
                    print "ERROR: Column name is already defined."
                if pos in col_pos:
                    print "ERROR: Column position has already been set."
                
            col_names.append(name)
            col_pos.append(pos)
        self.__col_names = col_names
        self.__col_pos = col_pos
        
    def __sort_columns(self, col_tuple, col_pos, col_names):
        """
        Sorts columns in increasing order by position
        """
        indexes = np.argsort(col_pos)
        pos = 0
        col_n = [None]*len(col_names)
        col_t = [None]*len(col_tuple)
        col_p = [None]*len(col_pos)
        for i in indexes:
            col_n[pos] = col_names[i]
            col_t[pos] = (col_tuple[i][0], col_pos[i])
            col_p[pos] = col_pos[i]
            pos +=1
        self.col_tuples = col_t
        self.__col_names = col_n
        self.__col_pos = col_p

class csv_dict_writer:
    def __init__(self, filename, fieldnames):
        self.csvfile = open(filename, 'w')
        self.w = csv.DictWriter(self.csvfile, fieldnames=fieldnames, delimiter='\t')
        self.w.writeheader()
            
    def write_row(self, dic):
        self.w.writerow(dic)
    
    def close_witer(self):
        self.w.cose()
        
if __name__ == '__main__':
    


