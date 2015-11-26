

import csv
import numpy as np
import matplotlib.pyplot as plt
import os

class csv_dict_writer:
    def __init__(self, filename, fieldnames=None):
        
        # if fieldnames should be default
        if fieldnames == None:        
            self.fieldnames = ['nr', 'P', 'Pdb', 'h', 'a', 'a_occ', 'R', 'time', 'R_ref']
        else:
            self.fieldnames = fieldnames
        
        # open csvfile
        self._csvfiledir = 'csvfiles/'
        self._csvfilename = self._csvfiledir + filename
        self.csvfile = open(self._csvfilename, 'wr')
        
        # make DictWriter
        self.w = csv.DictWriter(self.csvfile, fieldnames=self.fieldnames, delimiter='\t')
        
        # write header first
        self.w.writeheader()
        
        # some needed variables
        self._FLOATTYPES = [type(float()), type(np.float()), type(np.float128()), type(np.float16()), type(np.float32()), type(np.float64())]
        self._float_precision = "{: .2f}"
    
    def close(self):
        self.csvfile.close()
        
        
    def write_row(self, dic):
        shape_dic = {}  # holds the length of the values in dic (will be updated)
        max_dic = {}    # holds number of rows for each value in dic (won't be updated)
        """ generate shape_dic and max_dic with shape and maximum needed row information"""
        for key in dic.iterkeys():
            shape_dic[key] = np.shape(dic[key])
            if shape_dic[key] == ():
                shape_dic[key] = 1 # number of row is 1
                max_dic[key] = 1
            else:
                shape_dic[key] = shape_dic[key][0]  # number of rows, maximal rows
                max_dic[key] = shape_dic[key]    # maximum number of rows
                
        """ generate print_dic with row print information """
        print_dic = {}
        print_new_row = True
        while print_new_row:
            for key in dic.iterkeys():
                if shape_dic[key] == 0: # if already written
                    print_dic[key] = ''
                else:
                    position = max_dic[key] - shape_dic[key]    # calc position with nr of rows - current number
                    shape_dic[key] -= 1
                    try:
                        print_dic[key] = dic[key][position]
                    except TypeError:
#                        print "TypeError"
                        print_dic[key] = dic[key]
                    except IndexError:
#                        print "IndexError"
                        print_dic[key] = dic[key]
                    
                    # adjust precision of floats
                    if type(print_dic[key]) in self._FLOATTYPES:
                        print_dic[key] = self._float_precision.format(print_dic[key])
            if max(shape_dic.values()) == 0:
                print_new_row = False
                
            self.w.writerow(print_dic)
            
    def write_parameters(self, par_dict):
        """
        Writes settings into csvfile in a comment line
        """        
        settings_str = "#iterations={its} pstart={pstart} pend={pend} pnum={pnum} L={L} Ku={Ku} time_scale={time_scale}"
        settings_str = settings_str.format(**par_dict) + "\n"
        self.csvfile.write(settings_str)

class csv_dict_reader:
    def __init__(self):
#        self.csvfiledir = 'csvfiles/'        
        
        # find available files and choose one file to read
        csvfile = self.choose_csv_file('csvfiles/')

        # open file to read
        self.file = open(csvfile, 'r')
        
        self.r = csv.DictReader(self.file, delimiter='\t')
            
    def get_csv_files_avail(self, directory):
        """ Returns available .csv and .txt file in directory """
        files = os.listdir(directory)
        print files
        csvfiles_avail = []
        for file_str in files:
            file_format = file_str.split('.', 1)[1]
            if file_format in ['csv', 'txt']:
                csvfiles_avail.append(file_str)
                
        return csvfiles_avail
                
    def choose_csv_file(self, directory):
        """ chooses one file out of files in directory """
        csvfiles_avail = self.get_csv_files_avail(directory)
        print csvfiles_avail
        if not len(csvfiles_avail) == 0:
            for ind, file_name in enumerate(csvfiles_avail):
                print ind, file_name
            
            file_choose = int(raw_input('Which file should be read? ({}-{}): '.format(0, ind)))
            filename = csvfiles_avail[file_choose]
            csvfile = directory + filename
            print csvfile

        return csvfile
        
    def split_comment(self, comment_str):
        """ Returns a dict with parameters """
        print comment_str
#        par_keys = ['iterations', 'pstart', 'pend', 'pnum', 'L', 'Ku', 'time_scale']
        # get the string without leading #
        comment_str = comment_str.split('#')[1]
        # split all whitespaces and get a list with 'key=value'-pairs
        comment_list = comment_str.split()
        
        dic = {}
        for element in comment_list:
            split_list = element.split('=')
            key = split_list[0]
            value = float(split_list[1])
            dic[key] = value
        return dic
        
    def plot_r_p(self):
        """ Plots R (rate) vs p (power) """
        _undef_values = ['', None, ' ']
        first_run = True
        for row in self.r:
            # if comment line with new parameters
            if row['nr'].find('#') == 0:
                if not first_run:
                    if parameter_dic['L'] in [4, 8, 16]:
                        self.plot(Pdb, R, parameter_dic)
                parameter_dic = self.split_comment(row['nr'])
                P = []
                R = []
                Pdb = []
                first_run = False
                continue

            # iterate over all lines and grab P and R values
            for key, value in row.iteritems():
                if value in _undef_values:
                    continue
                if key == 'P':
                    P.append(float(value))
                elif key == 'Pdb':
                    Pdb.append(float(value))
                elif key == 'R':
                    R.append(float(value))
        # plot last setting
        if parameter_dic['L'] in [4, 8, 16]:
            self.plot(Pdb, R, parameter_dic)
                    
                    
    def plot(self, Pdb, R, parameter_dic):         
        title_str = "{iterations} iterations\nPdB={pstart}->{pend} ({pnum} steps)"
        plt.plot(Pdb, R, label="L={L:0} Ku={Ku:0}".format(**parameter_dic))
        plt.title(title_str.format(**parameter_dic))
        plt.xlabel('Power (P / dB)')
        plt.ylabel('average computation rate (bits/channel)')
        plt.legend(loc='upper left')
        
    def read_row(self):
        for row in self.r:
            print row   

if __name__ == '__main__':
    reader = csv_dict_reader()
    reader.plot_r_p()
    
    
    
    
    
    
    
    
    






