# -*- coding: utf-8 -*-


import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

class csv_dict_writer:
    def __init__(self, filename, fieldnames=None, filedir='csvfiles'):
        
        # if fieldnames should be default
        if fieldnames == None:        
            self.fieldnames = ['nr', 'P', 'Pdb', 'h', 'a', 'a_occ', 'Rmean', 'Rstd', 'time', 'R_ref_av']
        else:
            self.fieldnames = fieldnames
        
        # open csvfile
        self._csvfiledir = filedir
        self._csvfilename = filename
        self.csvfile = open(self._csvfiledir + '/' + self._csvfilename, 'w')
        
        # make DictWriter
        self.w = csv.DictWriter(self.csvfile, fieldnames=self.fieldnames, delimiter='\t')
        
        # write header first
        self.w.writeheader()
        
        # some needed variables
        self._FLOATTYPES = [type(float()), type(np.float()), type(np.float128()), type(np.float16()), type(np.float32()), type(np.float64())]
        self._FLOAT_PRECISION = "{: .6f}"
    
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
                        print_dic[key] = self._FLOAT_PRECISION.format(print_dic[key])
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
    def __init__(self, plot_par_dict):
        self.par = plot_par_dict
#        self.csvfiledir = 'csvfiles/'        
        
        # find available files and choose one file to read
        self.csvfilestring, self.csvfilename = self.choose_csv_file('csvfiles/')

        # open file for reading
        self.file = open(self.csvfilestring, 'r')
        # make new DictReader from csv
        self.r = csv.DictReader(self.file, delimiter='\t')
        # run correct routine to handl normal simulation and time measurement
        if self.csvfilename.find("qpr_run") == 0:
            self._plot_qpr_run()
        elif self.csvfilename.find("qpr_timeit") == 0:
            self._plot_qpr_timeit()
        else:
            print "Filename {} doesn't match qpr_run or qpr_timeit.".format(self.csvfilename)

        
    def _get_csv_files_avail(self, directory):
        """ Returns available .csv and .txt file in directory """
        files = os.listdir(directory)
#        print files
        csvfiles_avail = []
        for file_str in files:
            file_format = file_str.split('.', 1)[1]
            if file_format in ['csv', 'txt']:
                csvfiles_avail.append(file_str)
                
        return csvfiles_avail
                
    def choose_csv_file(self, directory):
        """ chooses one file out of files in directory """
        csvfiles_avail = self._get_csv_files_avail(directory)
#        print csvfiles_avail
        if not len(csvfiles_avail) == 0:
            for ind, file_name in enumerate(csvfiles_avail):
                print ind, file_name
            
            file_choose = int(raw_input('Which file should be read? ({}-{}): '.format(0, ind)))
            filename = csvfiles_avail[file_choose]
            csvfile = directory + filename
#            print csvfile

        return csvfile, filename
        
    def _split_comment(self, comment_str):
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
        
    def _plot_qpr_run(self):
        """ Plots R (rate) vs p (power) """
        _undef_values = ['', None, ' ']
        first_run = True
        P = []
        Rmean = []
        Rstd = []
        Pdb = []        
        Lrange = self.par['Lrange']
        nr = 0
        self.fig = plt.figure(num=1)
        for row in self.r:
            if self._row_is_comment(row):
                parameter_dic = self._split_comment(row['nr'])
                pnum = int(parameter_dic['pnum'])
                
            else:
                L = parameter_dic['L']
                if L in Lrange:
                    for key, value in row.iteritems():
                        if value in _undef_values:
                            continue
                        if key == 'P':
                            P.append(float(value))
                        elif key == 'Pdb':
                            Pdb.append(float(value))
                        elif (key == 'Rmean') or (key == 'R'):
                            Rmean.append(float(value))
                        elif (key == 'Rstd'):
                            Rstd.append(float(value))
                        elif key == 'nr':
                            nr = int(value)
                            
                if nr == pnum-1:
                    nr = 0
                    #print np.random.standard_normal(1)
                    #print len(Pdb), len(Rmean)
#                    self.fig = plt.figure(num=L)
                    self._qpr_run_plotter(Pdb, Rmean, Rstd, parameter_dic)
                    P = []
                    Rmean = []
                    Rstd = []
                    Pdb = []
                
                    
    def _read_iteration(self):
        """ """
        iteration_dict = {}        
        for row in self.r:
            if _row_is_comment(row):
                pass
            else:
                pass
    
    def _row_is_comment(self, row):
        is_comment = False
        if row['nr'].find('#') == 0: # '#' at position 0, >= 0 is also possible
            is_comment = True
        return is_comment

    def _qpr_run_plotter(self, Pdb, R, Rstd, parameter_dic):       
        ax = plt.gca()
        title_str = "{}\n".format(self.csvfilename) + "{iterations} iterations\nPdB={pstart}->{pend} ({pnum} steps)"
#        plt.plot(Pdb, R,
#                 label="L={L:0} Ku={Ku:0}".format(**parameter_dic),
#                 figure=self.fig,
#                 )
        
        if not self.par['errorbar']:
            Rstd = None
        plt.errorbar(Pdb, R,
                     yerr=Rstd,
                     label="L={L:0} Ku={Ku:0}".format(**parameter_dic),
                     figure=self.fig,
                     )
        plt.title(title_str.format(**parameter_dic))
        plt.xlabel('Power (P / dB)')
        plt.ylabel('average computation rate (bits/channel)')
        plt.legend(loc='upper left')
        ax.grid(b=True, axis='both', color='#858484', linestyle='--', linewidth=0.5)
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        ax.xaxis.set_ticks(np.arange(xmin, xmax+0.1, 2))
        ytick = np.round((ymax-ymin)/10, decimals=1)
        if ytick == 0:
            ytick = 0.05
        ax.yaxis.set_ticks(np.arange(ymin-ytick, ymax+ytick, ytick))
#        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(('%0.1f')))        
        
    def _plot_qpr_timeit(self):
        P = []
        T = []
        for row in self.r:
            P.append(row['P']); T.append(row['t'])
        plt.plot(P, T)
                       
    def _read_row(self):
        for row in self.r:
            print row   
    
    def _check_csv_filename(self, string):
        prefix = self.csvfilename.split('_')
        ret = False
        if prefix[0] == string:
            ret = True
        return ret
        

if __name__ == '__main__':
    
    plot_par_dict={
        "Prange": [0, 20, 2],
        "Lrange": range(2, 17, 2), # these L values will be plotted
        "errorbar": True,
    }    
    
    reader = csv_dict_reader(plot_par_dict)

    
    
    
    
    
    
    
    
    






