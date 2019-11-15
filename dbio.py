# funcinfo: 2D list of the function information 
#   [comp_num][0]:component name, [comp_num][1]: number of parameters
# params: 3d list of the parameters
#   [comp_num][param_num][0]: value 
#   [comp_num][param_num][1]: lower limit
#   [comp_num][param_num][2]: upper limit
#   [comp_num][param_num][3]: scale factor
#   [comp_num][param_num][4]: uncertainty  (not confirmed)
#   [comp_num][param_num][5]: fix:-1, optimize: 0, link: component number
# component_num: number of components
#

import datetime

class fit_database:
    def get_param(self, items):
        param = []
        param.append(float(items[0])) # value
        param.append(float(items[1])) # lower limit
        param.append(float(items[2])) # upper limit
        param.append(float(items[3])) # step or scale factor
        param.append(float(items[4])) # uncertainty  (not confirmed)
        param.append(int(items[5])) # fix:-1, optimize: 0, link: component number
        return(param)

    def dbread(self, fname):
        header_line_num=7
        self.component_num=0
        self.header=[]
        self.params=[]
        self.funcinfo=[]
        i = 0 # incremental number for lines of databese file
        j = 0 # flag: first line for each component?
        k = 0 # incremental number for parameters of each component
        f = open(fname, 'r')
        lines = f.readlines()
        f.close()
        for i in range(len(lines)):
            items = lines[i].split()
            if len(items) > 0 and items[0][0:1] != '#':
                if i < header_line_num:
                    self.header.append(lines[i])
                elif j == 0:
                    comment=list(items)
                    del(comment[0:2])
                    comment=' '.join(comment)
                    self.funcinfo.append([items[0],int(items[1]),comment])
                    #print(self.funcinfo)
                    j = 1
                    k = 0
                    self.component_num+=1
                    param=[]
                elif j == 1:
                    if k < self.funcinfo[len(self.funcinfo)-1][1]-1:
                        param.append(self.get_param(items))
                        k += 1
                    elif k == self.funcinfo[len(self.funcinfo)-1][1]-1:
                        param.append(self.get_param(items))
                        self.params.append(param)
                        j=0
                        
        if len(self.funcinfo) != self.component_num:
            print("ERROR in reading the database file, " + fname)
        
    def dbwrite(self, fname, area, sample_range, dbname, method,
                outdbfname, status):
        d = datetime.datetime.today()
        self.header = [d.strftime("%Y-%m-%d %H:%M:%S\n")]
        self.header.append('Input FITS file: ' + fname+'\n')
        self.header.append('Integrated area: ' + area + '\n')
        self.header.append('Wavelenght range: '+ sample_range+'\n')
        self.header.append('Input database file: '+ dbname+'\n')
        self.header.append('Fitting algorism: ' + method + '\n')
        self.header.append('Termination reson: '+status+'\n')

        f = open(outdbfname, "w")
        f.writelines(self.header)
        for i in range(self.component_num):
            f.write(self.funcinfo[i][0]+"\t"+str(self.funcinfo[i][1])+" "+self.funcinfo[i][2]+"\n")
            for j in range(self.funcinfo[i][1]):
                for k in range(5):
                    f.write("%13g" % self.params[i][j][k])
                f.write("%4d\n" % self.params[i][j][5])
        f.close()
