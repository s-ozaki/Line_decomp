#!/usr/bin/env pyhton

import warnings
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, least_squares
from scipy.linalg import svd
from scipy.optimize import OptimizeWarning
from matplotlib.widgets import Button
from astropy.io import fits

# Local modules
import dbio
import funcs
import fitstools as ft

db=dbio.fit_database()
func=funcs.funcs()

def model(x, param):
    i=0
    mod=0.0
    pp=[]
    vpnum=[]
    for k in range(db.component_num):
        pp.clear()
        vpnum1=[]
        for l in range(db.funcinfo[k][1]):
            if db.params[k][l][5] == 0:
                pp.append(param[i])
                vpnum1.append(i)
                i+=1
            elif db.params[k][l][5] == -1:
                pp.append(db.params[k][l][0])
                vpnum1.append(-1)
            else:
                pp.append(param[ vpnum[ db.params[k][l][5]-1 ] [l] ] * \
                          db.params[k][l][3])
                vpnum1.append(-1)
        vpnum.append(vpnum1)
        mod += getattr(func,db.funcinfo[k][0])(x, *pp)
    return(mod)

def calc_residuals(params, x, y):
    model_y = model(x, params)
    return model_y - y
    
def model_for_plot(x):
    mod=0.0
    for k in range(db.component_num):
        pp=[]
        for i in range(len(db.params[k])):
            pp.append(db.params[k][i][0])
        mod += getattr(func,db.funcinfo[k][0])(x, *pp)
    return(mod)


class Index(object):
    
    def accept(self, event):
        global key
        key='y'
        plt.close()
        return

    def redo(self, event):
        global key
        key = 'n'
        plt.close()
        return

def specplot(data, titleline, message_switch, message):
    # Plot the result
    callback = Index()

    fig = plt.figure() 
    ax1 = fig.add_axes((0.1 ,0.35, 0.75, 0.55))
    ax2 = fig.add_axes((0.1, 0.1, 0.75, 0.2), sharex=ax1)
    
    ax1.set_title(titleline)
    ax1.step(data[:,0],data[:,1],where='mid')
    ax1.plot(data[:,0],model_for_plot(data[:,0]))
    for k in range(db.component_num):
        pp=[]
        for i in range(len(db.params[k])):
            pp.append(db.params[k][i][0])
        ax1.plot(data[:,0], getattr(func, db.funcinfo[k][0])(data[:,0], *pp))
    ax1.grid()
    ax1.set_ylabel('Intensity')
    ax1.tick_params(labelbottom="off")
    
    ax2.plot(data[:,0],data[:,1] - model_for_plot(data[:,0]))
    ax2.grid()
    ax2.set_xlabel('Wavelength ($\AA$)')
    ax2.set_ylabel('Residual')

    if not message_switch:
        ax1.text(0.5,0.5,message, horizontalalignment='center', 
                 verticalalignment='center', transform=ax1.transAxes, 
                 color='red', weight='bold')    

    axredo = plt.axes([0.88, 0.8, 0.1, 0.05])
    bredo = Button(axredo, 'Redo')
    bredo.on_clicked(callback.redo)

    axaccept = plt.axes([0.88, 0.7, 0.1, 0.05])
    baccept = Button(axaccept, 'Accept')
    baccept.on_clicked(callback.accept)
    plt.show()
    return


def fitting(fnames, area, dbname, outdbfname, without_plot, sample_range, 
            without_fit, method, en, absolute_sigma=False):

    titleline = ''
    for fname in fnames:
        titleline += fname
        
    global key
    key = 'n'
    while key == 'n':
        # Reading the parameters
        db.dbread(dbname)
        #print(db.params)

        # Selecting the data within the sample rages
        
        if sample_range != "*":
            datatmp = GetData(fnames[0], area, en)
            if len(fnames) > 1:
                for fname in fnames[1:]:
                    datatmp3 = GetData(fname, area, en)
                    datatmp = np.vstack((datatmp, datatmp3))
            wrangetmp = sample_range.split(",")
            wrange = []
            datatmp2 = []
            for i in range(len(wrangetmp)):
                wrange.append(wrangetmp[i].split("-"))
            for i in range(len(wrangetmp)):   
                for j in range(len(datatmp[:,0])):
                    if datatmp[j][0] >= float(wrange[i][0]) and \
                       datatmp[j][0] <= float(wrange[i][1]):
                        datatmp2.append([datatmp[j][0],datatmp[j][1]])
            data = np.array(datatmp2)
        else:
            data = GetData(fnames[0], area, en)
            if len(fnames) > 1:
                for fname in fnames[1:]:
                    datatmp3 = GetData(fname, area, en)
                    data = np.vstack((data, datatmp3))
            
        #plt.plot(data[:,0], data[:,1])
        #plt.show()
        
        if without_fit == False:
            # Scaling the data
            scale = np.average(data[:,1])
            data[:,1]=data[:,1]/scale

            # Scaling the parameters if needed
            db.params = func.scale(scale, db.component_num, db.funcinfo,
                                   db.params)

            # Set the initial parameters (p0) and
            # boundaries (maxbound, minbound) for fitting
            p0=[]
            minbound=[]
            maxbound=[]
            for i in range(db.component_num):
                for j in range(db.funcinfo[i][1]):
                    if db.params[i][j][5] == 0:
                        p0.append(db.params[i][j][0])      
                        minbound.append(db.params[i][j][1])
                        maxbound.append(db.params[i][j][2])
            # Fitting
            if method == 'lm':
                res = least_squares(calc_residuals, p0,
                                    args=(data[:,0],data[:,1]),
                                    method = method, verbose=2,
                                    xtol=3.0e-16, ftol=3.0e-16,
                                    gtol=3.0e-16, max_nfev=1000)
            elif method == 'trf' or method == 'dogbox':
                res = least_squares(calc_residuals, p0,
                                    args=(data[:,0],data[:,1]),
                                    method = method, verbose=2,
                                    bounds=(minbound,maxbound),
                                    xtol=3.0e-16, ftol=3.0e-16,
                                    gtol=3.0e-16, max_nfev=1000)
            else:
                print('Method is lm, trf or dogbox.')
                break
            
            ####### Taken from curve_fit function in minpack.py
            #if not res.success:
            #     raise RuntimeError("Optimal parameters not found: " + res.message)

            cost = 2 * res.cost  # res.cost is half sum of squares!
            popt = res.x

            # Do Moore-Penrose inverse discarding zero singular values.
            _, s, VT = svd(res.jac, full_matrices=False)
            threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
            s = s[s > threshold]
            VT = VT[:s.size]
            pcov = np.dot(VT.T / s**2, VT)

            warn_cov = False
            if pcov is None:
                # indeterminate covariance
                pcov = np.zeros((len(popt), len(popt)), dtype=float)
                pcov.fill(np.inf)
                warn_cov = True
            elif not absolute_sigma:
                if data[:,1].size > len(p0):
                    s_sq = cost / (data[:,1].size - len(p0))
                    pcov = pcov * s_sq
                else:
                    pcov.fill(np.inf)
                    warn_cov = True

            if warn_cov:
                warnings.warn('Covariance of the parameters could not be estimated',
                              category=OptimizeWarning)
            ###### End of quate

            perr = np.sqrt(np.diag(pcov))        

            # Canceling scaling and Substituting the resultant parameters
            # to the parameter list
            k=0
            for i in range(db.component_num):
                for j in range(db.funcinfo[i][1]):
                    if db.params[i][j][5] == 0:
                        db.params[i][j][0] = popt[k]
                        db.params[i][j][4] = perr[k]
                        k+=1
                    elif db.params[i][j][5] == -1:
                        db.params[i][j][4] = 0.0
                    else:
                        db.params[i][j][0] = db.params[db.params[i][j][5]-1][j][0]*db.params[i][j][3]
                        db.params[i][j][4] = db.params[db.params[i][j][5]-1][j][4]*db.params[i][j][3]

            # Putting back before scaling
            scale = 1.0/scale
            db.params = func.scale(scale, db.component_num, db.funcinfo, db.params)

        else:
            #data=datatmp
            scale=1.0
            for i in range(db.component_num):
                for j in range(db.funcinfo[i][1]):
                    if db.params[i][j][5] != 0 and db.params[i][j][5] != -1:
                        db.params[i][j][0] = db.params[db.params[i][j][5]-1][j][0]*db.params[i][j][3]
                        db.params[i][j][4] = db.params[db.params[i][j][5]-1][j][4]*db.params[i][j][3]

        # Plot
        if without_plot == False:
            data[:,1]=data[:,1]/scale
            if without_fit == False:
                specplot(data, titleline, res.success, res.message)
            else:
                specplot(data, titleline, 1, '')
        else:
            key ='y'

    if without_fit == False:
        # Writing the results
        db.dbwrite(fname, area, sample_range, dbname, method, outdbfname,
                   res.message)

    return

def GetData(fname, area, en):
    # Creating a data array including wavelength and intensity.
    # fname: Input FITS file name
    # area: Spatial area (string like x1,x2,y1,y2)
    # en: Extension number of data used for ftting.

    hdl = fits.open(fname)
    lam = ft.GetLambdaArr(hdl[en].header)
    intensities = hdl[en].data

    # Creating array of 1D intensity data
    a = [int(s) for s in area.split(',')]
    # a is coodinate in DS9
    temp = np.sum(intensities[:, a[2]-1:a[3], a[0]-1:a[1]], axis=1)
    intensity = np.sum(temp, axis=1)

    # Stacking arrays of lambadas and 1D intensity data
    data = np.stack((lam, intensity)).T

    hdl.close()
    
    return data


if __name__ == '__main__' :
    import argparse
    
    # Analize argments
    parser=argparse.ArgumentParser()
    parser.add_argument('fnames', help='Input file names')
    parser.add_argument('area', help='Area')
    parser.add_argument('dbfname',
                        help='database file name for initial parameters')
    parser.add_argument('outfname', help='output file name')
    parser.add_argument('-without_plot',action='store_true',
                        help='Without result plot')
    parser.add_argument('-without_fit',action='store_true', help='Without fit')
    parser.add_argument('-sample_range',action='store',
                        help='Sample wavelength ranges',default="*")
    parser.add_argument('-method',action='store',
                        choices=['lm', 'trf', 'dogbox'],
                        help='Algorithm to perform minimization. Default:trf',
                        default='trf')
    parser.add_argument('-en', type=int, help='Extension number. Default:0')

    args = parser.parse_args()

    # Set the matplotlib figure saving directory to the current directory.
    plt.rcParams['savefig.directory'] = os.getcwd()
    
    print('')
    fitting(args.fnames.split(','), args.area, args.dbfname, args.outfname,
            args.without_plot,
            args.sample_range, args.without_fit, args.method, args.en)
    print('')

