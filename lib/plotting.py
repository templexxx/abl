from mpmath import *
from pylab import plot
from pylab import scatter
from pylab import savefig
from pylab import axis
from pylab import xlabel
from pylab import ylabel
from pylab import title
from pylab import xlim
from pylab import semilogx
from pylab import semilogy
from pylab import loglog
from pylab import ylim
from pylab import gca
from pylab import legend
from pylab import subplot
from pylab import errorbar
from pylab import LogFormatterExponent
from pylab import clf
from pylab import ScalarFormatter
import matplotlib

class DataPoints:
    
    def __init__(self):
        self.data_points = {}
        self.minx = None
        self.miny = None
        self.maxx = None
        self.maxy = None
        
    def addDataSet(self, x_coords, y_coords, name, errorbars=None):
        if len(x_coords) != len(y_coords):
            print "Coordinate dimension does not match!"
            exit(1)
            
        if max(x_coords) > self.maxx or self.maxx is None:
            self.maxx = max(x_coords)
        if min(x_coords) < self.minx or self.minx is None:
            self.minx = min(x_coords)
        if max(y_coords) > self.maxy or self.maxy is None:
            self.maxy = max(y_coords)
        if min(y_coords) < self.miny or self.miny is None:
            self.miny = min(y_coords)
        
        self.data_points[name] = (x_coords, y_coords, errorbars)
        
    def incrAddDataSet(self, x_coord, y_coord, name, errorbar=None):
        
        if not self.data_points.has_key(name):
            self.data_points[name] = ([], [], [])
        
        if x_coord > self.maxx or self.maxx is None:
            self.maxx = x_coord
        if x_coord < self.minx or self.minx is None:
            self.minx = x_coord
        if y_coord > self.maxy or self.maxy is None:
            self.maxy = y_coord
        if y_coord < self.miny or self.miny is None:
            self.miny = y_coord
        
        self.data_points[name][0].append(x_coord)
        self.data_points[name][1].append(y_coord)
        self.data_points[name][2].append(errorbar) 
    
    def __str__(self):
        return self.data_points.__str__()


class Plotting:
    LINE = "line"
    LOGLOG = "log log"
    XLOG = "x log"
    YLOG = "y log"
    ERRORBARS = "error bars"
    
    def __init__(self, data, filename, title="Title", xlbl="X Label", ylbl="Y Label", type="line", legend_loc=0):
        self.data = data
        self.filename = filename
        self.title = title
        self.xlbl = xlbl
        self.ylbl = ylbl
        self.type = type
        self.legend_loc = legend_loc
        self.linestyles=['o','^','v','<','>','s','+','x','D','d','1','2','3','4','h','H','p','|','_']
        self.errmarkers=['+', ',', '.', '1', '2', '3', '4', '+', ',', '.', '1', '2', '3', '4']
        self.colors=['r', 'b', 'g', 'm', 'k', 'c', 'y', 'c', 'k', 'm']


    def plot(self):
        if self.type == self.LINE:
            self.do_line()
        elif self.type == self.LOGLOG:
            self.do_loglog()
        elif self.type == self.XLOG:
            self.do_xlog()
        elif self.type == self.YLOG:
            self.do_ylog()
        elif self.type == self.ERRORBARS:
            self.do_errorbars()

    def do_line(self):
        i=0
        keys = self.data.data_points.keys()
        keys.sort()
        for name in keys:
            plot(self.data.data_points[name][0], self.data.data_points[name][1], "%c-%c" % (self.linestyles[i], self.colors[i]), label=name)
            i+=1
        legend(loc=self.legend_loc)
        title(self.title)
        xlabel(self.xlbl)
        ylabel(self.ylbl)
        xlim(self.data.minx, self.data.maxx)
        ylim(self.data.miny, self.data.maxy)
        savefig(self.filename)
        
    def do_loglog(self):
        i=0
        keys = self.data.data_points.keys()
        keys.sort()
        for name in keys:
            loglog(self.data.data_points[name][0], self.data.data_points[name][1], "%c-%c" % (self.linestyles[i], self.colors[i]), label=name)
            i+=1
        legend(loc=self.legend_loc)
        title(self.title)
        xlabel(self.xlbl)
        ylabel(self.ylbl)
        xlim(self.data.minx, self.data.maxx)
        ylim(self.data.miny, self.data.maxy)
        savefig(self.filename)
        
    def do_xlog(self):
        i=0
        keys = self.data.data_points.keys()
        keys.sort()
        for name in keys:
            semilogx(self.data.data_points[name][0], self.data.data_points[name][1], "%c-%c" % (self.linestyles[i], self.colors[i]), label=name)
            i+=1
        legend(loc=self.legend_loc)
        title(self.title)
        xlabel(self.xlbl)
        ylabel(self.ylbl)
        xlim(self.data.minx, self.data.maxx)
        ylim(self.data.miny, self.data.maxy)
        savefig(self.filename)
    
    def do_ylog(self):
        i=0

        keys = self.data.data_points.keys()
        keys.sort()
        for name in keys:
            semilogy(self.data.data_points[name][0], self.data.data_points[name][1], "%c-%c" % (self.linestyles[i], self.colors[i]), label=name)
            i+=1
        legend(loc=self.legend_loc)
        title(self.title)
        xlabel(self.xlbl)
        ylabel(self.ylbl)
        xlim(self.data.minx, self.data.maxx)
        ylim(self.data.miny, self.data.maxy)
        savefig(self.filename)
        
    def do_errorbars(self):
        i=0

        ax = gca()
        ax.set_yscale("log")
        #ax.xaxis.set_major_formatter(matplotlib.ticker.LogFormatterMathtext(10))
        #ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(10))

        keys = self.data.data_points.keys()
        keys.sort()
        for name in keys:
            errorbar(self.data.data_points[name][0], self.data.data_points[name][1], yerr=self.data.data_points[name][2], marker="%c" % (self.errmarkers[i]), label=name, markersize=20.0)
            i+=1
            
        legend(loc=self.legend_loc)
        title(self.title)
        xlabel(self.xlbl)
        ylabel(self.ylbl)
        xlim(self.data.minx, self.data.maxx)
        ylim(self.data.miny, self.data.maxy)
        savefig(self.filename)
        clf()
    
    