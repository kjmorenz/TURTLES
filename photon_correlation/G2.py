import csv
import fractions
import statistics
import matplotlib
import matplotlib.pyplot as plt

from .GN import GN
from .util import *

t3_center = (-0.5, 0.5)
t3_side = (0.5, 1.5)

class G2_T3(GN):
    @property
    def time_resolution(self):
        corr = next(iter(self._counts.keys()))
        return(numpy.diff(next(iter(self[corr][t3_center].keys())))[0])
               
    def from_stream(self, stream_in):
        self._counts = dict()

        for c0, c1, pulse_left, pulse_right, \
            time_left, time_right, counts in stream_in:
            correlation = (int(c0), int(c1))

            pulse_bin = (float(pulse_left), float(pulse_right))
            time_bin = (float(time_left), float(time_right))

            if int(counts) == float(counts):
                counts = int(counts)
            else:
                counts = float(counts)

            if not correlation in self._counts:
                self._counts[correlation] = dict()
                
            if not pulse_bin in self._counts[correlation]:
                self._counts[correlation][pulse_bin] = dict()

            self._counts[correlation][pulse_bin][time_bin] = counts

        return(self)

    def to_stream(self):
        for correlation in sorted(self):
            for pulse_bin in sorted(self[correlation]):
                for time_bin in sorted(self[correlation][pulse_bin]):
                    line = itertools.chain(
                        correlation,
                        pulse_bin,
                        time_bin,
                        [self[correlation][pulse_bin][time_bin]])
                    
                    yield(line)

    def pulse_bin_counts(self, correlation, pulse_bin):
        return(sum(self.counts()[correlation][pulse_bin].values()))
    
    def center_side_ratios(self, center=t3_center, side=t3_side):
        try:
            self._center_side_ratios
        except:
            self._center_side_ratios = dict()
                
            for correlation, g2 in self:
                center_counts = self.pulse_bin_counts(correlation, center)
                side_counts = self.pulse_bin_counts(correlation, side)

                self._center_side_ratios[correlation] = \
                                        (center_counts, side_counts)
                
        return(self._center_side_ratios)

    def center_side_ratio(self, center=t3_center, side=t3_side):
        """
        Return the center/side ratio formed by summing over
        all cross-correlations in the g2.
        """
        if self._center_side_ratio is None:
            center_total = 0
            side_total = 0
            
            for correlation, ratio in self.g2_ratios().items():               
                if is_cross_correlation(correlation):
                    center, side = ratio
                
                    center_total += center
                    side_total += side

            self._center_side_ratio = float(center_total)/side_total

        return(self._center_side_ratio)

    def autocorrelation(self):
        """
        Add together all of the counts from the cross-correlations to get the
        resulting approximate autocorrelation.
        """
        try:
            self._autocorrelation
        except:
            self._autocorrelation = None
            
        if self._autocorrelation is None:
            self._autocorrelation = dict()
            
            for correlation in self.cross_correlations():
                g2 = self[correlation]

                for pulse_bin in g2.keys():
                    if not pulse_bin in self._autocorrelation.keys():
                        self._autocorrelation[pulse_bin] = dict()

                    for time_bin, counts in g2[pulse_bin].items():
                        if time_bin in self._autocorrelation[pulse_bin].keys():
                            self._autocorrelation[pulse_bin][time_bin] += counts
                        else:
                            self._autocorrelation[pulse_bin][time_bin] = counts

        return(self._autocorrelation)            
    
##    def make_figure(self):
##        fig = plt.figure()
##        ax = fig.add_subplot(1, 1, 1)
##        self.add_to_axes(ax)
##        return(fig)

    def make_figure(self, xscale=0, xzoom = [-1,-1], yzoom = [-1,-1],fontsize = 12, normalize = 1, scale = 0, color = 'b'):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        g2 = self.autocorrelation()
        #print(g2)
        max_time = round(max(map(lambda x: x[0], g2[(-0.5, 0.5)]))*1e-3)

        for pulse_bin, color in [((-1.5, -0.5), "black"),
                                 ((0.5, 1.5), "black"),
                                 ((-0.5, 0.5), "red")]:
            my_g2 = g2[pulse_bin]

            times = list(map(lambda x: statistics.mean(x[0])*1e-3
                             +statistics.mean(pulse_bin)*max_time,
                             sorted(my_g2.items())))
            counts = list(map(lambda x: x[1],
                              sorted(my_g2.items())))
            
            if not normalize == 0:
                for i in range(len(counts)):
                    counts[i] = counts[i]/normalize
            ax.plot(times, counts, color=color)

        
        ax.set_xlim((-max_time*1.5, max_time*1.5))
        ax.set_xticks((-max_time, 0, max_time))
        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels[0] = "-1"
        labels[-1] = "1"
        ax.set_xticklabels(labels)
        ax.set_xlabel("Pulses elapsed")
        ax.set_ylabel(r"$g^{(2)}(\rho,\tau)$")


##        ax.set_ylabel("$g^{(2)}$")
##        ax.set_xlabel("Time/ns")
##        
        if xscale == 1:
            ax.set_xscale('log')

        if not xzoom == [-1,-1]:
            ax.set_xlim((xzoom[0],xzoom[1]))
        
        if not yzoom == [-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            
        #ax.set_xlim((-500, 500))
        #ax.set_ylim((0, 0.02))

        fig.tight_layout()
    
        return(fig)
    def make_figure2(self, xscale=0, xzoom = [-1,-1], yzoom = [-1,-1], fontsize = 12, normalize = 1, scale = 0, color = 'b'):
        fig = plt.figure()
        matplotlib.rcParams.update({'font.size': fontsize})
        ax = fig.add_subplot(1, 1, 1)

        g2 = sorted(self.autocorrelation().items())
        
        times = list(map(lambda x: statistics.mean(x[0])*1e-3, g2))
        counts = list(map(lambda x: x[1], g2))
        ax.plot(times, counts, color)
        ax.set_ylabel("$g^{(2)}$")
        ax.set_xlabel("Time/ns")
        
        if xscale == 1:
            ax.set_xscale('log')

        if xzoom == [-1,-1]:
            if len(times) > 0:
                ax.set_xlim((times[0], times[-1]))
            else:
                ax.set_xlim((-1,1))
        else:
            ax.set_xlim((xzoom[0],xzoom[1]))
        
        if not yzoom == [-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            
        #ax.set_xlim((-500, 500))
        #ax.set_ylim((0, 0.02))

        fig.tight_layout()
    
        return(fig)

    def add_to_axis(self, ax, color):
        g2 = self.autocorrelation()
        #print(g2)
        max_time = round(max(map(lambda x: x[0], g2[(-0.5, 0.5)]))*1e-3)

        for pulse_bin, color in [((-1.5, -0.5), "black"),
                                 ((0.5, 1.5), "black"),
                                 ((-0.5, 0.5), color)]:
            my_g2 = g2[pulse_bin]

            times = list(map(lambda x: statistics.mean(x[0])*1e-3
                             +statistics.mean(pulse_bin)*max_time,
                             sorted(my_g2.items())))
            counts = list(map(lambda x: x[1],
                              sorted(my_g2.items())))
            
            ax.plot(times, counts, color=color)

        ax.set_xlim((-max_time*1.5, max_time*1.5))
        ax.set_xticks((-max_time, 0, max_time))
        ax.set_xlabel("Time/ns")
        ax.set_ylabel(r"$g^{(2)}(\rho,\tau)$")
        

    def total_counts(self, pulse_bin):
        """
        Return all counts associated with all cross-correlations and the given
        pulse bin.
        """
        result = None
        
        for correlation, g2 in self:
            if not is_cross_correlation(correlation):
                continue

            if result is None:
                result = g2[pulse_bin]
            else:
                for time_bin, counts in g2[pulse_bin].items():
                    result[time_bin] += counts

        return(result)

    def unique_peaks(self):
        peaks = {"center": 0,
                 "side": 0}

        for correlation in self.cross_correlations():
            gn = self[correlation]

            for pulse_bin, peak in [((-0.5, 0.5), "center"),
                                    ((0.5, 1.5), "side")]:
                peaks[peak] += sum(gn[pulse_bin].values())
                    
        return(peaks)

    def to_time_resolution(self, resolution=None):
        if not resolution:
            return(self)

        binning = int(round(resolution / self.time_resolution))
        if binning < 2:
            return(self)
        else:
            return(self.rebin_time(n=binning))

    def rebin_time(self, n=2):
        result = G2_T3()
        
        for correlation in self:
            g2 = self[correlation]
            result[correlation] = dict()

            for pulse_bin in g2:
                time_bins, counts = zip(*sorted(g2[pulse_bin].items()))
                start, stop = zip(*time_bins)
                time_bins = zip(pc.util.smooth(start, n=n),
                                pc.util.smooth(stop, n=n))
                counts = pc.util.rebin(counts, n=n)
                
                result[correlation][pulse_bin] = {
                    time_bin: count for time_bin, count in \
                    zip(time_bins, counts)}

        return(result)
                                  
class G2_T2(GN):       
    def from_file(self, filename, int_counts=True):
        with open(filename) as stream_in:
            return(self.from_stream(csv.reader(stream_in),
                                    int_counts=int_counts))

    def from_stream(self, stream_in, int_counts=True):
        for c0, c1, time_left, time_right, counts in stream_in:
            correlation = (int(c0), int(c1))
            time_bin = (float(time_left), float(time_right))

            if int_counts:
                counts = float(counts)
            else:
                counts = float(counts)

            if not correlation in self._counts.keys():
                self._counts[correlation] = dict()

            self._counts[correlation][time_bin] = counts

        return(self)

    def autocorrelation(self):
        try:
            self._autocorrelation
        except:
            self._autocorrelation = None
            
        if self._autocorrelation is None:
            result = dict()
            
            for correlation in self:
                if not is_cross_correlation(correlation):
                    continue #changed

                for time_bin, counts in self[correlation].items():
                    if time_bin not in result.keys():
                        result[time_bin] = 0

                    result[time_bin] += counts

            self._autocorrelation = result
            
        return(self._autocorrelation)

    def to_resolution(self, resolution):
        result = G2_T2()

        for correlation, g2 in self:
            times = list(map(lambda x: x[0], sorted(g2.items())))
            counts = list(map(lambda x: x[1], sorted(g2.items())))

            current_resolution = int(times[0][1] - times[0][0])
            binning = resolution / current_resolution

            if binning != 1:
                times = list(zip(smooth(map(lambda x: x[0], times), n=binning),
                                 smooth(map(lambda x: x[1], times), n=binning)))
                             
                counts = rebin(counts, n=binning)

                my_g2 = dict()
                
                for time_bin, count in zip(times, counts):
                    my_g2[time_bin] = count
            else:
                my_g2 = g2
                
            result[correlation] = my_g2

        return(result)

    def make_figure(self, xscale=0, xzoom = [-1,-1], yzoom = [-1,-1], fontsize = 12, figsize = [-1,-1], normalize = 0, scale = 0, color = 'b', linewidth = 0.5, normlen = 50):
        fig = plt.figure()
        matplotlib.rcParams.update({'font.size': fontsize})
        matplotlib.rc('xtick', labelsize = fontsize)
        ax = fig.add_subplot(1, 1, 1)

        g2 = sorted(self.autocorrelation().items())
        
        times = list(map(lambda x: statistics.mean(x[0])*1e-3, g2))
        counts = list(map(lambda x: x[1], g2))
        '''print(times)
        print("\n")
        print(counts)
        crash = please'''

        if not normalize == 0 or not scale == 0:
            if normalize == 1:
                #print(counts)
                normalize = counts[0]
                a = 1
                for j in range(1,int(len(counts)/5)):
                    a = a + 1
                    normalize = normalize + counts[j]
                normalize = normalize/a
            if not normalize == 0:
                if normalize == 1:
                    normalize = 0
                    if normlen > len(counts)/6:
                        print("Err: normlen more than 1/6 counts")
                        normlen = int(len(counts)/6)
                    for i in range(normlen):
                        normalize = normalize + counts[i]
                    normalize = normalize/normlen
                for i in range(len(counts)):
                    counts[i] = (counts[i]-scale)/normalize

        ax.plot(times, counts, color, linewidth = linewidth)
        ax.set_ylabel("$g^{(2)}$")
        ax.set_xlabel("Time/ns")
        if fontsize > 12:
            ax.xaxis.set_tick_params(size = 7, width=2)
            ax.yaxis.set_tick_params(size = 7, width=2)
        if not figsize == [-1,-1]:
            fig.set_size_inches(figsize[0], figsize[1])

        
        if xscale == 1:
            ax.set_xscale('log')

        if xzoom == [-1,-1]:
            if len(times) > 0:
                ax.set_xlim((times[0], times[-1]))
            else:
                ax.set_xlim((-1,1))
        else:
            ax.set_xlim((xzoom[0],xzoom[1]))
        
        if not yzoom == [-1,-1]:
            ax.set_ylim((yzoom[0],yzoom[1]))
            
        #ax.set_xlim((-500, 500))
        #ax.set_ylim((0, 0.02))

        fig.tight_layout()
    
        return(fig)

if __name__ == "__main__":
    g2 = G2_T3(filename="/home/tsbischof/Documents/data/microscopy/analysis/triexciton/2014-09-04_oc2014-04-08/oc2014-04-08_1e-5_dot_009_250nW_000.ht3.g2.run/g2")
    g2.to_file("/home/tsbischof/tmp/blargh.g2")
