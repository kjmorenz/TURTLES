import csv
import statistics

import matplotlib.pyplot as plt
import numpy

from .GN import GN
from .util import is_cross_correlation

class G4_T3(GN):
    def from_stream(self, stream_in):
        for c0, \
            c1, p10, p11, t10, t11, \
            c2, p20, p21, t20, t21, \
            c3, p30, p31, t30, t31, \
            counts \
            in stream_in:
            correlation = tuple(map(int, (c0, c1, c2, c3)))

            floats = lambda *x: tuple(map(float, x))
            p1 = floats(p10, p11)
            t1 = floats(t10, t11)
            p2 = floats(p20, p21)
            t2 = floats(t20, t21)
            p3 = floats(p30, p31)
            t3 = floats(t30, t31)
            
            counts = int(counts)

            if correlation not in self._counts:
                self[correlation] = dict()

            if p1 not in self._counts[correlation]:
                self._counts[correlation][p1] = dict()

            if t1 not in self._counts[correlation][p1]:
                self._counts[correlation][p1][t1] = dict()

            if p2 not in self._counts[correlation][p1][t1]:
                self._counts[correlation][p1][t1][p2] = dict()

            if t2 not in self._counts[correlation][p1][t1][p2]:
                self._counts[correlation][p1][t1][p2][t2] = dict()

            if p3 not in self._counts[correlation][p1][t1][p2][t2]:
                self._counts[correlation][p1][t1][p2][t2][p3] = dict()

            if t3 not in self._counts[correlation][p1][t1][p2][t2][p3]:
                self._counts[correlation][p1][t1][p2][t2][p3][t3] = dict()

            self._counts[correlation][p1][t1][p2][t2][p3][t3] = counts

    def to_stream(self):
        for correlation in self:
            gn = self[correlation]
            
            for p1 in sorted(gn):
                for t1 in sorted(gn[p1]):
                    for p2 in sorted(gn[p1][t1]):
                        for t2 in sorted(gn[p1][t1][p2]):
                            for p3 in sorted(gn[p1][t1][p2][t2]):
                                for t3 in sorted(gn[p1][t1][p2][t2][p3]):
                                    line = list(correlation[:2]) + \
                                           list(p1) + list(t1) + \
                                           [correlation[2]] + \
                                           list(p2) + list(t2) + \
                                           [correlation[3]] + \
                                           list(p3) + list(t3) + \
                                           [gn[p1][t1][p2][t2][p3][t3]]

                                    yield(line)

    def unique_peaks(self):
        """
        Report the center, diagonal, and off-diagonal values.
        """
        peaks = {"g4": 0,
                 "g3g1": 0,
                 "g2g2": 0,
                 "g1g1g1g1": 0}

        for correlation in self.cross_correlations():
            gn = self[correlation]
            for p1, p2, p3, peak in \
                [((-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5), "g4"),
                 ((-0.5, 0.5), (-0.5, 0.5), (0.5, 1.5), "g3g1"),
                 ((-0.5, 0.5), (0.5, 1.5), (0.5, 1.5), "g2g2"),
                 ((0.5, 1.5), (1.5, 2.5), (3.5, 4.5), "g1g1g1g1")]:
                for t1 in gn[p1]:
                    for t2 in gn[p1][t1][p2]:
                        peaks[peak] += sum(gn[p1][t1][p2][t2][p3].values())
            
        return(peaks)
        
    def make_figure(self):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        pulses, g3 = self.combine()

        colorbar_ax = ax.imshow(g3,
                                origin="lower",
                                interpolation="none",
                                extent=[pulses[0]-0.5, pulses[-1]+0.5,
                                        pulses[0]-0.5, pulses[-1]+0.5])
        ax.set_xlabel(r"$\Delta p_{1}$")
        ax.set_ylabel(r"$\Delta p_{2}$")
        
        fig.colorbar(colorbar_ax)

        return(fig)

    def combine(self):
        total_g3 = None
        pulses = None
        
        for correlation in self.cross_correlations():
            g3 = self[correlation]

            if pulses is None:
                pulses = list(map(statistics.mean, sorted(g3.keys())))
                n_pulses = len(pulses)
##                print(pulses)
                total_g3 = numpy.zeros((n_pulses, n_pulses))

            for i0, p0 in enumerate(sorted(g3)):
                for t0 in g3[p0]:
                    for i1, p1 in enumerate(sorted(g3[p0][t0])):
                        total_g3[i0,i1] += sum(g3[p0][t0][p1].values())

        return(pulses, total_g3)

class G3_T2(GN):
    pass
