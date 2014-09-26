#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Plot scaling and wavelet functions for db, sym, coif, bior and rbio families

import itertools
import numpy, sys, os

from matplotlib import pyplot as plt

import pywt

iterations = 5

plot_data = [
    ('db', (4, 3)),
    ('sym', (4, 3)),
    ('coif', (3, 2))
]

for family, (rows, cols) in plot_data:
    f = pyplot.figure()
    f.subplots_adjust(
        hspace=0.2, wspace=0.2, bottom=.02, left=.06, right=.97, top=.94
    )
    colors = itertools.cycle('bgrcmyk')

    wnames = pywt.wavelist(family)
    print wnames
    i = iter(wnames)
    for col in xrange(cols):
        for row in xrange(rows):
            try:
                wavelet = pywt.Wavelet(i.next())
            except StopIteration:
                break
            phi, psi, x = wavelet.wavefun(iterations)

            color = colors.next()
            ax = pyplot.subplot(rows, 2 * cols, 1 + 2 * (col + row * cols))
            pyplot.title(wavelet.name + " phi")
            pyplot.plot(x, phi, color)
            pyplot.xlim(min(x), max(x))

            ax = pyplot.subplot(rows, 2 * cols, 1 + 2 * (col + row * cols) + 1)
            pyplot.title(wavelet.name + " psi")
            pyplot.plot(x, psi, color)
            pyplot.xlim(min(x), max(x))

for family, (rows, cols) in [('bior', (4, 3)), ('rbio', (4, 3))]:
    f = pyplot.figure()
    f.subplots_adjust(hspace=0.5, wspace=0.2, bottom=.02, left=.06, right=.97,
        top=.94)

    colors = itertools.cycle('bgrcmyk')
    wnames = pywt.wavelist(family)
    i = iter(wnames)
    for col in xrange(cols):
        for row in xrange(rows):
            try:
                wavelet = pywt.Wavelet(i.next())
            except StopIteration:
                break
            phi, psi, phi_r, psi_r, x = wavelet.wavefun(iterations)
            row *= 2

            color = colors.next()
            ax = pyplot.subplot(2 * rows, 2 * cols, 1 + 2 * (col + row * cols))
            pyplot.title(wavelet.name + " phi")
            pyplot.plot(x, phi, color)
            pyplot.xlim(min(x), max(x))

            ax = pyplot.subplot(2 * rows, 2 * cols,
                1 + 2 * (col + row * cols) + 1)
            pyplot.title(wavelet.name + " psi")
            pyplot.plot(x, psi, color)
            pyplot.xlim(min(x), max(x))

            row += 1
            ax = pyplot.subplot(2 * rows, 2 * cols, 1 + 2 * (col + row * cols))
            pyplot.title(wavelet.name + " phi_r")
            pyplot.plot(x, phi_r, color)
            pyplot.xlim(min(x), max(x))

            ax = pyplot.subplot(2 * rows, 2 * cols,
                1 + 2 * (col + row * cols) + 1)
            pyplot.title(wavelet.name + " psi_r")
            pyplot.plot(x, psi_r, color)
            pyplot.xlim(min(x), max(x))

pyplot.show()