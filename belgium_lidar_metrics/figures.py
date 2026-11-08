"""Shared figure style: a sequential blue ramp (dataviz reference palette, light = low), NoData in a neutral gray."""
from matplotlib.colors import LinearSegmentedColormap

BLUES = LinearSegmentedColormap.from_list("blues", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
BLUES.set_bad("#e1e0d9")
INK, MUTED, EDGE, BACKGROUND = "#0b0b0b", "#898781", "#c3c2b7", "#fcfcfb"
