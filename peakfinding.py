from scipy.signal import savgol_filter, find_peaks
from typing import Dict, List, Optional, Tuple
import numpy as np


def find_peaks_boundaries(y: np.ndarray,  # clipped transmission data for which to find peaks in
               full_data_size: int,  # size of the non-clipped y array, for width threshold
               # savgol_windows: List[int],  # savgol filtering windows, in samples [pre, post]
               # savgol_orders: List[int],  # savgol filtering orders [pre, post]
               peak_finder_width: float,  # peak finding width as fraction of full data length
               peak_finder_prominence: float  # peak finding prominence as span of clipped data
               ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:  # return filtered data and peak data
    """Filters data and carries out peak-finding routine, returning filtered data and peak indices and properties"""
    # filter data to enable differentiation, then differentiate data and filter again to prepare for peak-finding
    # initial_savgol: np.ndarray = savgol_filter(y, savgol_windows[0], savgol_orders[0])
    initial_savgol: np.ndarray = savgol_filter(y, 11, 3)
    indices, peak_data = find_peaks(  # submit prepared data to peak-finding routine
        initial_savgol,
        width=peak_finder_width * (full_data_size - 1),  # width as fraction of non-clipped data size
        prominence=peak_finder_prominence * (np.max(initial_savgol) - np.min(initial_savgol)))  # p as data span

    lefts, rights = [], []  # lists that will contain the peak boundaries
    first_order_diff = np.diff(initial_savgol)  # differentiate to second order to find turning points
    for peak in indices:  # iterate over found peaks
        for span, sign, limit, arr in zip(  # iterate over left/right associated properties
                [range(peak + 1, len(first_order_diff)), reversed(range(0, peak))],  # r: peak->end, l: 0->peak
                [-1, -1],                                    # r: incrementing idx, l: decrementing idc
                [len(first_order_diff) - 1, 0],  # r: limit to last idx, l: limit to 0
                [lefts, rights]):  # r: append to rights, l: append to lefts
            for idx in span:  # increment or decrement index from peak
                if (abs(idx - peak) > 100 or  # if search length is exceeded
                        (sign * first_order_diff[idx]) > 0 or  # or sign has changed as expected
                        idx == limit):  # or we have reached the limit of data range
                    arr.append(  # append boundary to list of boundaries
                        idx if abs(idx - peak) >= 10 else  # if the boundary is below the minimum size
                        (idx + (10 * sign)  # set boundary to minimum size instead
                         if 0 <= idx + (10 * sign) < len(first_order_diff)  # if that boundary exceeds limit
                         else limit))  # set boundary to limit instead
                    break  # stop iterating over range, a boundary has been found

    for label, data in zip(["lefts", "rights", "indices"],  # add boundaries and indices to peak data dictionary
                           [np.array(lefts), np.array(rights), indices]):  # convert boundary lists to numpy arrays
        peak_data[label] = data

    return initial_savgol, peak_data  # returned filtered data arrays and peak data dictionary


x, y = np.loadtxt("spring_dips_clipped_2.txt", dtype=float, skiprows=1, unpack=True)
data_x = x[2050:2300]
data_y = y[2050:2300]
data_size = len(data_x)
initial, peak_data = find_peaks_boundaries(data_y, data_size, 0.2, 0.2)
print(peak_data)
