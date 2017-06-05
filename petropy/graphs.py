# -*- coding: utf-8 -*-
"""
Graphs is a simple log viewer using matplotlib to create tracks of log data.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import xml.etree.ElementTree as ET


def create_log_viewer(log, template_xml_path = None, template_defaults = None, top = None, height = None):
    """
    Creates Log Viewer

    Uses matplotlib to create a figure and axes to display log data. An XML
    template is required to specify the tracks and curve colors. Default size
    is for 8.5" x 11" paper.

    Parameters
    ----------
    template_xml_path : str (default None)
        path to xml template.
    template_defaults : str {'raw', 'multimin_oil', 'multimin_gas', 'multimin_oil_sum', 'multimin_gas_sum'} (default None)
        name of default template options. Uses prebuilt template to display data
    top : float (default None)
        setting to set initial top of graph
    height : float (default None)
        setting to set amout of feet to display

    Returns
    -------
    fig : Figure
        matplotlib Figure object
    axes : ndarray
        numpy ndarray of AxesSubplot objects where each axes is a curve track

    Raises
    ------
    ValueError
        If template_defaults parameter is not in key of dictionary items

    ValueError
        If template_xml_path are both specified

    ValueError
        tick spacing must be specified in depth track

    ValueError
        number spacing must be specified in depth track

    ValueError
        left and right values for each curve must be specified in xml template

    ValueError
        curve must be in log

    ValueError
        curve display_name must be specified in xml template

    ValueError
        curve curve_name must be specified in xml template

    ValueError
        curve fill_color must be specified for cumulative tracks

    Example
    -------
    >>> # create plot for 11 x 17 paper with default raw template and show
    >>> import petropy as ptr
    >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
    >>> fig, axes = ptr.create_log_viewer(log, template_defaults = 'raw')
    >>> fig.set_size_inches(17, 11)
    >>> plt.show()

    >>> # create plot with custom template and save
    >>> import petropy as ptr
    >>> log = ptr.log_data('WFMP') # sample Wolfcamp log
    >>> t = 'path/to/my/custom/template.xml'
    >>> fig, axes = ptr.create_log_viewer(log, template_xml_path = t)
    >>> s = 'path/to/save/figure.png'
    >>> fig.savefig(s)

    """

    default_templates_paths = {
        'raw': 'default_raw_template.xml',
        'multimin_oil': 'default_multimin_oil_template.xml',
        'multimin_oil_sum': 'default_multimin_oil_sum_template.xml',
        'multimin_gas': 'default_multimin_gas_template.xml',
    }

    file_dir = os.path.dirname(__file__)
    if template_xml_path is None and template_defaults is None:
        template_xml_path = os.path.join(file_dir, '..', 'data', 'sample', 'default_raw_template.xml')

    elif template_xml_path is None and template_defaults is not None:
        if template_defaults in default_templates_paths:
            file_name = default_templates_paths[template_defaults]
        else:
            print('template_defaults paramter must be in:')
            for key in default_templates_paths:
                print(key)
            raise ValueError("%s is not valid template_defaults parameter" % template_defaults)
        template_xml_path = os.path.join(file_dir, '..', 'data', 'sample', file_name)

    elif template_xml_path is not None and template_defaults is None:
        template_xml_path = template_xml_path
    else:
        raise ValueError('template_xml_path and template_defaults cannot both be specified.')

    with open(template_xml_path, 'r') as f:
        root = ET.fromstring(f.read())

    if 'top' in root.attrib and top is None:
        top = float(root.attrib['top'])
    elif top is None:
        top = 1000

    if 'height' in root.attrib and height is None:
        bottom = top + float(root.attrib['height'])
    elif height is None:
        bottom = top + 500
    else:
        bottom = top + height

    tracks = root.findall('track')
    formations = root.findall('formation')
    num_tracks = len(tracks)

    fig, axes = plt.subplots(1, num_tracks + 2, sharey = True)

    ### add formation background color ###
    for formation in formations:

        if 'name' in formation.attrib:
            name = formation.attrib['name']
        else:
            raise ValueError('Formation name required in template %s' % template_xml_path)

        if 'color' in formation.attrib:
            color = formation.attrib['color']
        else:
            raise ValueError('Color required for formation %s in template %s' % (name, template_xml_path))

        if 'alpha' in formation.attrib:
            alpha = float(formation.attrib['alpha'])
        else:
            alpha = 0.5

        formation_top = log.tops[name]
        formation_bottom = log.next_formation_depth(name)
        for ax in axes:
            ax.axhspan(formation_top, formation_bottom, facecolor = color, alpha = alpha)

    numbers = None
    ticks = None
    track_widths = [0.5]
    track_height = 0.78

    depth_track_numbers = []
    track_names = []
    max_track_curves = 0
    c = 0
    for t, track in enumerate(tracks):

        ax = axes[t + 1]

        scale = 'linear'

        if 'display_name' in track.attrib:
            track_display_name = track.attrib['display_name']
            track_names.append(track_display_name)
        else:
            track_names.append(None)
            track_display_name = ''

        if 'scale' in track.attrib:
            scale = track.attrib['scale']
            ax.set_xscale(scale, nonposx = 'clip')

        track_width = 1
        if 'width' in track.attrib:
            track_width = float(track.attrib['width'])
        track_widths.append(track_width)


        if track_display_name == 'DEPTH':

            ax.set_xlim(0,1)
            depth_track_numbers.append(t + 1)

            if 'number_spacing' not in track.attrib:
                raise ValueError('number_spacing is required for depth track.')

            if 'tick_spacing' not in track.attrib:
                raise ValueError('tick_spacing is required for depth track.')

            font_size = 16
            if 'font_size' in track.attrib:
                font_size = float(track.attrib['font_size'])

            number_spacing = int(track.attrib['number_spacing'])
            tick_spacing = int(track.attrib['tick_spacing'])

            max_depth = log.curve_df.DEPTH.max()
            numbers = range(0, int(max_depth) + number_spacing, number_spacing)
            ticks = range(0, int(max_depth) + tick_spacing, tick_spacing)

            for n in numbers:
                ax.text(0.5, n, str(int(n)),
                        horizontalalignment='center',
                        verticalalignment = 'center',
                        clip_on = True,
                        fontsize = font_size)

        elif 'cumulative' in track.attrib:

            if 'left' not in track.attrib or 'right' not in track.attrib:
                raise ValueError('left and right values must be specified in cumulative tracks for template %s.' % template_xml_path)

            left = float(track.attrib['left'])
            right = float(track.attrib['right'])

            invert_axis = False
            if right < left:
                left, right = right, left
                invert_axis = True

            curve_sums = []
            colors = []
            summation = np.asarray(log.curve_df.DEPTH * 0)
            for c, curve in enumerate(track):
                ### names ###
                if 'curve_name' not in curve.attrib:
                    raise ValueError('Curve Name required in template at %s' % template_xml_path)
                curve_name = curve.attrib['curve_name']

                if curve_name not in log.curve_df.columns.values:
                    raise ValueError('Curve %s not found in log.' % curve_name)

                if 'fill_color' not in curve.attrib:
                    raise ValueError('Curve fill_color must be specificied for cumulative track in template at %s' % template_xml_path)

                fill_color = curve.attrib['fill_color']

                increase = summation + np.asarray(log.curve_df[curve_name])
                ax.fill_betweenx(log.curve_df.DEPTH,
                                 summation,
                                 increase,
                                 color = fill_color)

                summation = increase

                ### display ###
                if invert_axis:
                    x = (len(track) - float(c)) / len(track) - 1.0 / (2 * len(track))
                else:
                    x = float(c) / len(track) + 1.0 / (2 * len(track))

                ax.text(x, 1.03, curve_name, rotation = 90,
                        horizontalalignment = 'center',
                        verticalalignment = 'bottom',
                        transform = ax.transAxes,
                        fontsize = 12,
                        color = fill_color)


            if 'major_lines' in track.attrib:
                num_major_lines = int(track.attrib['major_lines']) + 1
                dist = abs(left - right) / num_major_lines
                major_lines = np.arange(left + dist, right, dist)
                for m in major_lines:
                    ax.plot((m, m), (0, log.curve_df.DEPTH.max()), color = '#c0c0c0', lw = 0.5)

            ax.set_xlim(left, right)
            if invert_axis:
                ax.invert_xaxis()

            c = int(c / 2) + 1
        else:

            for c, curve in enumerate(track):
                ### names ###
                if 'curve_name' not in curve.attrib:
                    raise ValueError('Curve Name required in template at %s' % template_xml_path)
                curve_name = curve.attrib['curve_name']

                if curve_name not in log.curve_df.columns.values:
                    raise ValueError('Curve %s not found in log.' % curve_name)

                ### style and scale ###

                left = None
                right = None
                if 'left' in curve.attrib:
                    left = float(curve.attrib['left'])
                    left_label = ' ' + curve.attrib['left']
                if 'right' in curve.attrib:
                    right = float(curve.attrib['right'])
                    right_label = curve.attrib['right'] + ' '

                if left is None:
                    print('Adjust template at %s.' % template_xml_path)
                    raise ValueError('Left X axis not found for curve %s.'\
                                     % curve_name)
                if right is None:
                    print('Adjust template at %s.' % template_xml_path)
                    raise ValueError('Right X axis not found for curve %s.'\
                                     % curve_name)

                line_style = '-'
                color = '#000000'
                width = 1
                alpha = 1
                marker = None
                marker_size = 0

                if 'line_style' in curve.attrib:
                    line_style = curve.attrib['line_style']
                if 'color' in curve.attrib:
                    color = curve.attrib['color']
                if 'width' in curve.attrib:
                    width = float(curve.attrib['width'])
                if 'alpha' in curve.attrib:
                    alpha = float(curve.attrib['alpha'])
                if 'marker' in curve.attrib:
                    marker = curve.attrib['marker']
                if 'marker_size' in curve.attrib:
                    marker_size = float(curve.attrib['marker_size'])


                if scale == 'log':
                    ax.set_xlim(left, right)
                    x = log.curve_df[curve_name]

                else:
                    ax.set_xlim(0,1)
                    m = (1 - 0) / (right - left)
                    b = -m * left
                    x = m * log.curve_df[curve_name] + b
                    left = 0
                    right = 1

                ### label ###
                if 'display_name' in curve.attrib:
                    ax.text(0.5, 0.98 + 0.035 * (c + 1),
                            curve.attrib['display_name'],
                            horizontalalignment = 'center',
                            verticalalignment = 'bottom',
                            transform = ax.transAxes,
                            fontsize = 12,
                            color = color)

                    ax.text(0, 0.98 + 0.035 * (c + 1),
                            left_label,
                            horizontalalignment = 'left',
                            verticalalignment = 'bottom',
                            transform = ax.transAxes,
                            fontsize = 12,
                            color = color)

                    ax.text(1, 0.98 + 0.035 * (c + 1),
                            right_label,
                            horizontalalignment = 'right',
                            verticalalignment = 'bottom',
                            transform = ax.transAxes,
                            fontsize = 12,
                            color = color)

                if 'fill' in curve.attrib:

                    fill_color = '#000000'
                    if 'fill_color' in curve.attrib:
                        fill_color = curve.attrib['fill_color']

                    if curve.attrib['fill'] == 'left':
                        ax.fill_betweenx(log.curve_df.DEPTH,
                                         left,
                                         x,
                                         color = fill_color)

                    elif curve.attrib['fill'] == 'right':
                         ax.fill_betweenx(log.curve_df.DEPTH,
                                          x,
                                          right,
                                          color = fill_color)

                if 'right_cutoff_fill' in curve.attrib:

                    fill_color = '#000000'
                    if 'right_cutoff_fill_color' in curve.attrib:
                        fill_color = curve.attrib['right_cutoff_fill_color']

                    cutoff_value = float(curve.attrib['right_cutoff_fill'])
                    if scale != 'log':
                        v = m * cutoff_value + b
                    else:
                        v = cutoff_value

                    ax.fill_betweenx(log.curve_df.DEPTH,
                                     v,
                                     x,
                                     color = fill_color,
                                     where = v < x)
                    ax.plot(v * np.ones(len(log.curve_df.DEPTH[v < x])),
                            log.curve_df.DEPTH[v < x], c = "#000000", lw = 0.5)

                if 'left_cutoff_fill' in curve.attrib:

                    fill_color = '#000000'
                    if 'left_cutoff_fill_color' in curve.attrib:
                        fill_color = curve.attrib['left_cutoff_fill_color']

                    cutoff_value = float(curve.attrib['left_cutoff_fill'])
                    if scale != 'log':
                        v = m * cutoff_value + b
                    else:
                        v = cutoff_value

                    ax.fill_betweenx(log.curve_df.DEPTH,
                                     v,
                                     x,
                                     color = fill_color,
                                     where = x < v)
                    ax.plot(v * np.ones(len(log.curve_df.DEPTH[x < v])),
                            log.curve_df.DEPTH[x < v], c = "#000000", lw = 0.5)

                if 'left_crossover' in curve.attrib:

                    fill_color = '#000000'
                    if 'left_crossover_fill_color' in curve.attrib:
                        fill_color = curve.attrib['left_crossover_fill_color']

                    left_curve = curve.attrib['left_crossover']
                    if left_curve not in log.curve_df.columns.values:
                        raise ValueError('Curve %s not found in log.' % left_curve)

                    if 'left_crossover_left' not in curve.attrib and \
                       'left_crossover_right' not in curve.attrib:
                       raise ValueError('left and right crossover values not found in template %s' \
                                        % template_xml_path)

                    left_crossover_left = float(curve.attrib['left_crossover_left'])
                    left_crossover_right = float(curve.attrib['left_crossover_right'])

                    if scale != 'log':
                        m = (1 - 0) / (left_crossover_right - left_crossover_left)
                        b = -m * left_crossover_left
                        v = m * log.curve_df[left_curve] + b
                    else:
                        v = log.curve_df[left_curve]

                    ax.fill_betweenx(log.curve_df.DEPTH,
                                     x,
                                     v,
                                     color = fill_color,
                                     where = v < x)

                if 'right_crossover' in curve.attrib:

                    fill_color = '#000000'
                    if 'right_crossover_fill_color' in curve.attrib:
                        fill_color = curve.attrib['right_crossover_fill_color']

                    left_curve = curve.attrib['right_crossover']
                    if left_curve not in log.curve_df.columns.values:
                        raise ValueError('Curve %s not found in log.' % left_curve)

                    if 'right_crossover_left' not in curve.attrib and \
                       'right_crossover_right' not in curve.attrib:
                       raise ValueError('left and right crossover values not found in template %s' \
                                        % template_xml_path)

                    left_crossover_left = float(curve.attrib['right_crossover_left'])
                    left_crossover_right = float(curve.attrib['right_crossover_right'])

                    if scale != 'log':
                        m = (1 - 0) / (left_crossover_right - left_crossover_left)
                        b = -m * left_crossover_left
                        v = m * log.curve_df[left_curve] + b
                    else:
                        v = log.curve_df[left_curve]

                    ax.fill_betweenx(log.curve_df.DEPTH,
                                     x,
                                     v,
                                     color = fill_color,
                                     where = v > x)

                ax.plot(x,
                        log.curve_df.DEPTH,
                        c = color,
                        lw = width,
                        ls = line_style,
                        marker = marker,
                        ms = marker_size)

            if scale == 'log':
                ax.xaxis.grid(True, which = 'both', color = '#e0e0e0')

            elif 'major_lines' in track.attrib:
                num_major_lines = int(track.attrib['major_lines']) + 1
                dist = abs(left - right) / num_major_lines
                major_lines = np.arange(left + dist, right, dist)
                for m in major_lines:
                    ax.plot((m, m), (0, log.curve_df.DEPTH.max()), color = '#c0c0c0', lw = 0.5)

        if c > max_track_curves:
            max_track_curves = c

    if max_track_curves < 5:
        max_track_curves = 5
    else:
        max_track_curves += 1

    ### adjust track widths ###
    track_widths.append(0.5)
    track_widths = np.asarray(track_widths)
    track_widths = track_widths / np.sum(track_widths)
    track_locations = [0]
    for t in range(1, len(track_widths)):
        track_locations.append(track_locations[t - 1] + track_widths[t - 1])

    for a, ax in enumerate(axes):
        post = ax.get_position()
        new_post = (track_locations[a], 0.01, track_widths[a], post.height)
        ax.set_position(new_post)
        if ticks is not None:
            ax.set_yticks(numbers, minor = False)
            ax.set_yticks(ticks, minor = True)
            ax.tick_params(axis = 'y', direction = 'inout', length = 6, width = 1, colors = '#000000', which = 'minor')
            if a not in depth_track_numbers:
                ax.yaxis.grid(True, which = 'major')
        else:
            ax.set_yticks([])
        ax.set_xticks([])

        if a > 0 and a < len(axes) - 1:
            track_title = track_names[a - 1]
            if track_title is not None:
                for i in range(max_track_curves):
                    track_title += '\n'
                ax.set_title(track_title, fontweight='bold')

    if top is not None and bottom is not None:
        plt.ylim((top, bottom))

    plt.gca().invert_yaxis()
    fig.set_size_inches(11, 8.5)

    return fig, axes
